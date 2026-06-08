"""
FAERS Pediatric Cardiac Safety Analysis — CardioSafe Pediatric
===============================================================
Downloads FDA FAERS quarterly ASCII files, filters to pediatric cases,
computes reporting odds ratios (ROR) for cardiac adverse events against
the 12-drug CardioSafe list, and aligns pharmacovigilance signal with
model-predicted delta-QTc.

Usage:
    python3 faers_pipeline.py --download          # fetch + cache FAERS zips
    python3 faers_pipeline.py --parse             # build unified parquet cache
    python3 faers_pipeline.py --analyze           # run ROR + model alignment
    python3 faers_pipeline.py --all               # full pipeline

Output (written to ../results/faers/):
    faers_drug_ror.csv              per-drug ROR with 95% CI
    faers_combo_ror.csv             pairwise combo ROR
    faers_model_alignment.csv       ROR vs model delta-QTc comparison

Data source:
    FDA FAERS ASCII quarterly files, 2015Q1–2024Q4
    https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html

Notes:
    - FAERS has known duplicate case issue; deduplication by primaryid
    - Pediatric filter: AGE_COD normalized to years, age < 18
    - ROR = (a/b) / (c/d); 95% CI from log-ROR ± 1.96*SE
    - Combo ROR uses cases where BOTH drugs appear in the same report
"""

import argparse
import io
import itertools
import re
import sys
import time
import warnings
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import requests
from scipy import stats
from tqdm import tqdm

warnings.filterwarnings("ignore")

# ── PATHS ──────────────────────────────────────────────────────────────────────
_SRC_DIR     = Path(__file__).resolve().parent
_ROOT_DIR    = _SRC_DIR.parent
_RESULTS_DIR = _ROOT_DIR / "results" / "faers"
_CACHE_DIR   = _ROOT_DIR / "data" / "faers_cache"
_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
_CACHE_DIR.mkdir(parents=True, exist_ok=True)

# ── QUARTERS ───────────────────────────────────────────────────────────────────
QUARTERS = [f"{yr}q{q}" for yr in range(2015, 2025) for q in range(1, 5)]
FAERS_BASE = "https://fis.fda.gov/content/Exports/faers_ascii_{quarter}.zip"

# ── DRUG ALIASES ───────────────────────────────────────────────────────────────
DRUG_ALIASES = {
    "Methylphenidate": ["methylphenidate", "ritalin", "concerta", "focalin",
                        "metadate", "quillivant", "daytrana"],
    "Amphetamine":     ["amphetamine", "adderall", "dextroamphetamine",
                        "lisdexamfetamine", "vyvanse", "dexedrine"],
    "Risperidone":     ["risperidone", "risperdal"],
    "Quetiapine":      ["quetiapine", "seroquel"],
    "Aripiprazole":    ["aripiprazole", "abilify"],
    "Sertraline":      ["sertraline", "zoloft"],
    "Fluoxetine":      ["fluoxetine", "prozac", "sarafem"],
    "Escitalopram":    ["escitalopram", "lexapro"],
    "Clonidine":       ["clonidine", "catapres", "kapvay"],
    "Guanfacine":      ["guanfacine", "intuniv", "tenex"],
    "Imipramine":      ["imipramine", "tofranil"],
    "Nortriptyline":   ["nortriptyline", "pamelor", "aventyl"],
}

# ── CARDIAC MedDRA PTs ─────────────────────────────────────────────────────────
CARDIAC_PTS = {
    "QT prolonged", "Electrocardiogram QT prolonged", "Torsade de pointes",
    "Ventricular tachycardia", "Ventricular fibrillation", "Cardiac arrest",
    "Sudden cardiac death", "Palpitations", "Heart rate increased",
    "Atrioventricular block", "Ventricular arrhythmia", "Long QT syndrome",
    "Electrocardiogram abnormal", "Syncope",
}
PRIMARY_CARDIAC_PTS = {
    "QT prolonged", "Electrocardiogram QT prolonged", "Torsade de pointes",
    "Ventricular tachycardia", "Ventricular fibrillation", "Cardiac arrest",
    "Sudden cardiac death", "Long QT syndrome", "Ventricular arrhythmia",
}

# ── MODEL PREDICTIONS (fallback if risk_grid_results.csv absent) ───────────────
MODEL_DQTC = {
    "Methylphenidate+Amphetamine":   23.8,
    "Nortriptyline+Methylphenidate": 21.7,
    "Methylphenidate+Aripiprazole":  17.9,
    "Amphetamine+Aripiprazole":      17.3,
    "Methylphenidate+Risperidone":   15.2,
    "Amphetamine+Risperidone":       14.8,
    "Methylphenidate+Quetiapine":    14.1,
    "Imipramine+Methylphenidate":    13.6,
    "Nortriptyline+Aripiprazole":    13.1,
    "Methylphenidate+Sertraline":    12.4,
    "Methylphenidate+Fluoxetine":    12.1,
    "Imipramine+Risperidone":        11.8,
    "Risperidone+Fluoxetine":         5.3,
    "Risperidone+Sertraline":         1.0,
    "Escitalopram+Clonidine":        -4.2,
    "Clonidine+Guanfacine":          -6.1,
}

# ── LOGGING HELPERS ────────────────────────────────────────────────────────────
def _ts():
    return datetime.now().strftime("%H:%M:%S")

def log(msg, level="INFO"):
    icons = {"INFO": "·", "OK": "✓", "WARN": "⚠", "ERR": "✗", "STEP": "▶"}
    print(f"[{_ts()}] {icons.get(level,'·')} {msg}", flush=True)

def log_section(title):
    bar = "─" * (70 - len(title) - 3)
    print(f"\n[{_ts()}] ══ {title} {bar}", flush=True)

def eta_str(elapsed_s, done, total):
    if done == 0:
        return "ETA: calculating..."
    rate = done / elapsed_s
    remaining = (total - done) / rate
    m, s = divmod(int(remaining), 60)
    return f"ETA: {m}m {s:02d}s"


# ─────────────────────────────────────────────────────────────────────────────
# 1. DOWNLOAD
# ─────────────────────────────────────────────────────────────────────────────

def download_faers(quarters=QUARTERS, force=False):
    log_section("DOWNLOAD")
    log(f"Target: {len(quarters)} quarters  →  {_CACHE_DIR}")
    log(f"Estimated download size: ~{len(quarters) * 50:.0f}–{len(quarters) * 80:.0f} MB total")
    log(f"Estimated time: {len(quarters) * 1.5 / 60:.0f}–{len(quarters) * 3 / 60:.0f} min "
        f"(FDA servers are slow; varies widely)")

    # Check what's already cached
    cached   = [q for q in quarters if (_CACHE_DIR / f"faers_ascii_{q}.zip").exists()]
    to_fetch = [q for q in quarters if q not in cached] if not force else quarters

    if cached and not force:
        log(f"Already cached: {len(cached)} quarters — skipping those", "OK")
    log(f"Will download: {len(to_fetch)} quarters")

    if not to_fetch:
        log("Nothing to download.", "OK")
        return []

    failed  = []
    t_start = time.time()

    pbar = tqdm(to_fetch, desc="Downloading FAERS", unit="quarter",
                bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]",
                file=sys.stdout)

    for i, q in enumerate(pbar, 1):
        dest = _CACHE_DIR / f"faers_ascii_{q}.zip"
        url  = FAERS_BASE.format(quarter=q)
        pbar.set_postfix_str(f"quarter={q}")

        try:
            r = requests.get(url, timeout=90, stream=True)
            r.raise_for_status()

            # Stream with inner byte-level progress bar
            total_bytes = int(r.headers.get("content-length", 0))
            chunk_size  = 1 << 20  # 1 MB

            with open(dest, "wb") as f, tqdm(
                total=total_bytes, unit="B", unit_scale=True,
                desc=f"  {q}", leave=False, file=sys.stdout,
            ) as byte_bar:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    f.write(chunk)
                    byte_bar.update(len(chunk))

            size_mb = dest.stat().st_size / 1e6
            elapsed = time.time() - t_start
            pbar.write(f"  [{_ts()}] ✓ {q}  ({size_mb:.1f} MB)  "
                       f"{eta_str(elapsed, i, len(to_fetch))}")

        except requests.HTTPError as e:
            pbar.write(f"  [{_ts()}] ✗ {q}: HTTP {e.response.status_code} — skipping")
            failed.append(q)
            if dest.exists():
                dest.unlink()
        except Exception as e:
            pbar.write(f"  [{_ts()}] ✗ {q}: {e} — skipping")
            failed.append(q)
            if dest.exists():
                dest.unlink()

    elapsed_total = time.time() - t_start
    m, s = divmod(int(elapsed_total), 60)
    log(f"Download complete in {m}m {s:02d}s", "OK")

    if failed:
        log(f"{len(failed)} quarters failed: {failed}", "WARN")
        log("Manual download: https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html", "WARN")
        log(f"Place zip files in: {_CACHE_DIR}", "WARN")
    else:
        log(f"All {len(to_fetch)} quarters downloaded successfully.", "OK")

    return failed


# ─────────────────────────────────────────────────────────────────────────────
# 2. PARSE
# ─────────────────────────────────────────────────────────────────────────────

def _read_table_from_zip(zf: zipfile.ZipFile, pattern: str) -> Optional[pd.DataFrame]:
    matches = [n for n in zf.namelist() if re.search(pattern, n, re.IGNORECASE)]
    if not matches:
        return None
    name = sorted(matches)[0]
    with zf.open(name) as f:
        raw = f.read().decode("latin-1", errors="replace")
    return pd.read_csv(io.StringIO(raw), sep="$", dtype=str,
                       low_memory=False, on_bad_lines="skip")


def _normalize_age_to_years(df: pd.DataFrame) -> pd.DataFrame:
    age = pd.to_numeric(df.get("age", df.get("AGE", pd.Series(dtype=str))),
                        errors="coerce")
    cod = df.get("age_cod", df.get("AGE_COD",
                 pd.Series(dtype=str, index=df.index)))
    cod  = cod.str.upper().fillna("YR")
    mult = cod.map({"YR": 1.0, "DY": 1/365.25, "WK": 7/365.25,
                    "MON": 1/12, "DEC": 10.0}).fillna(1.0)
    df = df.copy()
    df["age_years"] = age * mult
    return df


def parse_faers(quarters=QUARTERS):
    log_section("PARSE")

    zip_files = sorted(_CACHE_DIR.glob("faers_ascii_*.zip"))
    if not zip_files:
        log("No FAERS zip files found in cache. Run --download first.", "ERR")
        return None, None, None

    log(f"Found {len(zip_files)} zip files to parse")
    log("Each zip contains DEMO (demographics), DRUG, and REAC (reactions) tables")
    log("Estimated time: 10–20 min depending on disk speed")

    demo_parts, drug_parts, reac_parts = [], [], []
    parse_errors = []
    t_start = time.time()

    pbar = tqdm(zip_files, desc="Parsing quarters", unit="quarter",
                bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]",
                file=sys.stdout)

    total_cases = total_drug_rows = total_reac_rows = 0

    for zpath in pbar:
        quarter = zpath.stem.replace("faers_ascii_", "")
        pbar.set_postfix_str(f"quarter={quarter}")

        try:
            with zipfile.ZipFile(zpath) as zf:
                demo = _read_table_from_zip(zf, r"DEMO\w+\.txt")
                drug = _read_table_from_zip(zf, r"DRUG\w+\.txt")
                reac = _read_table_from_zip(zf, r"REAC\w+\.txt")

            if demo is None or drug is None or reac is None:
                missing = [t for t, d in [("DEMO",demo),("DRUG",drug),("REAC",reac)] if d is None]
                pbar.write(f"  [{_ts()}] ⚠ {quarter}: missing tables {missing} — skipping")
                parse_errors.append(quarter)
                continue

            demo.columns = demo.columns.str.lower()
            drug.columns = drug.columns.str.lower()
            reac.columns = reac.columns.str.lower()

            demo = _normalize_age_to_years(demo)
            demo["quarter"] = quarter

            demo_keep = [c for c in ["primaryid","age_years","sex","event_dt","quarter"] if c in demo.columns]
            drug_keep = [c for c in ["primaryid","drugname"] if c in drug.columns]
            reac_keep = [c for c in ["primaryid","pt"] if c in reac.columns]

            demo_parts.append(demo[demo_keep])
            drug_parts.append(drug[drug_keep])
            reac_parts.append(reac[reac_keep])

            total_cases     += len(demo)
            total_drug_rows += len(drug)
            total_reac_rows += len(reac)

            pbar.write(f"  [{_ts()}] ✓ {quarter}: "
                       f"{len(demo):,} cases | {len(drug):,} drug rows | {len(reac):,} reactions")

        except zipfile.BadZipFile:
            pbar.write(f"  [{_ts()}] ✗ {quarter}: corrupt zip — skipping")
            parse_errors.append(quarter)
        except Exception as e:
            pbar.write(f"  [{_ts()}] ✗ {quarter}: {type(e).__name__}: {e} — skipping")
            parse_errors.append(quarter)

    if not demo_parts:
        log("No data parsed successfully.", "ERR")
        return None, None, None

    log("Concatenating and deduplicating...", "STEP")
    demo_df = pd.concat(demo_parts, ignore_index=True)
    drug_df = pd.concat(drug_parts, ignore_index=True)
    reac_df = pd.concat(reac_parts, ignore_index=True)

    n_before = len(demo_df)
    demo_df  = demo_df.drop_duplicates("primaryid")
    n_dupes  = n_before - len(demo_df)
    log(f"Removed {n_dupes:,} duplicate case records", "OK")

    log("Matching drug names to canonical list...", "STEP")
    drug_df["drugname_lower"]  = drug_df["drugname"].str.lower().fillna("")
    drug_df["drug_canonical"]  = None

    for canonical, aliases in tqdm(DRUG_ALIASES.items(), desc="  Drug matching",
                                   leave=False, file=sys.stdout):
        pattern = "|".join(re.escape(a) for a in aliases)
        mask    = drug_df["drugname_lower"].str.contains(pattern, na=False)
        drug_df.loc[mask, "drug_canonical"] = canonical

    n_matched = drug_df["drug_canonical"].notna().sum()
    drug_df   = drug_df[drug_df["drug_canonical"].notna()].copy()
    log(f"Drug rows matched to 12-drug list: {n_matched:,}", "OK")

    log("Writing parquet cache...", "STEP")
    demo_df.to_parquet(_CACHE_DIR / "demo.parquet")
    drug_df.to_parquet(_CACHE_DIR / "drug.parquet")
    reac_df.to_parquet(_CACHE_DIR / "reac.parquet")

    elapsed = time.time() - t_start
    m, s    = divmod(int(elapsed), 60)

    log_section("PARSE SUMMARY")
    log(f"Quarters parsed:   {len(zip_files) - len(parse_errors)} / {len(zip_files)}", "OK")
    log(f"Total cases:       {len(demo_df):,}  (after dedup)")
    log(f"Drug rows (12):    {len(drug_df):,}")
    log(f"Reaction rows:     {len(reac_df):,}")
    log(f"Parse errors:      {len(parse_errors)}  {parse_errors if parse_errors else ''}")
    log(f"Runtime:           {m}m {s:02d}s", "OK")
    log(f"Cache written to:  {_CACHE_DIR}", "OK")

    return demo_df, drug_df, reac_df


def load_parsed():
    paths = [_CACHE_DIR / p for p in ["demo.parquet","drug.parquet","reac.parquet"]]
    if not all(p.exists() for p in paths):
        return None, None, None
    log("Loading parquet cache...", "STEP")
    demo = pd.read_parquet(paths[0])
    drug = pd.read_parquet(paths[1])
    reac = pd.read_parquet(paths[2])
    log(f"Loaded: {len(demo):,} cases | {len(drug):,} drug rows | {len(reac):,} reaction rows", "OK")
    return demo, drug, reac


# ─────────────────────────────────────────────────────────────────────────────
# 3. ANALYZE
# ─────────────────────────────────────────────────────────────────────────────

def compute_ror(a, b, c, d):
    a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    if b == 0 or c == 0:
        return np.nan, np.nan, np.nan
    ror     = (a / b) / (c / d)
    log_ror = np.log(ror)
    se      = np.sqrt(1/a + 1/b + 1/c + 1/d)
    ci_lo   = np.exp(log_ror - 1.96 * se)
    ci_hi   = np.exp(log_ror + 1.96 * se)
    return round(ror, 3), round(ci_lo, 3), round(ci_hi, 3)


def analyze(demo_df=None, drug_df=None, reac_df=None):
    log_section("ANALYZE")
    t_start = time.time()

    if demo_df is None:
        demo_df, drug_df, reac_df = load_parsed()
    if demo_df is None:
        log("No parsed data found. Run --parse first.", "ERR")
        return

    # ── Pediatric filter ──────────────────────────────────────────────────────
    log("Applying pediatric filter (age < 18)...", "STEP")
    ped_ids   = set(demo_df[demo_df["age_years"] < 18]["primaryid"])
    n_age_unk = demo_df["age_years"].isna().sum()
    log(f"Pediatric cases:       {len(ped_ids):,} / {len(demo_df):,} total")
    log(f"Age unknown/missing:   {n_age_unk:,}  (excluded conservatively)")

    drug_ped = drug_df[drug_df["primaryid"].isin(ped_ids)].copy()
    reac_ped = reac_df[reac_df["primaryid"].isin(ped_ids)].copy()

    # ── Cardiac events ────────────────────────────────────────────────────────
    log("Identifying cardiac adverse events...", "STEP")
    primary_cardiac_ids = set(reac_ped[reac_ped["pt"].isin(PRIMARY_CARDIAC_PTS)]["primaryid"])
    broad_cardiac_ids   = set(reac_ped[reac_ped["pt"].isin(CARDIAC_PTS)]["primaryid"])
    log(f"Ped cases with primary cardiac PT (QTc/TdP/arrest): {len(primary_cardiac_ids):,}")
    log(f"Ped cases with any cardiac PT (broad set):          {len(broad_cardiac_ids):,}")

    total_ped   = len(ped_ids)
    drug_case_map = (drug_ped.groupby("drug_canonical")["primaryid"].apply(set).to_dict())
    log(f"Drugs with ≥1 pediatric report: {len(drug_case_map)}/12")

    # ── Per-drug ROR ──────────────────────────────────────────────────────────
    log_section("PER-DRUG ROR")
    log("ROR > 1 with 95% CI lower bound > 1 = pharmacovigilance signal")
    log(f"{'Drug':22s}  {'n':>6s}  {'cardiac':>7s}  {'ROR':>6s}  {'95% CI':14s}  signal")
    log("─" * 70)

    drug_rows = []
    for drug in tqdm(sorted(drug_case_map), desc="Drug ROR", unit="drug",
                     leave=False, file=sys.stdout):
        drug_cases = drug_case_map[drug]
        n_drug = len(drug_cases)
        a = len(drug_cases & primary_cardiac_ids)
        b = n_drug - a
        c = len(primary_cardiac_ids - drug_cases)
        d = total_ped - n_drug - c

        ror, ci_lo, ci_hi = compute_ror(a, b, c, d)
        sig = (not np.isnan(ci_lo)) and (ci_lo > 1.0)

        a2 = len(drug_cases & broad_cardiac_ids)
        ror2, ci_lo2, ci_hi2 = compute_ror(a2, n_drug-a2,
                                            len(broad_cardiac_ids-drug_cases),
                                            total_ped-n_drug-len(broad_cardiac_ids-drug_cases))

        drug_rows.append({
            "drug": drug, "n_cases_with_drug": n_drug, "n_cardiac_events": a,
            "ROR_primary": ror, "CI_lo_primary": ci_lo, "CI_hi_primary": ci_hi,
            "signal_primary": sig, "ROR_broad": ror2,
            "CI_lo_broad": ci_lo2, "CI_hi_broad": ci_hi2,
        })

        flag   = "⚠ SIGNAL" if sig else "  ─"
        ci_str = f"{ci_lo:.2f}–{ci_hi:.2f}" if not np.isnan(ci_lo) else "n/a"
        log(f"  {drug:22s}  {n_drug:>6,}  {a:>7,}  "
            f"{ror if not np.isnan(ror) else 'n/a':>6}  {ci_str:14s}  {flag}")

    drug_ror_df = pd.DataFrame(drug_rows).sort_values("ROR_primary", ascending=False)
    drug_ror_df.to_csv(_RESULTS_DIR / "faers_drug_ror.csv", index=False)
    log(f"Saved → faers_drug_ror.csv", "OK")

    # ── Pairwise combo ROR ────────────────────────────────────────────────────
    log_section("PAIRWISE COMBO ROR")
    drugs      = list(drug_case_map.keys())
    all_pairs  = list(itertools.combinations(drugs, 2))
    log(f"Computing ROR for {len(all_pairs)} drug pairs (min n=3 for inclusion)")

    combo_rows = []
    skipped    = 0

    pbar = tqdm(all_pairs, desc="Combo ROR", unit="pair",
                bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]",
                file=sys.stdout)

    for drug_a, drug_b in pbar:
        combo_cases = drug_case_map[drug_a] & drug_case_map[drug_b]
        n_combo     = len(combo_cases)
        pbar.set_postfix_str(f"{drug_a[:6]}+{drug_b[:6]} n={n_combo}")

        if n_combo < 3:
            skipped += 1
            continue

        a = len(combo_cases & primary_cardiac_ids)
        b = n_combo - a
        c = len(primary_cardiac_ids - combo_cases)
        d = total_ped - n_combo - c

        ror, ci_lo, ci_hi = compute_ror(a, b, c, d)
        sig = (not np.isnan(ci_lo)) and (ci_lo > 1.0)

        combo_rows.append({
            "combo": f"{drug_a}+{drug_b}", "drug_A": drug_a, "drug_B": drug_b,
            "n_co_reported": n_combo, "n_cardiac_events": a,
            "ROR": ror, "CI_lo": ci_lo, "CI_hi": ci_hi, "signal": sig,
        })

    combo_ror_df = pd.DataFrame(combo_rows).sort_values("ROR", ascending=False)
    combo_ror_df.to_csv(_RESULTS_DIR / "faers_combo_ror.csv", index=False)

    n_sig = combo_ror_df["signal"].sum() if not combo_ror_df.empty else 0
    log(f"Pairs evaluated: {len(combo_rows)}  |  skipped (n<3): {skipped}  |  signals: {n_sig}", "OK")

    if not combo_ror_df.empty:
        log("Top 10 combos by ROR:")
        for _, row in combo_ror_df.head(10).iterrows():
            flag   = "⚠" if row["signal"] else " "
            ci_str = f"{row['CI_lo']:.2f}–{row['CI_hi']:.2f}" if not np.isnan(row['CI_lo']) else "n/a"
            log(f"  {flag} {row['combo']:42s}  n={row['n_co_reported']:4,}  "
                f"ROR={row['ROR']:.2f} [{ci_str}]")

    log(f"Saved → faers_combo_ror.csv", "OK")

    # ── Model alignment ───────────────────────────────────────────────────────
    log_section("MODEL vs FAERS ALIGNMENT")

    # risk_grid_results.csv uses abbreviations (MPH+ARI); combo_ror_df uses full names
    ABBREV_TO_FULL = {
        "MPH": "Methylphenidate", "AMP": "Amphetamine",
        "RIS": "Risperidone",     "QUE": "Quetiapine",   "ARI": "Aripiprazole",
        "SER": "Sertraline",      "FLU": "Fluoxetine",   "ESC": "Escitalopram",
        "CLO": "Clonidine",       "GUA": "Guanfacine",
        "IMI": "Imipramine",      "NOR": "Nortriptyline",
    }

    results_path = _ROOT_DIR / "results" / "risk_grid_results.csv"
    if results_path.exists():
        model_df  = pd.read_csv(results_path)
        model_map = dict(zip(model_df["combination"], model_df["\u0394QTc_ms"]))
        log(f"Loaded live model results from risk_grid_results.csv ({len(model_map)} combos)", "OK")
    else:
        model_map = MODEL_DQTC
        log("risk_grid_results.csv not found — using hardcoded predictions", "WARN")

    align_rows = []
    for combo_key, dqtc in model_map.items():
        parts = combo_key.split("+")
        if len(parts) != 2:
            continue
        abbr_a, abbr_b = parts
        da = ABBREV_TO_FULL.get(abbr_a, abbr_a)
        db = ABBREV_TO_FULL.get(abbr_b, abbr_b)

        match = combo_ror_df[
            ((combo_ror_df["drug_A"]==da) & (combo_ror_df["drug_B"]==db)) |
            ((combo_ror_df["drug_A"]==db) & (combo_ror_df["drug_B"]==da))
        ] if not combo_ror_df.empty else pd.DataFrame()

        if match.empty:
            ror = ci_lo = ci_hi = n_co = np.nan
            sig = False
        else:
            r   = match.iloc[0]
            ror, ci_lo, ci_hi, n_co, sig = r["ROR"], r["CI_lo"], r["CI_hi"], r["n_co_reported"], r["signal"]

        tier = ("HIGH" if dqtc>=20 else "MODERATE" if dqtc>=10
                else "LOW-MOD" if dqtc>=5 else "LOW")
        concordant = (tier in ("HIGH","MODERATE") and sig) or (tier in ("LOW","LOW-MOD") and not sig)

        align_rows.append({
            "combination": combo_key, "model_dQTc_ms": dqtc, "model_risk_tier": tier,
            "faers_n_coreported": n_co, "faers_ROR": ror,
            "faers_CI_lo": ci_lo, "faers_CI_hi": ci_hi,
            "faers_signal": sig, "concordant": concordant,
        })

    align_df = pd.DataFrame(align_rows).sort_values("model_dQTc_ms", ascending=False)
    align_df.to_csv(_RESULTS_DIR / "faers_model_alignment.csv", index=False)

    # Print alignment table
    header = f"{'Combination':42s} {'dQTc':>6s}  {'Tier':8s}  {'ROR':>6s}  {'95% CI':14s}  {'Sig':3s}  {'Match':5s}"
    log(header)
    log("─" * 80)
    for _, row in align_df.iterrows():
        ci_str = (f"{row['faers_CI_lo']:.2f}–{row['faers_CI_hi']:.2f}"
                  if not np.isnan(row['faers_CI_lo']) else "n/a         ")
        ror_str  = f"{row['faers_ROR']:.2f}" if not np.isnan(row['faers_ROR']) else "n/a "
        sig_str  = "YES" if row["faers_signal"] else "no "
        match_ch = "✓" if row["concordant"] else "✗"
        log(f"  {row['combination']:40s} {row['model_dQTc_ms']:+6.1f}  "
            f"{row['model_risk_tier']:8s}  {ror_str:>6s}  {ci_str:14s}  "
            f"{sig_str:3s}  {match_ch}")

    n_with_faers = len(align_df.dropna(subset=["faers_ROR"]))
    n_concordant = int(align_df["concordant"].sum())
    log("─" * 80)
    log(f"Concordance: {n_concordant}/{n_with_faers} pairs with FAERS data", "OK")

    # PT breakdown
    log_section("CARDIAC PT BREAKDOWN (high-risk drugs)")
    high_drugs = {"Methylphenidate","Amphetamine","Aripiprazole",
                  "Risperidone","Imipramine","Nortriptyline"}
    high_ids   = set().union(*[drug_case_map[d] for d in high_drugs if d in drug_case_map])
    pt_counts  = (
        reac_ped[reac_ped["primaryid"].isin(high_ids) & reac_ped["pt"].isin(CARDIAC_PTS)]
        ["pt"].value_counts()
    )
    for pt, cnt in pt_counts.items():
        log(f"  {pt:40s}  {cnt:,}")

    elapsed = time.time() - t_start
    m, s    = divmod(int(elapsed), 60)
    log_section("DONE")
    log(f"Analyze runtime: {m}m {s:02d}s", "OK")
    log(f"Outputs written to: {_RESULTS_DIR}", "OK")
    log("Files: faers_drug_ror.csv | faers_combo_ror.csv | faers_model_alignment.csv", "OK")

    return drug_ror_df, combo_ror_df, align_df


# ─────────────────────────────────────────────────────────────────────────────
# 4. CLI
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="FAERS pediatric cardiac safety analysis — CardioSafe Pediatric",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 faers_pipeline.py --all
  python3 faers_pipeline.py --download --quarters 2022q1 2022q2 2022q3 2022q4
  python3 faers_pipeline.py --parse
  python3 faers_pipeline.py --analyze
  python3 faers_pipeline.py --download --force    # re-download even if cached
        """
    )
    parser.add_argument("--download", action="store_true", help="Download FAERS quarterly zips")
    parser.add_argument("--parse",    action="store_true", help="Parse zips into parquet cache")
    parser.add_argument("--analyze",  action="store_true", help="Run ROR analysis + model alignment")
    parser.add_argument("--all",      action="store_true", help="Full pipeline: download + parse + analyze")
    parser.add_argument("--quarters", nargs="+", default=None,
                        help="Override quarter list, e.g. --quarters 2022q1 2022q2")
    parser.add_argument("--force",    action="store_true", help="Re-download even if cached")
    args = parser.parse_args()

    quarters = args.quarters if args.quarters else QUARTERS

    if not any([args.download, args.parse, args.analyze, args.all]):
        parser.print_help()
        return

    t0 = time.time()

    if args.all or args.download:
        download_faers(quarters, force=args.force)
    if args.all or args.parse:
        parse_faers(quarters)
    if args.all or args.analyze:
        analyze()

    total = time.time() - t0
    m, s  = divmod(int(total), 60)
    print(f"\n[{_ts()}] ══ PIPELINE COMPLETE — total runtime {m}m {s:02d}s")


if __name__ == "__main__":
    main()