"""
FAERS Temporal Trend + Age Stratification Analysis — CardioSafe Pediatric
=========================================================================
Two secondary analyses on the existing parquet cache:

1. TEMPORAL TREND — Year-by-year ROR for key combinations (MPH+SER, MPH+ARI,
   MPH+QUE, ARI+GUA) to test whether cardiac AE signal strength tracks the
   growth of pediatric polypharmacy over 2015–2024.

2. AGE STRATIFICATION — Split pediatric cases into children (6–12y) vs.
   adolescents (13–17y) and compute ROR separately. Tests whether stimulant
   cardiac signal is age-specific.

Input:  data/faers_cache/demo.parquet, drug.parquet, reac.parquet
Output: results/faers/temporal_trend.csv
        results/faers/age_stratification.csv
        results/faers/temporal_age_memo.md

Usage:
    python3 src/faers_secondary.py
"""

import sys
import time
import warnings
import numpy as np
import pandas as pd
from pathlib import Path

warnings.filterwarnings("ignore")

_SRC_DIR   = Path(__file__).resolve().parent
_ROOT_DIR  = _SRC_DIR.parent
_CACHE_DIR = _ROOT_DIR / "data" / "faers_cache"
_FAERS_DIR = _ROOT_DIR / "results" / "faers"
_FAERS_DIR.mkdir(parents=True, exist_ok=True)

def ts(): return time.strftime("%H:%M:%S")
def log(msg): print(f"[{ts()}] {msg}", flush=True)

# ── COMBOS OF INTEREST ────────────────────────────────────────────────────────
COMBOS = {
    "MPH+SER":  ("Methylphenidate", "Sertraline"),
    "MPH+ARI":  ("Methylphenidate", "Aripiprazole"),
    "MPH+QUE":  ("Methylphenidate", "Quetiapine"),
    "ARI+GUA":  ("Aripiprazole",    "Guanfacine"),
    "QUE+GUA":  ("Quetiapine",      "Guanfacine"),
    "QUE+FLU":  ("Quetiapine",      "Fluoxetine"),
    "MPH+RIS":  ("Methylphenidate", "Risperidone"),
    "AMP+SER":  ("Amphetamine",     "Sertraline"),
}

PRIMARY_CARDIAC_PTS = {
    "QT prolonged", "Electrocardiogram QT prolonged", "Torsade de pointes",
    "Ventricular tachycardia", "Ventricular fibrillation", "Cardiac arrest",
    "Sudden cardiac death", "Long QT syndrome", "Ventricular arrhythmia",
}

def compute_ror(a, b, c, d):
    a, b, c, d = a+0.5, b+0.5, c+0.5, d+0.5
    if b == 0 or c == 0: return np.nan, np.nan, np.nan
    ror    = (a/b)/(c/d)
    se     = np.sqrt(1/a + 1/b + 1/c + 1/d)
    log_r  = np.log(ror)
    return round(ror,3), round(np.exp(log_r-1.96*se),3), round(np.exp(log_r+1.96*se),3)

def quarter_to_year(q_str):
    """'2019q3' -> 2019"""
    return int(str(q_str)[:4])

def load_cache():
    log("Loading parquet cache...")
    demo = pd.read_parquet(_CACHE_DIR / "demo.parquet")
    drug = pd.read_parquet(_CACHE_DIR / "drug.parquet")
    reac = pd.read_parquet(_CACHE_DIR / "reac.parquet")
    log(f"Loaded: {len(demo):,} cases | {len(drug):,} drug rows | {len(reac):,} reactions")
    return demo, drug, reac

# ─────────────────────────────────────────────────────────────────────────────
# 1. TEMPORAL TREND
# ─────────────────────────────────────────────────────────────────────────────

def temporal_trend(demo, drug, reac):
    log("=" * 65)
    log("TEMPORAL TREND ANALYSIS")
    log("=" * 65)

    # Pediatric filter
    ped_ids = set(demo[demo["age_years"] < 18]["primaryid"])
    log(f"Pediatric cases: {len(ped_ids):,}")

    # Add year to demo
    demo_ped = demo[demo["primaryid"].isin(ped_ids)].copy()
    if "quarter" in demo_ped.columns:
        demo_ped["year"] = demo_ped["quarter"].apply(quarter_to_year)
    else:
        log("WARNING: no quarter column — cannot do temporal analysis")
        return pd.DataFrame()

    # Cardiac event flag
    cardiac_ids = set(reac[reac["pt"].isin(PRIMARY_CARDIAC_PTS)]["primaryid"])

    # Drug case maps per year
    drug_ped = drug[drug["primaryid"].isin(ped_ids)].copy()

    rows = []
    years = sorted(demo_ped["year"].dropna().unique().astype(int))

    log(f"\nYears available: {years[0]}–{years[-1]}")
    log(f"\n{'Combo':12s} {'Year':>5s}  {'n_combo':>8s}  {'cardiac':>7s}  {'ROR':>6s}  {'95% CI':14s}  signal")
    log("─" * 65)

    for combo_key, (drug_a, drug_b) in COMBOS.items():
        for year in years:
            year_ids = set(demo_ped[demo_ped["year"] == year]["primaryid"])
            drug_year = drug_ped[drug_ped["primaryid"].isin(year_ids)]

            cases_a = set(drug_year[drug_year["drug_canonical"]==drug_a]["primaryid"])
            cases_b = set(drug_year[drug_year["drug_canonical"]==drug_b]["primaryid"])
            combo_cases = cases_a & cases_b
            n_combo = len(combo_cases)

            if n_combo < 3:
                rows.append({"combo":combo_key,"year":year,"n_combo":n_combo,
                             "n_cardiac":0,"ROR":np.nan,"CI_lo":np.nan,
                             "CI_hi":np.nan,"signal":False})
                continue

            year_cardiac = cardiac_ids & year_ids
            a = len(combo_cases & year_cardiac)
            b = n_combo - a
            c = len(year_cardiac - combo_cases)
            d = len(year_ids) - n_combo - c

            ror, ci_lo, ci_hi = compute_ror(a, b, c, d)
            sig = (not np.isnan(ci_lo)) and (ci_lo > 1.0)

            rows.append({"combo":combo_key,"year":year,"n_combo":n_combo,
                         "n_cardiac":a,"ROR":ror,"CI_lo":ci_lo,
                         "CI_hi":ci_hi,"signal":sig})

            sig_str = "⚠ YES" if sig else "  no"
            ci_str  = f"{ci_lo:.2f}–{ci_hi:.2f}" if not np.isnan(ci_lo) else "n/a"
            ror_str = f"{ror:.2f}" if not np.isnan(ror) else "n/a"
            log(f"{combo_key:12s} {year:>5d}  {n_combo:>8,}  {a:>7,}  "
                f"{ror_str:>6s}  {ci_str:14s}  {sig_str}")

    df = pd.DataFrame(rows)
    df.to_csv(_FAERS_DIR / "temporal_trend.csv", index=False)
    log(f"\nSaved → results/faers/temporal_trend.csv")

    # Trend summary
    log("\nTREND SUMMARY — signal years per combo:")
    for combo_key in COMBOS:
        sub = df[df["combo"]==combo_key]
        sig_years = sub[sub["signal"]==True]["year"].tolist()
        n_max = sub["n_combo"].max()
        log(f"  {combo_key:12s}  signal years: {sig_years}  |  peak n={n_max:,}")

    return df


# ─────────────────────────────────────────────────────────────────────────────
# 2. AGE STRATIFICATION
# ─────────────────────────────────────────────────────────────────────────────

def age_stratification(demo, drug, reac):
    log("\n" + "=" * 65)
    log("AGE STRATIFICATION ANALYSIS")
    log("Children (6–12y) vs. Adolescents (13–17y)")
    log("=" * 65)

    cardiac_ids = set(reac[reac["pt"].isin(PRIMARY_CARDIAC_PTS)]["primaryid"])

    age_groups = {
        "children (6-12y)":    (6,  12),
        "adolescents (13-17y)":(13, 17),
        "all pediatric (<18y)":(0,  17),
    }

    rows = []

    for group_label, (age_lo, age_hi) in age_groups.items():
        group_ids = set(demo[
            (demo["age_years"] >= age_lo) &
            (demo["age_years"] <= age_hi)
        ]["primaryid"])
        n_group = len(group_ids)

        drug_grp = drug[drug["primaryid"].isin(group_ids)]
        drug_case_map = drug_grp.groupby("drug_canonical")["primaryid"].apply(set).to_dict()

        log(f"\n── {group_label} (n={n_group:,}) ──")
        log(f"{'Combo':12s}  {'n_combo':>8s}  {'cardiac':>7s}  {'ROR':>6s}  {'95% CI':14s}  signal")
        log("─" * 60)

        for combo_key, (drug_a, drug_b) in COMBOS.items():
            cases_a = drug_case_map.get(drug_a, set())
            cases_b = drug_case_map.get(drug_b, set())
            combo_cases = cases_a & cases_b
            n_combo = len(combo_cases)

            if n_combo < 3:
                rows.append({"age_group":group_label,"combo":combo_key,
                             "n_group":n_group,"n_combo":n_combo,
                             "n_cardiac":0,"ROR":np.nan,"CI_lo":np.nan,
                             "CI_hi":np.nan,"signal":False})
                continue

            group_cardiac = cardiac_ids & group_ids
            a = len(combo_cases & group_cardiac)
            b = n_combo - a
            c = len(group_cardiac - combo_cases)
            d = n_group - n_combo - c

            ror, ci_lo, ci_hi = compute_ror(a, b, c, d)
            sig = (not np.isnan(ci_lo)) and (ci_lo > 1.0)

            rows.append({"age_group":group_label,"combo":combo_key,
                         "n_group":n_group,"n_combo":n_combo,"n_cardiac":a,
                         "ROR":ror,"CI_lo":ci_lo,"CI_hi":ci_hi,"signal":sig})

            sig_str = "⚠ YES" if sig else "  no"
            ci_str  = f"{ci_lo:.2f}–{ci_hi:.2f}" if not np.isnan(ci_lo) else "n/a"
            ror_str = f"{ror:.2f}" if not np.isnan(ror) else "n/a"
            log(f"{combo_key:12s}  {n_combo:>8,}  {a:>7,}  "
                f"{ror_str:>6s}  {ci_str:14s}  {sig_str}")

    df = pd.DataFrame(rows)
    df.to_csv(_FAERS_DIR / "age_stratification.csv", index=False)
    log(f"\nSaved → results/faers/age_stratification.csv")

    # Cross-group comparison
    log("\nCROSS-GROUP ROR COMPARISON:")
    log(f"{'Combo':12s}  {'Children ROR':>13s}  {'Adolescents ROR':>16s}  {'All Ped ROR':>12s}  Direction")
    log("─" * 70)
    for combo_key in COMBOS:
        sub = df[df["combo"]==combo_key]
        def get_ror(label):
            r = sub[sub["age_group"]==label]
            if r.empty or np.isnan(r.iloc[0]["ROR"]): return "n/a"
            return f"{r.iloc[0]['ROR']:.2f}"
        def get_ror_val(label):
            r = sub[sub["age_group"]==label]
            if r.empty: return np.nan
            return r.iloc[0]["ROR"]

        ch = get_ror("children (6-12y)")
        ad = get_ror("adolescents (13-17y)")
        al = get_ror("all pediatric (<18y)")
        ch_v = get_ror_val("children (6-12y)")
        ad_v = get_ror_val("adolescents (13-17y)")

        if np.isnan(ch_v) or np.isnan(ad_v):
            direction = "insufficient data"
        elif ad_v > ch_v * 1.3:
            direction = "↑ HIGHER in adolescents"
        elif ch_v > ad_v * 1.3:
            direction = "↑ HIGHER in children"
        else:
            direction = "similar across groups"

        log(f"{combo_key:12s}  {ch:>13s}  {ad:>16s}  {al:>12s}  {direction}")

    return df


# ─────────────────────────────────────────────────────────────────────────────
# 3. WRITE MEMO
# ─────────────────────────────────────────────────────────────────────────────

def write_memo(trend_df, age_df):
    if trend_df.empty and age_df.empty:
        return

    memo = """# FAERS Secondary Analyses — Temporal Trend & Age Stratification
**CardioSafe Pediatric | June 2026**

---

## 1. Temporal Trend Analysis

Year-by-year ROR for key drug combinations in pediatric FAERS cases (2015–2024).
Tests whether pharmacovigilance signal strength has changed over time as
pediatric polypharmacy has grown.

### Methods
Per-year ROR computed for 8 key combinations using the same pediatric filter
(age < 18) and primary cardiac PT definition as the main FAERS analysis.
Minimum n=3 co-reported cases required for ROR computation.

### Key Findings
"""

    if not trend_df.empty:
        for combo in COMBOS:
            sub = trend_df[trend_df["combo"]==combo].dropna(subset=["ROR"])
            if sub.empty: continue
            sig_years = sub[sub["signal"]==True]["year"].tolist()
            n_trend = sub.set_index("year")["n_combo"]
            first_yr = n_trend.index.min()
            last_yr  = n_trend.index.max()
            n_first  = n_trend.get(first_yr, 0)
            n_last   = n_trend.get(last_yr, 0)
            if n_first > 0:
                pct = (n_last - n_first) / n_first * 100
                trend_str = f"{pct:+.0f}% change in co-reports {first_yr}→{last_yr}"
            else:
                trend_str = "insufficient early data"
            memo += f"\n**{combo}**: Signal years: {sig_years or 'none'}. {trend_str}.\n"

    memo += """
---

## 2. Age Stratification Analysis

Comparison of ROR between children (6–12y) and adolescents (13–17y) for key
combinations to test whether cardiac AE signal is age-specific.

### Methods
Same pediatric cases split into two age bands. Cases with unknown age excluded
(as in main analysis). Same primary cardiac PT definition.

### Key Findings
"""

    if not age_df.empty:
        for combo in COMBOS:
            sub = age_df[age_df["combo"]==combo]
            def get(label):
                r = sub[sub["age_group"]==label]
                if r.empty or np.isnan(r.iloc[0]["ROR"]): return None
                return r.iloc[0]
            ch = get("children (6-12y)")
            ad = get("adolescents (13-17y)")
            if ch is None and ad is None: continue
            ch_str = f"ROR {ch['ROR']:.2f} [{ch['CI_lo']:.2f}–{ch['CI_hi']:.2f}] n={int(ch['n_combo'])}" if ch is not None else "n/a"
            ad_str = f"ROR {ad['ROR']:.2f} [{ad['CI_lo']:.2f}–{ad['CI_hi']:.2f}] n={int(ad['n_combo'])}" if ad is not None else "n/a"
            memo += f"\n**{combo}**: Children: {ch_str}. Adolescents: {ad_str}.\n"

    memo += """
---

## Limitations

- Age is unknown for ~43% of FAERS cases; unknown-age cases excluded
  conservatively and may be enriched for pediatric patients
- Small n per year-group for rare combinations limits temporal analysis
  statistical power
- Reporting patterns may have changed over time (FDA MedWatch campaigns,
  EHR integration) independent of true risk changes
- Age band definitions (6–12 vs 13–17) are arbitrary; puberty onset varies
"""

    (_FAERS_DIR / "temporal_age_memo.md").write_text(memo)
    log(f"Saved → results/faers/temporal_age_memo.md")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    t0 = time.time()
    demo, drug, reac = load_cache()

    trend_df = temporal_trend(demo, drug, reac)
    age_df   = age_stratification(demo, drug, reac)
    write_memo(trend_df, age_df)

    elapsed = time.time() - t0
    m, s = divmod(int(elapsed), 60)
    log(f"\nComplete in {m}m {s:02d}s")
    log("Outputs: temporal_trend.csv | age_stratification.csv | temporal_age_memo.md")
