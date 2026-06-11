#!/usr/bin/env python3
"""
CardioSafe Pediatric - FAERS concordance recomputation vs IKs risk grid.

Tests whether the model-FAERS validation survives the IKs sweep, with the
key breakdown: stimulant-containing vs non-stimulant combinations.

Usage:
    python3 concordance_iks.py \
        --grid results/risk_grid_results.csv \
        --faers results/faers/faers_combo_ror.csv \
        [--old-grid path/to/pre_iks_grid.csv]   # optional, for old-vs-new
"""
import argparse
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata

ABBREV = {
    "methylphenidate": "MPH", "amphetamine": "AMP", "risperidone": "RIS",
    "quetiapine": "QUE", "aripiprazole": "ARI", "sertraline": "SER",
    "fluoxetine": "FLU", "escitalopram": "ESC", "clonidine": "CLO",
    "guanfacine": "GUA", "imipramine": "IMI", "nortriptyline": "NOR",
}
KNOWN = set(ABBREV.values())
STIMULANTS = {"MPH", "AMP"}


def to_abbrev(token):
    t = token.strip().lower()
    if t.upper() in KNOWN:
        return t.upper()
    return ABBREV.get(t, token.strip().upper()[:3])


def norm_combo(s):
    """Order-independent canonical key, e.g. 'Methylphenidate + Aripiprazole' -> 'ARI+MPH'."""
    parts = [p for p in str(s).replace("+", " ").replace(",", " ").split() if p]
    codes = sorted({to_abbrev(p) for p in parts})
    return "+".join(codes)


def find_col(df, candidates, required=True, what=""):
    low = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in low:
            return low[cand.lower()]
    # fuzzy contains
    for c in df.columns:
        cl = c.lower()
        if any(cand.lower() in cl for cand in candidates):
            return c
    if required:
        sys.exit(f"ERROR: could not find {what} column. Looked for {candidates}. "
                 f"Columns present: {list(df.columns)}")
    return None


def auc_mw(scores, labels):
    """AUC via Mann-Whitney U (no sklearn dependency). labels in {0,1}."""
    scores = np.asarray(scores, float)
    labels = np.asarray(labels, int)
    pos = scores[labels == 1]
    neg = scores[labels == 0]
    n1, n0 = len(pos), len(neg)
    if n1 == 0 or n0 == 0:
        return np.nan, n1, n0
    r = rankdata(np.concatenate([pos, neg]))
    auc = (r[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
    return auc, n1, n0


def report(label, df, dqtc, ikr, signal, ror):
    print(f"\n{'='*60}\n{label}  (n={len(df)})\n{'='*60}")
    if len(df) < 3:
        print("  too few rows for stats")
        return
    n_sig = int(df[signal].sum())
    print(f"  signals: {n_sig} / {len(df)}")

    # AUC: predictor discriminates signal vs non-signal
    auc_d, n1, n0 = auc_mw(df[dqtc], df[signal])
    print(f"  AUC (delta-QTc -> signal):  {auc_d:.3f}   [{n1} sig / {n0} non]")
    if ikr is not None:
        auc_i, _, _ = auc_mw(df[ikr], df[signal])
        print(f"  AUC (IKr%   -> signal):     {auc_i:.3f}")

    # Spearman: continuous delta-QTc vs continuous ROR
    if ror is not None and df[ror].notna().sum() >= 3:
        rho, p = spearmanr(df[dqtc], df[ror], nan_policy="omit")
        print(f"  Spearman (delta-QTc vs ROR): rho={rho:+.3f}  p={p:.3g}")


def load_grid(path):
    g = pd.read_csv(path)
    combo_c = find_col(g, ["combo", "combination", "drugs", "pair", "name"],
                       required=False, what="grid combo")
    if combo_c is None:
        d1 = find_col(g, ["drug1", "drug_1", "a"], what="grid drug1")
        d2 = find_col(g, ["drug2", "drug_2", "b"], what="grid drug2")
        g["_combo_raw"] = g[d1].astype(str) + "+" + g[d2].astype(str)
        combo_c = "_combo_raw"
    dqtc_c = find_col(g, ["dqtc", "delta_qtc", "qtc_delta", "ΔQTc", "deltaqtc", "dQTc"],
                      what="grid delta-QTc")
    ikr_c = find_col(g, ["ikr", "ikr_pct", "ikr_block", "ikr_percent"],
                     required=False, what="grid IKr")
    g["key"] = g[combo_c].map(norm_combo)
    g = g.rename(columns={dqtc_c: "dqtc"})
    if ikr_c:
        g = g.rename(columns={ikr_c: "ikr"})
    keep = ["key", "dqtc"] + (["ikr"] if ikr_c else [])
    return g[keep].drop_duplicates("key"), bool(ikr_c)


def load_faers(path):
    f = pd.read_csv(path)
    combo_c = find_col(f, ["combo", "combination", "drugs", "pair", "name"],
                       what="faers combo")
    ror_c = find_col(f, ["ror", "reporting_odds_ratio", "ror_value"],
                     required=False, what="faers ROR")
    sig_c = find_col(f, ["signal", "is_signal", "ror_signal", "flag"],
                     required=False, what="faers signal")
    f["key"] = f[combo_c].map(norm_combo)
    ren = {}
    if ror_c:
        ren[ror_c] = "ror"
    f = f.rename(columns=ren)
    if sig_c:
        f = f.rename(columns={sig_c: "signal"})
        # coerce to 0/1
        f["signal"] = (f["signal"].astype(str).str.lower()
                       .isin(["1", "true", "yes", "signal", "y"])).astype(int)
    elif "ror" in f.columns:
        # define signal as ROR lower-CI > 1 proxy: ROR >= 2 if no CI available
        f["signal"] = (f["ror"] >= 2.0).astype(int)
        print("NOTE: no signal column found; using ROR>=2 as signal proxy.")
    keep = ["key"] + [c for c in ["ror", "signal"] if c in f.columns]
    return f[keep].drop_duplicates("key")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", required=True)
    ap.add_argument("--faers", required=True)
    ap.add_argument("--old-grid", default=None)
    args = ap.parse_args()

    grid, has_ikr = load_grid(args.grid)
    faers = load_faers(args.faers)
    m = grid.merge(faers, on="key", how="inner")
    if "signal" not in m.columns:
        sys.exit("No signal column derivable from FAERS file; cannot compute AUC.")

    m["is_stim"] = m["key"].apply(lambda k: bool(set(k.split("+")) & STIMULANTS))

    ikr_col = "ikr" if has_ikr else None
    ror_col = "ror" if "ror" in m.columns else None

    print(f"\nMatched {len(m)} combinations (grid {len(grid)} x FAERS {len(faers)}).")
    report("ALL COMBINATIONS", m, "dqtc", ikr_col, "signal", ror_col)
    report("STIMULANT-CONTAINING (MPH or AMP)", m[m.is_stim], "dqtc", ikr_col, "signal", ror_col)
    report("NON-STIMULANT", m[~m.is_stim], "dqtc", ikr_col, "signal", ror_col)

    # Old-vs-new delta-QTc on the FAERS-signal combos
    if args.old_grid:
        old, _ = load_grid(args.old_grid)
        cmp = (m.merge(old.rename(columns={"dqtc": "dqtc_old"}), on="key", how="left"))
        sig = cmp[cmp.signal == 1].sort_values("dqtc_old", ascending=False)
        print(f"\n{'='*60}\nOLD vs NEW delta-QTc on FAERS-SIGNAL combos\n{'='*60}")
        print(f"  {'combo':<14}{'old':>8}{'new':>8}{'  stim':>7}")
        for _, r in sig.iterrows():
            print(f"  {r.key:<14}{r.dqtc_old:>8.1f}{r.dqtc:>8.1f}{('  yes' if r.is_stim else '  no'):>7}")

    print("\nInterpretation guide:")
    print("  - If STIM AUC ~0.5 and STIM Spearman flat/negative -> stimulant validation broke.")
    print("  - If NON-STIM AUC stays high -> the hERG/ARI story now carries the concordance.")


if __name__ == "__main__":
    main()
