#!/usr/bin/env python3
"""
CardioSafe Pediatric - IKs upregulation sensitivity sweep.

Sweeps the IKs upregulation factor in the sympathomimetic pathway and, at each
value, recomputes (a) the risk structure for tracked combos and (b) the full
FAERS concordance (delta-QTc AUC, IKr AUC, Spearman vs ROR) broken out by
stimulant / non-stimulant. Tells you whether the 0.632 AUC / +0.503 stimulant
Spearman are stable across plausible IKs values or a knife-edge artifact of one
setting.

REQUIRES a 2-line hook in ord_model.py (see apply_iks_factor docstring).

Usage (run in tmux, this is many sims):
    python3 src/iks_sensitivity.py \
        --faers results/faers/faers_combo_ror.csv \
        --params data/herg_master_params.csv \
        --factors 1.0,1.25,1.5,1.75,2.0,2.5 \
        --beats 200 \
        --out results/iks_sensitivity.csv
"""
import argparse
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata

sys.path.insert(0, "src")
import ord_model  # noqa: E402

ABBREV = {
    "methylphenidate": "MPH", "amphetamine": "AMP", "risperidone": "RIS",
    "quetiapine": "QUE", "aripiprazole": "ARI", "sertraline": "SER",
    "fluoxetine": "FLU", "escitalopram": "ESC", "clonidine": "CLO",
    "guanfacine": "GUA", "imipramine": "IMI", "nortriptyline": "NOR",
}
FULLNAME = {v: k.capitalize() for k, v in ABBREV.items()}
FULLNAME.update({"MPH": "Methylphenidate"})  # capitalize() handles the rest
KNOWN = set(ABBREV.values())
STIMULANTS = {"MPH", "AMP"}

# Combos to plot as markers across the sweep (key form, order-independent)
TRACKED = ["AMP+MPH", "ARI+MPH", "ARI+NOR", "ARI+QUE", "ARI+RIS", "ARI+IMI"]


def to_abbrev(tok):
    t = tok.strip().lower()
    if t.upper() in KNOWN:
        return t.upper()
    return ABBREV.get(t, tok.strip().upper()[:3])


def norm_combo(s):
    parts = [p for p in str(s).replace("+", " ").replace(",", " ").split() if p]
    return "+".join(sorted({to_abbrev(p) for p in parts}))


def key_to_drugdict(key, dose="therapeutic"):
    out = {}
    for code in key.split("+"):
        full = FULLNAME.get(code)
        if full is None:
            return None
        out[full] = dose
    return out


def auc_mw(scores, labels):
    scores, labels = np.asarray(scores, float), np.asarray(labels, int)
    pos, neg = scores[labels == 1], scores[labels == 0]
    n1, n0 = len(pos), len(neg)
    if n1 == 0 or n0 == 0:
        return np.nan
    r = rankdata(np.concatenate([pos, neg]))
    return (r[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def find_col(df, cands, required=True, what=""):
    low = {c.lower(): c for c in df.columns}
    for c in cands:
        if c.lower() in low:
            return low[c.lower()]
    for c in df.columns:
        if any(x.lower() in c.lower() for x in cands):
            return c
    if required:
        sys.exit(f"ERROR: no {what} column. Tried {cands}. Have {list(df.columns)}")
    return None


def apply_iks_factor(factor):
    """
    Set the IKs upregulation factor for the next run_simulation calls.

    REQUIRED ord_model.py change (2 lines). At module level add:
        IKS_UPREG_OVERRIDE = None
    Then wherever the sympathomimetic pathway applies IKs upregulation, replace
    your hardcoded factor, e.g.:
        gks_scale = 1.5
    with:
        gks_scale = IKS_UPREG_OVERRIDE if IKS_UPREG_OVERRIDE is not None else 1.5

    If your IKs term is not a simple GKs multiplier (e.g. you scale gating or add
    a fixed conductance bump), send me those lines and I will match the injection.
    """
    ord_model.IKS_UPREG = factor


def run_qtc_ikr(drugs, params, beats):
    """Return (QTc, IKr_block_pct_or_nan). drugs=None for baseline."""
    res = ord_model.run_simulation(drugs, params, n_beats=beats, verbose=False)
    qtc = res.get("QTc", np.nan)
    ikr = res.get("IKr_block_pct", np.nan)
    return float(qtc), float(ikr) if ikr is not None else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--faers", required=True)
    ap.add_argument("--params", default="data/herg_master_params.csv")
    ap.add_argument("--factors", default="1.0,1.25,1.5,1.75,2.0,2.5")
    ap.add_argument("--beats", type=int, default=200)
    ap.add_argument("--out", default="results/iks_sensitivity.csv")
    ap.add_argument("--fig", default="docs/figures/figureS_iks_sensitivity.png")
    ap.add_argument("--dose", default="therapeutic")
    args = ap.parse_args()

    factors = [float(x) for x in args.factors.split(",")]

    f = pd.read_csv(args.faers)
    combo_c = find_col(f, ["combo", "combination", "drugs", "pair", "name"], what="faers combo")
    ror_c = find_col(f, ["ror", "reporting_odds_ratio"], required=False, what="ror")
    sig_c = find_col(f, ["signal", "is_signal", "ror_signal", "flag"], required=False, what="signal")
    f["key"] = f[combo_c].map(norm_combo)
    if sig_c:
        f["signal"] = (f[sig_c].astype(str).str.lower()
                       .isin(["1", "true", "yes", "signal", "y"])).astype(int)
    elif ror_c:
        f["signal"] = (f[ror_c] >= 2.0).astype(int)
        print("NOTE: no signal column; using ROR>=2 proxy.")
    else:
        sys.exit("FAERS file has neither signal nor ROR column.")
    f["ror"] = f[ror_c] if ror_c else np.nan
    f["is_stim"] = f["key"].apply(lambda k: bool(set(k.split("+")) & STIMULANTS))

    # parse combos -> drug dicts once
    f["drugs"] = f["key"].apply(lambda k: key_to_drugdict(k, args.dose))
    bad = f[f["drugs"].isna()]
    if len(bad):
        print(f"WARNING: {len(bad)} combos did not map to model drugs, skipping:",
              list(bad["key"]))
    f = f[f["drugs"].notna()].reset_index(drop=True)
    print(f"Running {len(f)} combos x {len(factors)} IKs factors at {args.beats} beats.\n")

    summary_rows = []
    detail_rows = []
    for fi, factor in enumerate(factors, 1):
        apply_iks_factor(factor)
        base_qtc, _ = run_qtc_ikr(None, args.params, args.beats)
        print(f"[{fi}/{len(factors)}] IKs={factor:>4}  baseline QTc={base_qtc:.1f} ms")

        dqtc, ikr_vals = [], []
        for i, row in f.iterrows():
            qtc, ikr = run_qtc_ikr(row["drugs"], args.params, args.beats)
            dq = qtc - base_qtc
            dqtc.append(dq)
            ikr_vals.append(ikr)
            detail_rows.append({"iks_factor": factor, "key": row["key"],
                                "dqtc": dq, "ikr": ikr, "ror": row["ror"],
                                "signal": row["signal"], "is_stim": row["is_stim"]})
        f["_dq"] = dqtc
        f["_ikr"] = ikr_vals

        def block(sub):
            if len(sub) < 3 or sub["signal"].nunique() < 2:
                return dict(auc_dq=np.nan, auc_ikr=np.nan, rho=np.nan)
            return dict(
                auc_dq=auc_mw(sub["_dq"], sub["signal"]),
                auc_ikr=auc_mw(sub["_ikr"], sub["signal"]) if sub["_ikr"].notna().all() else np.nan,
                rho=spearmanr(sub["_dq"], sub["ror"], nan_policy="omit")[0]
                if sub["ror"].notna().sum() >= 3 else np.nan,
            )

        allb = block(f)
        stb = block(f[f.is_stim])
        nsb = block(f[~f.is_stim])
        tracked = {k: round(float(f.loc[f.key == k, "_dq"].iloc[0]), 1)
                   for k in TRACKED if (f.key == k).any()}
        print(f"        AUC dQTc all={allb['auc_dq']:.3f} stim={stb['auc_dq']:.3f} "
              f"nonstim={nsb['auc_dq']:.3f} | stim rho={stb['rho']:+.3f} | "
              f"MPH+AMP={tracked.get('AMP+MPH','?')} ARI+NOR={tracked.get('ARI+NOR','?')}")

        summary_rows.append({
            "iks_factor": factor, "baseline_qtc": base_qtc,
            "auc_dq_all": allb["auc_dq"], "auc_dq_stim": stb["auc_dq"],
            "auc_dq_nonstim": nsb["auc_dq"], "auc_ikr_all": allb["auc_ikr"],
            "rho_all": allb["rho"], "rho_stim": stb["rho"], "rho_nonstim": nsb["rho"],
            **{f"dq_{k}": tracked.get(k, np.nan) for k in TRACKED},
        })

    summ = pd.DataFrame(summary_rows)
    summ.to_csv(args.out, index=False)
    pd.DataFrame(detail_rows).to_csv(args.out.replace(".csv", "_detail.csv"), index=False)
    print(f"\nSaved -> {args.out}")

    # ---- figure ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))

        for k in TRACKED:
            col = f"dq_{k}"
            if col in summ:
                ax[0].plot(summ.iks_factor, summ[col], marker="o", label=k)
        ax[0].axhline(0, color="gray", lw=0.8, ls="--")
        ax[0].set(title="Tracked combo delta-QTc", xlabel="IKs factor", ylabel="delta-QTc (ms)")
        ax[0].legend(fontsize=8)

        ax[1].plot(summ.iks_factor, summ.auc_dq_all, marker="o", label="all")
        ax[1].plot(summ.iks_factor, summ.auc_dq_stim, marker="o", label="stimulant")
        ax[1].plot(summ.iks_factor, summ.auc_dq_nonstim, marker="o", label="non-stim")
        ax[1].axhline(0.5, color="gray", lw=0.8, ls="--")
        ax[1].set(title="FAERS concordance AUC (delta-QTc)", xlabel="IKs factor",
                  ylabel="AUC", ylim=(0.4, 0.8))
        ax[1].legend(fontsize=8)

        ax[2].plot(summ.iks_factor, summ.rho_stim, marker="o", label="stimulant")
        ax[2].plot(summ.iks_factor, summ.rho_all, marker="o", label="all")
        ax[2].axhline(0, color="gray", lw=0.8, ls="--")
        ax[2].set(title="Spearman delta-QTc vs ROR", xlabel="IKs factor", ylabel="rho")
        ax[2].legend(fontsize=8)

        plt.tight_layout()
        plt.savefig(args.fig, dpi=140)
        print(f"Saved -> {args.fig}")
    except Exception as e:
        print(f"(figure skipped: {e})")

    print("\nRead the figure like this:")
    print("  - Panel 1: if MPH+AMP slides smoothly to ~0 while ARI combos stay flat,")
    print("    the risk inversion is a gradual, robust function of IKs, not a cliff.")
    print("  - Panel 2/3: if AUC and stim-rho stay roughly flat across the range,")
    print("    the 0.632 / +0.503 concordance is stable and safe to build the paper on.")
    print("  - Sanity check: at your CURRENT IKs factor the AUC should read ~0.632")
    print("    and stim rho ~+0.503. If not, the script isn't hitting your IKs term.")


if __name__ == "__main__":
    main()
