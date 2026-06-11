#!/usr/bin/env python3
"""
CardioSafe Pediatric — rate-correction comparison (Bazett vs Fridericia vs raw APD90).

Pure post-processing: recomputes both heart-rate corrections from the grid's
APD90 and CL_eff columns (no simulation). Shows that rate correction distorts
QTc risk for any psychotropic that moves heart rate -- INFLATING stimulant
combos (HR up) and MASKING alpha-2 agonist combos (HR down) -- while the
uncorrected interval (ΔAPD90) does not.

NOTE: ΔAPD90 is the uncorrected interval and still carries physiological
rate-dependence (restitution); isolating pure channel pharmacology from
autonomic rate needs a fixed-CL re-run on the engine. This script characterizes
the *correction artifact*, which it does without any simulation.

Usage:
    python3 src/rate_correction.py --grid results/risk_grid_results.csv
"""
import argparse
import numpy as np, pandas as pd

STIM = ["MPH", "AMP"]
A2 = ["CLO", "GUA"]

def classify(combo):
    parts = set(str(combo).replace("+", " ").split())
    if parts & set(STIM): return "stimulant"
    if parts & set(A2):   return "alpha-2"
    return "hERG/other"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", required=True)
    ap.add_argument("--out", default="results/rate_correction_comparison.csv")
    ap.add_argument("--fig", default="docs/figures/figure_rate_correction.png")
    a = ap.parse_args()

    df = pd.read_csv(a.grid)
    # column resolution (tolerant of the Δ-prefixed unicode names)
    def col(*cands):
        for c in cands:
            for x in df.columns:
                if c.lower() in x.lower(): return x
        raise SystemExit(f"missing column among {cands}; have {list(df.columns)}")
    APD = col("APD90_ms", "apd90"); CL = col("CL_eff", "cl_eff")
    dAPDc = col("ΔAPD90", "dapd90", "delta_apd")
    df = df.dropna(subset=[dAPDc]).copy()

    # baseline APD90 (CL=1000): back out from any row
    base = float((df[APD] - df[dAPDc]).median())
    RR = df[CL] / 1000.0
    df["QTcB"] = df[APD] / np.sqrt(RR)
    df["QTcF"] = df[APD] / RR ** (1/3)
    df["dQTcB"] = df["QTcB"] - base
    df["dQTcF"] = df["QTcF"] - base
    df["dAPD"] = df[dAPDc]
    df["class"] = df["combination"].map(classify)

    out = df[["combination", "class", "dAPD", "dQTcF", "dQTcB", APD, CL]].copy()
    out = out.rename(columns={APD: "APD90_ms", CL: "CL_eff_ms"})
    out["_distortion"] = (out["dQTcB"] - out["dAPD"]).abs()
    out = out.sort_values("_distortion", ascending=False).drop(columns="_distortion")
    out.to_csv(a.out, index=False)
    print(f"baseline APD90 = {base:.1f} ms   (Δ measured vs this)\n")

    print("Most distorted by rate correction:")
    print(f"  {'combo':<15}{'class':<11}{'dAPD':>7}{'QTcF':>7}{'QTcB':>7}")
    for _, r in out.head(12).iterrows():
        print(f"  {r.combination:<15}{r['class']:<11}{r.dAPD:>+7.1f}{r.dQTcF:>+7.1f}{r.dQTcB:>+7.1f}")

    print("\nClass means (ms):")
    for cl, sub in df.groupby("class"):
        print(f"  {cl:<12} dAPD={sub.dAPD.mean():+5.1f}  QTcF={sub.dQTcF.mean():+5.1f}  "
              f"QTcB={sub.dQTcB.mean():+5.1f}   (n={len(sub)})")
    print(f"\nSaved -> {a.out}")

    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        cmap = {"stimulant": "#d1495b", "alpha-2": "#3a86ff", "hERG/other": "#8d99ae"}
        fig, ax = plt.subplots(figsize=(7, 6))
        lim = [df[["dAPD", "dQTcB"]].min().min() - 2, df[["dAPD", "dQTcB"]].max().max() + 2]
        ax.plot(lim, lim, "k--", lw=1, label="no distortion (QTc = ΔAPD90)")
        for cl, sub in df.groupby("class"):
            ax.scatter(sub.dAPD, sub.dQTcB, c=cmap[cl], label=f"{cl} (Bazett)", s=42, edgecolor="white")
            ax.scatter(sub.dAPD, sub.dQTcF, c=cmap[cl], marker="^", s=34, alpha=0.55)
        ax.set_xlabel("uncorrected ΔAPD90 (ms)"); ax.set_ylabel("rate-corrected ΔQTc (ms)")
        ax.set_title("Rate correction inflates stimulants, masks alpha-2 agonists\n"
                     "(circles = Bazett, triangles = Fridericia)")
        ax.legend(fontsize=8); ax.set_xlim(lim); ax.set_ylim(lim)
        plt.tight_layout(); plt.savefig(a.fig, dpi=140)
        print(f"Saved -> {a.fig}")
    except Exception as e:
        print(f"(figure skipped: {e})")

if __name__ == "__main__":
    main()
