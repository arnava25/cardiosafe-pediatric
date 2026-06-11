#!/usr/bin/env python3
"""
CardioSafe Pediatric — supratherapeutic exposure sweep.

Where the mechanistic model earns its keep: extrapolating genuine repolarization
risk into exposures you can't get clinical data for (overdose, CYP2D6 poor
metabolizers, PK drug-drug interactions). For each hERG-relevant combination it
scales free Cmax 1x -> 10x and reports the GENUINE metric (ΔAPD90, not Bazett
ΔQTc), then computes the safety margin: the fold-over-therapeutic exposure at
which genuine APD90 prolongation crosses a clinically concerning line.

NOTE: only the hERG pathway is concentration-dependent in the current model
(Hill block scales with free Cmax). The autonomic effects are presence-based and
do NOT scale, so this sweep specifically probes hERG-mediated risk under high
exposure -- which is the TdP-relevant, dose-dependent mechanism.

Usage (tmux; ~2-3 h with --grid filter at 150 beats):
    python3 src/supratherapeutic_sweep.py \
        --grid results/risk_grid_results.csv \
        --params data/herg_master_params.csv \
        --multipliers 1,2,3,5,10 \
        --beats 150 \
        --threshold 10
"""
import argparse, sys
import numpy as np, pandas as pd
sys.path.insert(0, "src")
from ord_model import run_simulation, load_drug_params

def parse_drugs(row):
    a, b = row.get("drug_A"), row.get("drug_B")
    if isinstance(a, str) and a and isinstance(b, str) and b:
        return [a, b]
    s = str(row.get("drugs", ""))
    return [d.strip() for d in s.split("+") if d.strip()]

def fold_to(mults, dapd, thr):
    for i in range(len(mults)):
        if dapd[i] >= thr:
            if i == 0:
                return mults[0]
            x0, x1, y0, y1 = mults[i-1], mults[i], dapd[i-1], dapd[i]
            return x0 + (thr - y0) * (x1 - x0) / (y1 - y0) if y1 != y0 else x1
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", required=True, help="risk_grid_results.csv (combo source)")
    ap.add_argument("--params", default="data/herg_master_params.csv")
    ap.add_argument("--multipliers", default="1,2,3,5,10")
    ap.add_argument("--beats", type=int, default=150)
    ap.add_argument("--min-ikr", type=float, default=0.5,
                    help="only sweep combos with therapeutic IKr_block_pct >= this (others don't scale)")
    ap.add_argument("--threshold", type=float, default=10.0,
                    help="genuine ΔAPD90 (ms) considered clinically concerning")
    ap.add_argument("--out", default="results/supratherapeutic_sweep.csv")
    ap.add_argument("--fig", default="docs/figures/figureS_supratherapeutic.png")
    a = ap.parse_args()

    mults = [float(x) for x in a.multipliers.split(",")]
    params = load_drug_params(a.params)
    grid = pd.read_csv(a.grid)

    # focus on combos whose risk is concentration-dependent (have hERG block)
    grid = grid[grid["IKr_block_pct"].fillna(0) >= a.min_ikr].copy()
    print(f"Sweeping {len(grid)} hERG-relevant combos x {len(mults)} exposures "
          f"at {a.beats} beats. Concern line: genuine ΔAPD90 >= {a.threshold} ms.\n")

    base = run_simulation(None, a.params, n_beats=a.beats, verbose=False)
    bAPD = base["APD90"]
    print(f"Baseline APD90 = {bAPD:.1f} ms\n")

    rows = []
    for _, r in grid.iterrows():
        drugs = parse_drugs(r)
        if any(d not in params for d in drugs):
            print(f"  skip {r['combination']} (drug not in params)"); continue
        dapd = []
        for mlt in mults:
            combo = {d: float(params[d]["cmax_free_nM"]) * mlt for d in drugs}
            res = run_simulation(combo, a.params, n_beats=a.beats, verbose=False)
            dapd.append(res["APD90"] - bAPD)
        ft = fold_to(mults, dapd, a.threshold)
        rec = {"combination": r["combination"],
               **{f"dAPD_{int(m)}x": round(dapd[i], 1) for i, m in enumerate(mults)},
               "max_dAPD": round(max(dapd), 1),
               "fold_to_concern": (round(ft, 1) if ft is not None else np.nan)}
        rows.append(rec)
        ftxt = f"{ft:.1f}x" if ft is not None else f">{int(max(mults))}x"
        print(f"  {r['combination']:16s} 1x={dapd[0]:+5.1f} -> {int(max(mults))}x={dapd[-1]:+6.1f} ms"
              f"  | margin: {ftxt}")

    df = pd.DataFrame(rows).sort_values(
        "fold_to_concern", ascending=True, na_position="last")
    df.to_csv(a.out, index=False)
    print(f"\nSaved -> {a.out}")

    print("\n── NARROWEST SAFETY MARGINS (genuine ΔAPD90 reaches concern soonest) ──")
    for _, r in df.head(10).iterrows():
        m = r["fold_to_concern"]
        mtxt = f"{m:.1f}x over therapeutic" if pd.notna(m) else f">{int(max(mults))}x (safe in range)"
        print(f"  {r['combination']:16s} concern at {mtxt}   (max {r['max_dAPD']:+.1f} ms)")

    # figure
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        top = df.head(8)["combination"].tolist()
        plt.figure(figsize=(8, 5))
        for _, r in df[df.combination.isin(top)].iterrows():
            ys = [r[f"dAPD_{int(m)}x"] for m in mults]
            plt.plot(mults, ys, marker="o", label=r["combination"])
        plt.axhline(a.threshold, color="red", ls="--", lw=1, label=f"concern ({a.threshold:.0f} ms)")
        plt.xlabel("free Cmax (x therapeutic)"); plt.ylabel("genuine ΔAPD90 (ms)")
        plt.title("Supratherapeutic genuine repolarization risk")
        plt.legend(fontsize=8); plt.tight_layout(); plt.savefig(a.fig, dpi=140)
        print(f"Saved -> {a.fig}")
    except Exception as e:
        print(f"(figure skipped: {e})")

    print("\nRead it: combos with a NARROW margin (concern at 2-3x) are the ones that")
    print("turn dangerous under CYP2D6 poor-metabolism or modest overdose -- exactly the")
    print("regime FAERS and therapeutic-dose data can't see. Wide margin (>10x) = robust.")

if __name__ == "__main__":
    main()
