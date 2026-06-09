"""
sympathomimetic_sensitivity.py — CardioSafe Pediatric
======================================================
Stress-tests the sympathomimetic pathway parameters.

Base case: CL_reduction=10%, GCaL_upregulation=15%
Sweep:     CL_reduction in {5%, 7.5%, 10%, 12.5%, 15%}
           GCaL_upregulation in {10%, 12.5%, 15%, 17.5%, 20%}

For each parameter combination, reruns all stimulant-containing pairs
and reports the range of delta-QTc predictions.

Key question: does MPH+ARI stay in a clinically meaningful risk range
(>10ms MODERATE) across the full parameter space?

Usage:
    python3 src/sympathomimetic_sensitivity.py
    python3 src/sympathomimetic_sensitivity.py --beats 100   # faster
"""

import sys
import argparse
import itertools
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ord_model import run_simulation

PARAMS_PATH = str(ROOT / "data" / "herg_master_params.csv")
OUT_CSV     = ROOT / "results" / "sympathomimetic_sensitivity.csv"

# Parameter grid
CL_REDUCTIONS    = [0.05, 0.075, 0.10, 0.125, 0.15]   # fraction
GCAL_UPREGULATIONS = [0.10, 0.125, 0.15, 0.175, 0.20] # fraction

# Focus combinations — stimulant-containing pairs
FOCUS_COMBOS = [
    ("Methylphenidate", "Amphetamine"),
    ("Methylphenidate", "Aripiprazole"),
    ("Methylphenidate", "Quetiapine"),
    ("Methylphenidate", "Sertraline"),
    ("Methylphenidate", "Nortriptyline"),
    ("Amphetamine",     "Aripiprazole"),
    ("Amphetamine",     "Sertraline"),
]

DRUG_CODES = {
    "Methylphenidate":"MPH","Amphetamine":"AMP","Risperidone":"RIS",
    "Quetiapine":"QUE","Aripiprazole":"ARI","Sertraline":"SER",
    "Fluoxetine":"FLU","Escitalopram":"ESC","Clonidine":"CLO",
    "Guanfacine":"GUA","Imipramine":"IMI","Nortriptyline":"NOR",
}

def tier(dq):
    if dq >= 20:  return "HIGH"
    if dq >= 10:  return "MODERATE"
    if dq >= 5:   return "LOW-MOD"
    if dq >= 0:   return "LOW"
    return "PROTECTIVE"

def run_with_params(drug_a, drug_b, cl_red, gcal_up, n_beats, baseline_qtc):
    """Run simulation with custom sympathomimetic parameters."""
    # Patch the sympathomimetic parameters into the combo dict
    # ord_model reads cl_reduction and gcal_mult from combo metadata
    combo = {
        drug_a: "therapeutic",
        drug_b: "therapeutic",
        "__cl_reduction__":  cl_red,
        "__gcal_mult__":     1.0 + gcal_up,
    }
    try:
        res = run_simulation(combo, PARAMS_PATH, n_beats=n_beats, verbose=False)
        return round(res["QTc"] - baseline_qtc, 1)
    except Exception as e:
        return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--beats", type=int, default=200)
    args = parser.parse_args()
    N_BEATS = args.beats

    print("CardioSafe Pediatric — Sympathomimetic Parameter Sensitivity")
    print(f"CL reductions:     {[f'{x*100:.1f}%' for x in CL_REDUCTIONS]}")
    print(f"GCaL upregulations: {[f'{x*100:.1f}%' for x in GCAL_UPREGULATIONS]}")
    print(f"Combos: {len(FOCUS_COMBOS)}  |  Grid: {len(CL_REDUCTIONS)*len(GCAL_UPREGULATIONS)} parameter sets")
    print(f"Total simulations: {len(FOCUS_COMBOS)*len(CL_REDUCTIONS)*len(GCAL_UPREGULATIONS)}")
    print()

    # Check if ord_model supports __cl_reduction__ override
    # If not, patch it directly
    import ord_model
    src = open(ROOT / "src" / "ord_model.py").read()
    has_override = "__cl_reduction__" in src

    if not has_override:
        print("NOTE: ord_model does not support parameter override via combo dict.")
        print("Patching sympathomimetic modifier function directly for each run.")
        print()

    # Baseline
    print("Running baseline...")
    baseline = run_simulation(None, PARAMS_PATH, n_beats=N_BEATS, verbose=True)
    baseline_qtc = baseline["QTc"]
    print(f"Baseline QTc: {baseline_qtc:.1f} ms\n")

    results = []
    total = len(FOCUS_COMBOS) * len(CL_REDUCTIONS) * len(GCAL_UPREGULATIONS)
    done = 0

    # Since ord_model uses its own internal params, we'll monkey-patch
    # compute_autonomic_modifiers to inject custom CL/GCaL for stimulants
    import ord_model as om

    original_autonomic = om.compute_autonomic_modifiers

    for cl_red, gcal_up in itertools.product(CL_REDUCTIONS, GCAL_UPREGULATIONS):

        # Monkey-patch sympathomimetic modifiers
        def make_modifier(cl_r, gcal_u):
            def compute_autonomic_modifiers_patched(drug_combo, drug_params):
                hr_mult = 1.0; gcal_mult = 1.0; gna_mult = 1.0
                for drug in drug_combo:
                    if drug.startswith("__"):
                        continue
                    if drug not in drug_params:
                        continue
                    dp = drug_params[drug]
                    mech = dp.get("primary_mechanism", "")
                    if mech == "sympathomimetic":
                        hr_mult   *= (1.0 - cl_r)
                        gcal_mult *= (1.0 + gcal_u)
                    elif mech == "autonomic":
                        hr_mult  *= 1.10
                        gna_mult *= 0.95
                return hr_mult, gcal_mult, gna_mult
            return compute_autonomic_modifiers_patched

        om.compute_autonomic_modifiers = make_modifier(cl_red, gcal_up)

        for drug_a, drug_b in FOCUS_COMBOS:
            combo = {drug_a: "therapeutic", drug_b: "therapeutic"}
            try:
                res = om.run_simulation(combo, PARAMS_PATH,
                                        n_beats=N_BEATS, verbose=False)
                dqtc = round(res["QTc"] - baseline_qtc, 1)
            except Exception as e:
                dqtc = None

            done += 1
            a = DRUG_CODES[drug_a]; b = DRUG_CODES[drug_b]
            is_base = (abs(cl_red - 0.10) < 0.001 and abs(gcal_up - 0.15) < 0.001)
            marker = " ← BASE" if is_base else ""
            if dqtc is not None:
                print(f"  [{done:3d}/{total}] {a}+{b} | "
                      f"CL-{cl_red*100:.1f}% GCaL+{gcal_up*100:.1f}%: "
                      f"dQTc={dqtc:+.1f}ms ({tier(dqtc)}){marker}")

            results.append({
                "combination":    f"{a}+{b}",
                "cl_reduction":   cl_red,
                "gcal_up":        gcal_up,
                "dQTc_ms":        dqtc,
                "risk_tier":      tier(dqtc) if dqtc is not None else None,
                "is_base_case":   is_base,
            })

    # Restore
    om.compute_autonomic_modifiers = original_autonomic

    df = pd.DataFrame(results).dropna(subset=["dQTc_ms"])
    df.to_csv(OUT_CSV, index=False)
    print(f"\nSaved -> {OUT_CSV}")

    # ── SUMMARY ──────────────────────────────────────────────────────────────
    print("\n── RANGE SUMMARY BY COMBINATION ──")
    print(f"{'Combo':12s}  {'Min dQTc':10s}  {'Max dQTc':10s}  "
          f"{'Base dQTc':10s}  {'All MODERATE+?':15s}  Tier range")
    print("-" * 80)

    for combo in [f"{DRUG_CODES[a]}+{DRUG_CODES[b]}" for a,b in FOCUS_COMBOS]:
        sub = df[df["combination"] == combo]
        if sub.empty:
            continue
        mn   = sub["dQTc_ms"].min()
        mx   = sub["dQTc_ms"].max()
        base = sub[sub["is_base_case"]]["dQTc_ms"].values
        base_val = base[0] if len(base) > 0 else None
        all_mod  = all(sub["dQTc_ms"] >= 10)
        tiers    = sorted(sub["risk_tier"].unique())
        print(f"{combo:12s}  {mn:+10.1f}  {mx:+10.1f}  "
              f"{str(base_val) if base_val else 'N/A':10s}  "
              f"{'YES' if all_mod else 'NO':15s}  {'/'.join(tiers)}")

    # ── FIGURE ───────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 4, figsize=(15, 8))
    axes = axes.flatten()

    combos = [f"{DRUG_CODES[a]}+{DRUG_CODES[b]}" for a, b in FOCUS_COMBOS]
    cmap   = plt.cm.RdYlGn_r

    for idx, combo in enumerate(combos):
        if idx >= len(axes):
            break
        ax   = axes[idx]
        sub  = df[df["combination"] == combo]
        if sub.empty:
            ax.set_visible(False)
            continue

        # Build grid
        cl_vals   = sorted(sub["cl_reduction"].unique())
        gcal_vals = sorted(sub["gcal_up"].unique())
        grid = np.zeros((len(gcal_vals), len(cl_vals)))
        for i, g in enumerate(gcal_vals):
            for j, c in enumerate(cl_vals):
                v = sub[(sub["cl_reduction"]==c) & (sub["gcal_up"]==g)]["dQTc_ms"]
                grid[i,j] = v.values[0] if len(v) > 0 else np.nan

        im = ax.imshow(grid, cmap=cmap, vmin=0, vmax=25, aspect="auto",
                       origin="lower")
        ax.set_xticks(range(len(cl_vals)))
        ax.set_yticks(range(len(gcal_vals)))
        ax.set_xticklabels([f"{x*100:.0f}%" for x in cl_vals], fontsize=7)
        ax.set_yticklabels([f"{x*100:.0f}%" for x in gcal_vals], fontsize=7)
        ax.set_xlabel("CL reduction", fontsize=8)
        ax.set_ylabel("GCaL upregulation", fontsize=8)
        ax.set_title(combo, fontsize=10, fontweight="bold", color="#1B3A6B")

        # Mark base case
        bi = gcal_vals.index(0.15) if 0.15 in gcal_vals else 2
        bj = cl_vals.index(0.10)   if 0.10 in cl_vals   else 2
        ax.add_patch(plt.Rectangle((bj-0.5, bi-0.5), 1, 1,
            fill=False, edgecolor="#1B3A6B", linewidth=2.5, zorder=5))

        # Annotate values
        for i in range(len(gcal_vals)):
            for j in range(len(cl_vals)):
                v = grid[i, j]
                color = "white" if v > 15 else "#1F2937"
                ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                        fontsize=7, color=color, fontweight="bold")

        # MODERATE threshold line
        ax.axhline(1.5, color="white", linewidth=0.8, linestyle="--", alpha=0.6)

    # Hide unused subplots
    for idx in range(len(combos), len(axes)):
        axes[idx].set_visible(False)

    plt.suptitle(
        "Figure 5. Sympathomimetic Parameter Sensitivity — ΔQTc across CL reduction and GCaL upregulation\n"
        "Blue box = base case (CL-10%, GCaL+15%). Color scale: green=low risk, red=high risk.",
        fontsize=10, fontweight="bold", color="#1B3A6B", y=1.01
    )
    plt.colorbar(im, ax=axes[:len(combos)], label="ΔQTc (ms)",
                 fraction=0.015, pad=0.04)
    plt.tight_layout()
    plt.savefig("/mnt/user-data/outputs/figure5_sympathomimetic_sensitivity.png",
                dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print("Saved -> figure5_sympathomimetic_sensitivity.png")

if __name__ == "__main__":
    main()
