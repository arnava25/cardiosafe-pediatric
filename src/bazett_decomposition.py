"""
bazett_decomposition.py — CardioSafe Pediatric
================================================
Decomposes model delta-QTc for stimulant combinations into:
  (1) Bazett artifact from CL reduction alone (no GCaL change)
  (2) Genuine ICaL-mediated APD90 change from GCaL upregulation
  (3) Full model (current: CL + GCaL together)

Runs each stimulant combination three ways at N_BEATS.
Output: table showing Bazett artifact vs genuine ICaL component.

Usage:
    python3 src/bazett_decomposition.py
    python3 src/bazett_decomposition.py --beats 100
"""

import sys, argparse
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import ord_model as om

PARAMS_PATH = str(ROOT / "data" / "herg_master_params.csv")
OUT_CSV     = ROOT / "results" / "bazett_decomposition.csv"
OUT_FIG     = ROOT / "docs" / "figures" / "figureS_bazett.png"

STIM_COMBOS = [
    ("Methylphenidate","Amphetamine"),
    ("Methylphenidate","Aripiprazole"),
    ("Methylphenidate","Quetiapine"),
    ("Methylphenidate","Sertraline"),
    ("Methylphenidate","Nortriptyline"),
    ("Methylphenidate","Risperidone"),
    ("Methylphenidate","Fluoxetine"),
    ("Methylphenidate","Imipramine"),
    ("Methylphenidate","Escitalopram"),
    ("Methylphenidate","Clonidine"),
    ("Methylphenidate","Guanfacine"),
    ("Amphetamine","Aripiprazole"),
    ("Amphetamine","Sertraline"),
    ("Amphetamine","Quetiapine"),
    ("Amphetamine","Nortriptyline"),
    ("Amphetamine","Risperidone"),
    ("Amphetamine","Fluoxetine"),
]

DRUG_CODES = {
    "Methylphenidate":"MPH","Amphetamine":"AMP","Risperidone":"RIS",
    "Quetiapine":"QUE","Aripiprazole":"ARI","Sertraline":"SER",
    "Fluoxetine":"FLU","Escitalopram":"ESC","Clonidine":"CLO",
    "Guanfacine":"GUA","Imipramine":"IMI","Nortriptyline":"NOR",
}

CL_REDUCTION     = 0.10   # 10% — current model
GCAL_UPREGULATION= 0.15   # 15% — current model

original_autonomic = om.compute_autonomic_modifiers

def make_modifier(use_cl=True, use_gcal=True):
    def patched(drug_combo, drug_params):
        hr_mult=1.0; gcal_mult=1.0; gna_mult=1.0
        for drug in drug_combo:
            if drug.startswith("__"): continue
            if drug not in drug_params: continue
            dp = drug_params[drug]
            mech = dp.get("primary_mechanism","")
            if mech == "sympathomimetic":
                if use_cl:   hr_mult   *= (1.0 - CL_REDUCTION)
                if use_gcal: gcal_mult *= (1.0 + GCAL_UPREGULATION)
            elif mech == "autonomic":
                hr_mult  *= 1.10
                gna_mult *= 0.95
        return hr_mult, gcal_mult, gna_mult
    return patched

def tier(dq):
    if dq >= 20: return "HIGH"
    if dq >= 10: return "MODERATE"
    if dq >= 5:  return "LOW-MOD"
    if dq >= 0:  return "LOW"
    return "PROTECTIVE"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--beats", type=int, default=200)
    args = parser.parse_args()
    N = args.beats

    print("CardioSafe Pediatric — Bazett Decomposition Analysis")
    print(f"CL_reduction={CL_REDUCTION*100:.0f}%  GCaL_upregulation={GCAL_UPREGULATION*100:.0f}%")
    print(f"Beats: {N} | Combinations: {len(STIM_COMBOS)}\n")

    # Baseline
    om.compute_autonomic_modifiers = original_autonomic
    baseline = om.run_simulation(None, PARAMS_PATH, n_beats=N, verbose=False)
    baseline_qtc = baseline["QTc"]
    print(f"Baseline QTc: {baseline_qtc:.1f} ms\n")

    results = []
    total = len(STIM_COMBOS) * 3
    done = 0

    for drug_a, drug_b in STIM_COMBOS:
        a = DRUG_CODES[drug_a]; b = DRUG_CODES[drug_b]
        combo = {drug_a: "therapeutic", drug_b: "therapeutic"}

        # Run A: CL only (Bazett artifact, no GCaL change)
        om.compute_autonomic_modifiers = make_modifier(use_cl=True, use_gcal=False)
        r_cl = om.run_simulation(combo, PARAMS_PATH, n_beats=N, verbose=False)
        dq_cl = round(r_cl["QTc"] - baseline_qtc, 1)
        done += 1

        # Run B: GCaL only (no CL change)
        om.compute_autonomic_modifiers = make_modifier(use_cl=False, use_gcal=True)
        r_gcal = om.run_simulation(combo, PARAMS_PATH, n_beats=N, verbose=False)
        dq_gcal = round(r_gcal["QTc"] - baseline_qtc, 1)
        done += 1

        # Run C: Full model (current)
        om.compute_autonomic_modifiers = make_modifier(use_cl=True, use_gcal=True)
        r_full = om.run_simulation(combo, PARAMS_PATH, n_beats=N, verbose=False)
        dq_full = round(r_full["QTc"] - baseline_qtc, 1)
        ikr = round(r_full.get("IKr_block_pct", 0), 2)
        done += 1

        # Genuine ICaL component = full - CL_only
        genuine_ical = round(dq_full - dq_cl, 1)
        bazett_frac  = round(dq_cl  / dq_full * 100, 1) if dq_full > 0 else 0.0
        ical_frac    = round(genuine_ical / dq_full * 100, 1) if dq_full > 0 else 0.0

        print(f"  {a}+{b}:  CL-only={dq_cl:+.1f}ms  GCaL-only={dq_gcal:+.1f}ms  "
              f"Full={dq_full:+.1f}ms  Genuine-ICaL={genuine_ical:+.1f}ms  "
              f"Bazett%={bazett_frac:.0f}%  IKr={ikr:.2f}%")

        results.append({
            "combination":         f"{a}+{b}",
            "dQTc_full_ms":        dq_full,
            "tier_full":           tier(dq_full),
            "dQTc_CL_only_ms":     dq_cl,
            "dQTc_GCaL_only_ms":   dq_gcal,
            "genuine_ICaL_ms":     genuine_ical,
            "bazett_artifact_pct": bazett_frac,
            "genuine_ICaL_pct":    ical_frac,
            "IKr_block_pct":       ikr,
        })

    om.compute_autonomic_modifiers = original_autonomic

    df = pd.DataFrame(results).sort_values("dQTc_full_ms", ascending=False)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nSaved -> {OUT_CSV}")

    # ── Summary ───────────────────────────────────────────────────────────────
    print("\n── SUMMARY ──")
    print(f"Mean Bazett artifact:     {df['dQTc_CL_only_ms'].mean():.1f} ms")
    print(f"Mean genuine ICaL:        {df['genuine_ICaL_ms'].mean():.1f} ms")
    print(f"Mean full model dQTc:     {df['dQTc_full_ms'].mean():.1f} ms")
    print(f"Mean Bazett fraction:     {df['bazett_artifact_pct'].mean():.1f}%")
    print(f"Mean genuine ICaL frac:   {df['genuine_ICaL_pct'].mean():.1f}%")

    # ── Figure ────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    NAVY="#1B3A6B"; AMBER="#D97706"; TEAL="#0D5C6B"; MGRAY="#CDD3DE"

    # Stacked bar: Bazett artifact vs genuine ICaL
    ax = axes[0]
    combos = df["combination"].tolist()
    x = range(len(combos))
    ax.bar(x, df["dQTc_CL_only_ms"], label="Bazett artifact (CL reduction)",
           color=AMBER, alpha=0.85, width=0.6)
    ax.bar(x, df["genuine_ICaL_ms"],
           bottom=df["dQTc_CL_only_ms"],
           label="Genuine ICaL component (GCaL upregulation)",
           color=NAVY, alpha=0.85, width=0.6)
    ax.axhline(10, color="#B01C1C", lw=1, linestyle="--", alpha=0.7,
               label="MODERATE threshold (10 ms)")
    ax.set_xticks(list(x))
    ax.set_xticklabels(combos, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("ΔQTc contribution (ms)", fontsize=10)
    ax.set_title("A.  Bazett Artifact vs Genuine ICaL Component", fontsize=11,
                 fontweight="bold", color=NAVY, loc="left")
    ax.legend(fontsize=8, framealpha=0.9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.yaxis.grid(True, color=MGRAY, linewidth=0.5)
    ax.set_axisbelow(True)

    # Scatter: Bazett artifact vs genuine ICaL per combination
    ax2 = axes[1]
    ax2.scatter(df["dQTc_CL_only_ms"], df["genuine_ICaL_ms"],
                color=NAVY, s=60, alpha=0.85, zorder=3)
    for _, r in df.iterrows():
        ax2.annotate(r["combination"], (r["dQTc_CL_only_ms"], r["genuine_ICaL_ms"]),
                     fontsize=7, xytext=(3, 2), textcoords="offset points",
                     color=TEAL)
    ax2.set_xlabel("Bazett artifact from CL reduction (ms)", fontsize=10)
    ax2.set_ylabel("Genuine ICaL component from GCaL upregulation (ms)", fontsize=10)
    ax2.set_title("B.  Per-combination decomposition", fontsize=11,
                  fontweight="bold", color=NAVY, loc="left")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.yaxis.grid(True, color=MGRAY, linewidth=0.5)
    ax2.set_axisbelow(True)

    plt.suptitle(
        "Supplementary Figure S_Bazett. Bazett Correction Decomposition for Stimulant Combinations\n"
        "Total ΔQTc (full model) decomposes into Bazett artifact from CL reduction "
        "and genuine ICaL-mediated APD90 change from GCaL upregulation.",
        fontsize=9, y=1.01, color=NAVY, fontweight="bold")
    plt.tight_layout()
    plt.savefig(str(OUT_FIG), dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved -> {OUT_FIG}")

if __name__ == "__main__":
    main()
