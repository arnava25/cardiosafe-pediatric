"""
pk_pediatric_sensitivity.py — CardioSafe Pediatric
====================================================
Sensitivity analysis comparing adult-reference vs pediatric-specific
free Cmax values for all 12 drugs across the top combinations.

Adult reference: current herg_master_params.csv values
Pediatric values: derived from published pediatric PK literature
  RIS: Boellner/Aman Clin Ther 2007
  SER: Alderman 1998 JAACAP / FDA Zoloft label
  GUA: FDA Intuniv label children 6-12
  MPH: IR 0.5mg/kg/day estimate from OROS/TDM literature
  AMP: McGough 2003 Adderall XR scaled to 10mg/day
  QUE: Winter 2008 scaled to 150mg/day
  FLU: Wilens 2002 / FDA Prozac label avg SS
  ESC: FDA Lexapro label adolescent single dose
  IMI: Sallee 1986 estimated parent at 3mg/kg/day
  ARI: Adult benchmark (no clean pediatric data; retained)
  CLO: No pediatric data (retained)
  NOR: Therapeutic window midpoint (retained)
"""

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ord_model import run_simulation

PARAMS_PATH = str(ROOT / "data" / "herg_master_params.csv")
OUT_CSV     = ROOT / "results" / "pk_pediatric_sensitivity.csv"
OUT_FIG     = ROOT / "docs" / "figures" / "figureS3_pediatric_pk_sensitivity.png"
N_BEATS     = 200

# ── Pediatric free Cmax values (nM) ──────────────────────────────────────────
PEDIATRIC_CMAX = {
    "Methylphenidate": 21.9,   # IR 0.5mg/kg/day; MEDIUM confidence
    "Amphetamine":    148.0,   # d-AMP scaled from Adderall XR 10mg/day; MEDIUM
    "Risperidone":     3.87,   # Boellner/Aman 2007; HIGH
    "Quetiapine":     86.6,    # Winter 2008 at 150mg/day; MEDIUM
    "Aripiprazole":   11.9,    # Adult benchmark; LOW (retained)
    "Sertraline":     10.8,    # Alderman 1998 / FDA label; HIGH
    "Fluoxetine":     30.4,    # Wilens 2002 avg SS; MEDIUM
    "Escitalopram":   13.2,    # FDA Lexapro label adolescent; MEDIUM
    "Clonidine":       0.5,    # No pediatric data; LOW (retained)
    "Guanfacine":     12.2,    # FDA Intuniv label children 6-12; HIGH
    "Imipramine":     53.5,    # Sallee 1986; MEDIUM
    "Nortriptyline":  30.4,    # Therapeutic window midpoint; LOW
}

CONFIDENCE = {
    "Methylphenidate": "MEDIUM",
    "Amphetamine":     "MEDIUM",
    "Risperidone":     "HIGH",
    "Quetiapine":      "MEDIUM",
    "Aripiprazole":    "LOW",
    "Sertraline":      "HIGH",
    "Fluoxetine":      "MEDIUM",
    "Escitalopram":    "MEDIUM",
    "Clonidine":       "LOW",
    "Guanfacine":      "HIGH",
    "Imipramine":      "MEDIUM",
    "Nortriptyline":   "LOW",
}

DRUG_CODES = {
    "Methylphenidate":"MPH","Amphetamine":"AMP","Risperidone":"RIS",
    "Quetiapine":"QUE","Aripiprazole":"ARI","Sertraline":"SER",
    "Fluoxetine":"FLU","Escitalopram":"ESC","Clonidine":"CLO",
    "Guanfacine":"GUA","Imipramine":"IMI","Nortriptyline":"NOR",
}

def tier(dq):
    if dq >= 20: return "HIGH"
    if dq >= 10: return "MODERATE"
    if dq >= 5:  return "LOW-MOD"
    if dq >= 0:  return "LOW"
    return "PROTECTIVE"

# Top 30 combinations to test (by adult-reference composite score / risk)
FOCUS_COMBOS = [
    ("Methylphenidate","Aripiprazole"),
    ("Amphetamine","Aripiprazole"),
    ("Methylphenidate","Amphetamine"),
    ("Methylphenidate","Quetiapine"),
    ("Methylphenidate","Sertraline"),
    ("Amphetamine","Sertraline"),
    ("Methylphenidate","Nortriptyline"),
    ("Amphetamine","Nortriptyline"),
    ("Aripiprazole","Nortriptyline"),
    ("Aripiprazole","Imipramine"),
    ("Aripiprazole","Sertraline"),
    ("Aripiprazole","Quetiapine"),
    ("Aripiprazole","Fluoxetine"),
    ("Quetiapine","Nortriptyline"),
    ("Methylphenidate","Risperidone"),
    ("Amphetamine","Risperidone"),
    ("Methylphenidate","Imipramine"),
    ("Amphetamine","Imipramine"),
    ("Fluoxetine","Nortriptyline"),
    ("Risperidone","Nortriptyline"),
    ("Quetiapine","Aripiprazole"),
    ("Quetiapine","Fluoxetine"),
    ("Aripiprazole","Escitalopram"),
    ("Methylphenidate","Fluoxetine"),
    ("Quetiapine","Sertraline"),
    ("Methylphenidate","Escitalopram"),
    ("Risperidone","Aripiprazole"),
    ("Imipramine","Nortriptyline"),
    ("Aripiprazole","Guanfacine"),
    ("Quetiapine","Guanfacine"),
]

def main():
    print("CardioSafe Pediatric — Pediatric PK Sensitivity Analysis")
    print(f"Beats: {N_BEATS} | Combinations: {len(FOCUS_COMBOS)}")
    print()

    # Baseline
    print("Running baseline...")
    baseline = run_simulation(None, PARAMS_PATH, n_beats=N_BEATS, verbose=False)
    baseline_qtc = baseline["QTc"]
    print(f"Baseline QTc: {baseline_qtc:.1f} ms\n")

    results = []
    total = len(FOCUS_COMBOS) * 2
    done = 0

    for drug_a, drug_b in FOCUS_COMBOS:
        a = DRUG_CODES[drug_a]
        b = DRUG_CODES[drug_b]

        # ── Adult reference (therapeutic Cmax from params CSV) ────────────────
        combo_adult = {drug_a: "therapeutic", drug_b: "therapeutic"}
        res_adult = run_simulation(combo_adult, PARAMS_PATH, n_beats=N_BEATS, verbose=False)
        dqtc_adult = round(res_adult["QTc"] - baseline_qtc, 1)
        ikr_adult  = round(res_adult.get("IKr_block_pct", 0), 2)
        done += 1
        print(f"  [{done:3d}/{total}] {a}+{b} adult: {dqtc_adult:+.1f}ms ({tier(dqtc_adult)})")

        # ── Pediatric Cmax ────────────────────────────────────────────────────
        combo_peds = {
            drug_a: PEDIATRIC_CMAX[drug_a],
            drug_b: PEDIATRIC_CMAX[drug_b],
        }
        res_peds = run_simulation(combo_peds, PARAMS_PATH, n_beats=N_BEATS, verbose=False)
        dqtc_peds = round(res_peds["QTc"] - baseline_qtc, 1)
        ikr_peds  = round(res_peds.get("IKr_block_pct", 0), 2)
        done += 1
        print(f"  [{done:3d}/{total}] {a}+{b} peds:  {dqtc_peds:+.1f}ms ({tier(dqtc_peds)})")

        conf_a = CONFIDENCE[drug_a]
        conf_b = CONFIDENCE[drug_b]
        overall_conf = "HIGH" if conf_a=="HIGH" and conf_b=="HIGH" else \
                       "LOW"  if conf_a=="LOW"  or  conf_b=="LOW"  else "MEDIUM"

        results.append({
            "combination":      f"{a}+{b}",
            "drug_a":           a,
            "drug_b":           b,
            "dQTc_adult_ms":    dqtc_adult,
            "tier_adult":       tier(dqtc_adult),
            "IKr_adult_pct":    ikr_adult,
            "dQTc_peds_ms":     dqtc_peds,
            "tier_peds":        tier(dqtc_peds),
            "IKr_peds_pct":     ikr_peds,
            "delta_dQTc_ms":    round(dqtc_peds - dqtc_adult, 1),
            "tier_changed":     tier(dqtc_adult) != tier(dqtc_peds),
            "tier_escalated":   tier(dqtc_peds) in ["HIGH","MODERATE"] and
                                tier(dqtc_adult) not in ["HIGH","MODERATE"] or
                                (tier(dqtc_peds)=="HIGH" and tier(dqtc_adult)!="HIGH"),
            "conf_a":           conf_a,
            "conf_b":           conf_b,
            "overall_conf":     overall_conf,
        })

    df = pd.DataFrame(results).sort_values("dQTc_peds_ms", ascending=False)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nSaved -> {OUT_CSV}")

    # ── Summary ───────────────────────────────────────────────────────────────
    escalated = df[df["tier_escalated"]==True]
    changed   = df[df["tier_changed"]==True]

    print(f"\n── TIER ESCALATIONS ({len(escalated)}) ──")
    for _, r in escalated.iterrows():
        print(f"  {r['combination']:12s}: {r['tier_adult']} → {r['tier_peds']} "
              f"(adult {r['dQTc_adult_ms']:+.1f}ms → peds {r['dQTc_peds_ms']:+.1f}ms, "
              f"conf: {r['overall_conf']})")

    print(f"\n── TIER CHANGES (up or down) ({len(changed)}) ──")
    for _, r in changed[~changed["tier_escalated"]].iterrows():
        print(f"  {r['combination']:12s}: {r['tier_adult']} → {r['tier_peds']} "
              f"(Δ{r['delta_dQTc_ms']:+.1f}ms)")

    # ── Distribution ─────────────────────────────────────────────────────────
    print("\n── DISTRIBUTION ──")
    for label, col in [("Adult ref", "tier_adult"), ("Pediatric", "tier_peds")]:
        dist = df[col].value_counts()
        print(f"  {label}: " + " | ".join(
            f"{t}:{dist.get(t,0)}" for t in ["HIGH","MODERATE","LOW-MOD","LOW","PROTECTIVE"]))

    # ── Figure ────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(13, 6))

    # Scatter: adult vs pediatric dQTc
    ax = axes[0]
    colors = {"HIGH":"#8B1A1A","MODERATE":"#8B4A0A","LOW-MOD":"#0D5C6B","LOW":"#27682e","PROTECTIVE":"#0D5C6B"}
    for _, r in df.iterrows():
        c = colors.get(r["tier_peds"],"#7A7468")
        ax.scatter(r["dQTc_adult_ms"], r["dQTc_peds_ms"], color=c,
                   s=60, alpha=0.85, zorder=3)
        if r["tier_changed"]:
            ax.annotate(r["combination"], (r["dQTc_adult_ms"], r["dQTc_peds_ms"]),
                        fontsize=7, xytext=(4,3), textcoords="offset points", color="#1B3A6B")
    ax.plot([-5,35],[-5,35], color="#CDD3DE", lw=1, linestyle="--", zorder=1)
    ax.axhline(10, color="#8B4A0A", lw=0.8, linestyle=":", alpha=0.6)
    ax.axvline(10, color="#8B4A0A", lw=0.8, linestyle=":", alpha=0.6)
    ax.set_xlabel("Adult-reference ΔQTc (ms)", fontsize=11)
    ax.set_ylabel("Pediatric-Cmax ΔQTc (ms)", fontsize=11)
    ax.set_title("A.  Adult vs Pediatric ΔQTc Predictions", fontsize=11,
                 fontweight="bold", color="#1B3A6B", loc="left")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    # Legend
    for tier_lbl, col in colors.items():
        ax.scatter([],[], color=col, s=40, label=tier_lbl)
    ax.legend(fontsize=8, framealpha=0.9)

    # Bar chart: delta dQTc for top 15
    ax2 = axes[1]
    top15 = df.head(15)
    bar_colors = [colors.get(t,"#7A7468") for t in top15["tier_peds"]]
    bars = ax2.barh(range(len(top15)), top15["delta_dQTc_ms"],
                    color=bar_colors, height=0.6, edgecolor="white")
    ax2.set_yticks(range(len(top15)))
    ax2.set_yticklabels(top15["combination"], fontsize=9)
    ax2.invert_yaxis()
    ax2.axvline(0, color="#CDD3DE", lw=1)
    ax2.set_xlabel("ΔQTc change: pediatric minus adult (ms)", fontsize=10)
    ax2.set_title("B.  ΔQTc increase with pediatric Cmax (top 15)", fontsize=11,
                  fontweight="bold", color="#1B3A6B", loc="left")
    ax2.spines["top"].set_visible(False); ax2.spines["right"].set_visible(False)
    for i, (_, r) in enumerate(top15.iterrows()):
        if r["tier_changed"]:
            ax2.text(r["delta_dQTc_ms"]+0.3, i,
                     f"→{r['tier_peds']}", fontsize=7, va="center", color="#1B3A6B")

    plt.suptitle("Supplementary Figure S3. Pediatric PK Sensitivity Analysis\n"
                 "Comparison of adult-reference vs pediatric-specific free Cmax predictions.\n"
                 "Dashed line = identity (no change). Dotted lines = MODERATE threshold (10 ms).",
                 fontsize=9, y=1.01, color="#1B3A6B", fontweight="bold")
    plt.tight_layout()
    plt.savefig(str(OUT_FIG), dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved -> {OUT_FIG}")

if __name__ == "__main__":
    main()
