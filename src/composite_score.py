"""
composite_score.py — CardioSafe Pediatric
==========================================
Computes a composite cardiac risk score (0-100) for each drug combination
integrating three independent signal components:

  1. Delta-QTc component    (weight 0.50) — mechanistic AP model prediction
  2. IKr block component    (weight 0.20) — hERG channel occupancy
  3. FAERS ROR component    (weight 0.30) — real-world pharmacovigilance signal

Plus binary flag penalties:
  - CYP2D6 interaction flag (+15 pts): fluoxetine + CYP2D6 substrate
  - Conduction flag         (+10 pts): guanfacine in combination

Score is capped at 100. Risk labels:
  75-100: HIGH
  50-74:  MODERATE
  25-49:  LOW-MOD
  0-24:   LOW

Usage:
    python3 src/composite_score.py
    python3 src/composite_score.py --show-top 20

Output:
    results/composite_scores.csv
    results/composite_score_memo.md
"""

import sys
import argparse
import math
from pathlib import Path
from datetime import datetime

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import pandas as pd
import numpy as np

GRID_CSV   = ROOT / "results" / "risk_grid_results.csv"
FAERS_CSV  = ROOT / "results" / "faers" / "faers_combo_ror.csv"
OUT_CSV    = ROOT / "results" / "composite_scores.csv"
OUT_MEMO   = ROOT / "results" / "composite_score_memo.md"

# ── WEIGHTS ──────────────────────────────────────────────────────────────────
W_DQTC  = 0.50
W_IKR   = 0.20
W_FAERS = 0.30

# ── SCALING PARAMETERS ───────────────────────────────────────────────────────
# Delta-QTc: scale so that 20ms = 100, <0 = 0
DQTC_MAX = 20.0

# IKr block: scale so that 10% block = 100
IKR_MAX  = 10.0

# FAERS ROR: log scale, ROR=1 = 0, ROR=30 = 100
# log(30)/log(30) = 1.0; log(1) = 0
FAERS_ROR_REF = 30.0   # ROR at which FAERS component = 100

# ── FLAGS ────────────────────────────────────────────────────────────────────
CYP2D6_PENALTY    = 15
CONDUCTION_PENALTY = 10

# CYP2D6 substrate drugs (fluoxetine inhibits their metabolism)
CYP2D6_SUBSTRATES = {"QUE", "ARI", "RIS", "NOR", "IMI"}
CYP2D6_INHIBITOR  = "FLU"
CONDUCTION_DRUG   = "GUA"

DRUG_CODES = {
    "Methylphenidate":"MPH","Amphetamine":"AMP","Risperidone":"RIS",
    "Quetiapine":"QUE","Aripiprazole":"ARI","Sertraline":"SER",
    "Fluoxetine":"FLU","Escitalopram":"ESC","Clonidine":"CLO",
    "Guanfacine":"GUA","Imipramine":"IMI","Nortriptyline":"NOR",
}


def score_label(s):
    if s >= 75: return "HIGH"
    if s >= 50: return "MODERATE"
    if s >= 25: return "LOW-MOD"
    return "LOW"


def dqtc_component(dqtc):
    """Scale delta-QTc to 0-100. Negative values = 0."""
    if dqtc <= 0:
        return 0.0
    return min(100.0, (dqtc / DQTC_MAX) * 100.0)


def ikr_component(ikr_pct):
    """Scale IKr block % to 0-100."""
    return min(100.0, (ikr_pct / IKR_MAX) * 100.0)


def faers_component(ror, has_signal):
    """
    Log-scale FAERS ROR to 0-100.
    No signal (CI lower bound <= 1) = 0.
    ROR = FAERS_ROR_REF = 100.
    """
    if not has_signal or ror <= 1.0:
        return 0.0
    return min(100.0, (math.log(ror) / math.log(FAERS_ROR_REF)) * 100.0)


def cyp2d6_flag(drug_a, drug_b):
    """True if fluoxetine + CYP2D6 substrate."""
    drugs = {drug_a, drug_b}
    return CYP2D6_INHIBITOR in drugs and bool(drugs & CYP2D6_SUBSTRATES)


def conduction_flag(drug_a, drug_b):
    """True if guanfacine is in the combination."""
    return CONDUCTION_DRUG in (drug_a, drug_b)


def compute_score(dqtc, ikr_pct, ror, has_signal, drug_a, drug_b):
    """Compute composite score 0-100."""
    c_dqtc  = dqtc_component(dqtc)
    c_ikr   = ikr_component(ikr_pct)
    c_faers = faers_component(ror, has_signal)

    base = W_DQTC * c_dqtc + W_IKR * c_ikr + W_FAERS * c_faers

    penalty = 0
    if cyp2d6_flag(drug_a, drug_b):
        penalty += CYP2D6_PENALTY
    if conduction_flag(drug_a, drug_b):
        penalty += CONDUCTION_PENALTY

    score = min(100.0, base + penalty)
    return round(score, 1), round(c_dqtc, 1), round(c_ikr, 1), round(c_faers, 1), penalty


def canonical_key(a, b):
    return "+".join(sorted([a, b]))


def load_grid():
    df = pd.read_csv(GRID_CSV)
    grid = {}
    for _, row in df.iterrows():
        combo = str(row.get("combination", "")).strip()
        parts = combo.split("+")
        if len(parts) != 2:
            continue
        key = canonical_key(parts[0], parts[1])
        try:
            grid[key] = {
                "dQTc":     float(row.get("ΔQTc_ms", row.get("delta_qtc", 0))),
                "IKr_pct":  float(row.get("IKr_block_pct", 0)),
                "drug_a":   parts[0],
                "drug_b":   parts[1],
            }
        except (ValueError, TypeError):
            continue
    return grid


def load_faers():
    faers = {}
    if not FAERS_CSV.exists():
        print(f"WARNING: {FAERS_CSV} not found")
        return faers
    df = pd.read_csv(FAERS_CSV)
    for _, row in df.iterrows():
        combo = str(row.get("combo", row.get("combination", ""))).strip()
        parts = combo.split("+")
        if len(parts) != 2:
            continue
        a = DRUG_CODES.get(parts[0].strip(), parts[0].strip())
        b = DRUG_CODES.get(parts[1].strip(), parts[1].strip())
        key = canonical_key(a, b)
        try:
            ror    = float(row.get("ROR", row.get("ror", 1.0)))
            ci_lo  = float(row.get("CI_lo", row.get("ci_lo", 0.0)))
            signal = ci_lo > 1.0
            faers[key] = {"ROR": ror, "signal": signal}
        except (ValueError, TypeError):
            continue
    return faers


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--show-top", type=int, default=15,
                        help="Number of top combinations to print")
    args = parser.parse_args()

    print("CardioSafe Pediatric — Composite Risk Score")
    print(f"Weights: dQTc={W_DQTC}  IKr={W_IKR}  FAERS={W_FAERS}")
    print(f"Flags:   CYP2D6 +{CYP2D6_PENALTY}pts  Conduction +{CONDUCTION_PENALTY}pts")
    print()

    grid  = load_grid()
    faers = load_faers()
    print(f"Loaded {len(grid)} combinations from risk grid")
    print(f"Loaded {len(faers)} FAERS ROR values")
    print()

    results = []
    for key, g in grid.items():
        drug_a = g["drug_a"]
        drug_b = g["drug_b"]
        f = faers.get(key, {"ROR": 1.0, "signal": False})

        score, c_dqtc, c_ikr, c_faers, penalty = compute_score(
            g["dQTc"], g["IKr_pct"], f["ROR"], f["signal"], drug_a, drug_b
        )

        results.append({
            "combination":      f"{drug_a}+{drug_b}",
            "composite_score":  score,
            "score_label":      score_label(score),
            "dQTc_ms":          round(g["dQTc"], 1),
            "IKr_block_pct":    round(g["IKr_pct"], 2),
            "FAERS_ROR":        round(f["ROR"], 2),
            "FAERS_signal":     f["signal"],
            "component_dQTc":   c_dqtc,
            "component_IKr":    c_ikr,
            "component_FAERS":  c_faers,
            "flag_penalty":     penalty,
            "cyp2d6_flag":      cyp2d6_flag(drug_a, drug_b),
            "conduction_flag":  conduction_flag(drug_a, drug_b),
        })

    df = pd.DataFrame(results).sort_values("composite_score", ascending=False)
    df.to_csv(OUT_CSV, index=False)
    print(f"Saved -> {OUT_CSV}")
    print()

    # Print top N
    print(f"── TOP {args.show_top} COMBINATIONS BY COMPOSITE SCORE ──")
    print(f"{'Combination':25s}  {'Score':6s}  {'Label':8s}  {'dQTc':6s}  {'IKr%':5s}  "
          f"{'ROR':6s}  {'Flags':12s}  Components")
    print("-" * 105)
    for _, row in df.head(args.show_top).iterrows():
        flags = ""
        if row["cyp2d6_flag"]:   flags += "CYP2D6 "
        if row["conduction_flag"]: flags += "COND"
        print(f"{row['combination']:25s}  {row['composite_score']:6.1f}  "
              f"{row['score_label']:8s}  {row['dQTc_ms']:+6.1f}  "
              f"{row['IKr_block_pct']:5.2f}  {row['FAERS_ROR']:6.2f}  "
              f"{flags:12s}  "
              f"dQTc={row['component_dQTc']:.0f} IKr={row['component_IKr']:.0f} "
              f"FAERS={row['component_FAERS']:.0f} pen={row['flag_penalty']}")

    # Distribution
    print()
    print("── SCORE DISTRIBUTION ──")
    for label in ["HIGH", "MODERATE", "LOW-MOD", "LOW"]:
        n = len(df[df["score_label"] == label])
        print(f"  {label:8s}: {n:3d} combinations")

    # Key comparisons showing FAERS lifting scores
    print()
    print("── NOTABLE SCORES (FAERS effect) ──")
    notable = ["MPH+SER", "SER+MPH", "MPH+ARI", "ARI+MPH",
               "ARI+GUA", "GUA+ARI", "QUE+FLU", "FLU+QUE",
               "QUE+GUA", "GUA+QUE"]
    shown = set()
    for combo in notable:
        parts = combo.split("+")
        key = canonical_key(parts[0], parts[1])
        match = df[df["combination"].apply(lambda x: canonical_key(*x.split("+")) == key)]
        if not match.empty and key not in shown:
            r = match.iloc[0]
            shown.add(key)
            flags = ""
            if r["cyp2d6_flag"]:    flags += "CYP2D6 "
            if r["conduction_flag"]: flags += "COND"
            print(f"  {r['combination']:25s}  score={r['composite_score']:.1f} ({r['score_label']:8s})  "
                  f"dQTc={r['dQTc_ms']:+.1f}ms  ROR={r['FAERS_ROR']:.2f}  flags={flags.strip() or 'none'}")

    write_memo(df)
    print(f"\nSaved -> {OUT_MEMO}")


def write_memo(df):
    lines = []
    lines.append("# Composite Risk Score — CardioSafe Pediatric")
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    lines.append("")
    lines.append("## Formula")
    lines.append(f"score = {W_DQTC} * dQTc_component + {W_IKR} * IKr_component + {W_FAERS} * FAERS_component")
    lines.append(f"      + CYP2D6_flag({CYP2D6_PENALTY}pts) + conduction_flag({CONDUCTION_PENALTY}pts)")
    lines.append("      capped at 100")
    lines.append("")
    lines.append("## Component Scaling")
    lines.append(f"- dQTc: (dQTc / {DQTC_MAX}ms) * 100, floor 0")
    lines.append(f"- IKr: (IKr_pct / {IKR_MAX}%) * 100")
    lines.append(f"- FAERS: (log(ROR) / log({FAERS_ROR_REF})) * 100, 0 if no signal")
    lines.append("")
    lines.append("## Score Labels")
    lines.append("75-100: HIGH | 50-74: MODERATE | 25-49: LOW-MOD | 0-24: LOW")
    lines.append("")
    lines.append("## Distribution")
    for label in ["HIGH", "MODERATE", "LOW-MOD", "LOW"]:
        n = len(df[df["score_label"] == label])
        lines.append(f"- {label}: {n} combinations")
    lines.append("")
    lines.append("## Top 20 Combinations")
    lines.append("| Combination | Score | Label | dQTc | IKr% | ROR | Flags |")
    lines.append("|---|---|---|---|---|---|---|")
    for _, row in df.head(20).iterrows():
        flags = []
        if row["cyp2d6_flag"]:    flags.append("CYP2D6")
        if row["conduction_flag"]: flags.append("COND")
        lines.append(f"| {row['combination']} | {row['composite_score']:.1f} | {row['score_label']} | "
                     f"{row['dQTc_ms']:+.1f} | {row['IKr_block_pct']:.2f}% | {row['FAERS_ROR']:.2f} | "
                     f"{', '.join(flags) or 'none'} |")
    lines.append("")
    lines.append("## Manuscript Note")
    lines.append("The composite score integrates mechanistic and epidemiological evidence.")
    lines.append("Key property: MPH+SER scores higher than its delta-QTc alone would suggest,")
    lines.append("because the FAERS ROR 12.79 contributes 0.30 weight. This correctly reflects")
    lines.append("the convergent evidence from two independent methods.")

    OUT_MEMO.write_text("\n".join(lines))


if __name__ == "__main__":
    main()
