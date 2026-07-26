"""
pk_pediatric_genuine.py — CardioSafe Pediatric
================================================
Pediatric PK sensitivity analysis in the GENUINE delta-APD90 domain.

This is a corrected fork of pk_pediatric_sensitivity.py. The original computed
res["QTc"] (Bazett-corrected), which is the metric the manuscript argues against.
This version reads res["APD90"] so the escalations are expressed in the same
genuine repolarization currency the paper uses everywhere else. Bazett dQTc is
retained as a clearly-labeled reference column only.

Tiers match the manuscript genuine-domain definition:
  elevated   : >= 10 ms
  borderline : 5 to 9 ms
  minimal    : 0 to 5 ms
  shortening : < 0 ms

Run from repo root:
    python3 src/pediatric_genuine.py
Default 200 beats to scout; set N_BEATS=500 for final numbers that line up with
the production grid.

NOTE (July 2026): PARAMS_PATH below points at data/herg_master_params.csv, which
was withdrawn in the June 2026 audit and moved to
archive/params_invalidated_202606/. This script will not run until it is
repointed at the rebuilt table in params/.
"""

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent   # file lives in src/
sys.path.insert(0, str(ROOT / "src"))

import pandas as pd
import numpy as np
from ord_model import run_simulation

PARAMS_PATH = str(ROOT / "data" / "herg_master_params.csv")
OUT_CSV     = ROOT / "results" / "pk_pediatric_genuine.csv"
N_BEATS     = 200   # set to 500 for final, to match the production grid

# ── Pediatric free Cmax values (nM) — verbatim from original analysis ─────────
PEDIATRIC_CMAX = {
    "Methylphenidate": 21.9, "Amphetamine": 148.0, "Risperidone": 3.87,
    "Quetiapine": 86.6, "Aripiprazole": 11.9, "Sertraline": 10.8,
    "Fluoxetine": 30.4, "Escitalopram": 13.2, "Clonidine": 0.5,
    "Guanfacine": 12.2, "Imipramine": 53.5, "Nortriptyline": 30.4,
}
CONFIDENCE = {
    "Methylphenidate": "MEDIUM", "Amphetamine": "MEDIUM", "Risperidone": "HIGH",
    "Quetiapine": "MEDIUM", "Aripiprazole": "LOW", "Sertraline": "HIGH",
    "Fluoxetine": "MEDIUM", "Escitalopram": "MEDIUM", "Clonidine": "LOW",
    "Guanfacine": "HIGH", "Imipramine": "MEDIUM", "Nortriptyline": "LOW",
}
DRUG_CODES = {
    "Methylphenidate":"MPH","Amphetamine":"AMP","Risperidone":"RIS",
    "Quetiapine":"QUE","Aripiprazole":"ARI","Sertraline":"SER",
    "Fluoxetine":"FLU","Escitalopram":"ESC","Clonidine":"CLO",
    "Guanfacine":"GUA","Imipramine":"IMI","Nortriptyline":"NOR",
}

# ── Genuine-domain tiers (manuscript definition) ─────────────────────────────
RANK = {"shortening": 0, "minimal": 1, "borderline": 2, "elevated": 3}
def tier(d):
    if d >= 10: return "elevated"
    if d >= 5:  return "borderline"
    if d >= 0:  return "minimal"
    return "shortening"

FOCUS_COMBOS = [
    ("Methylphenidate","Aripiprazole"),("Amphetamine","Aripiprazole"),
    ("Methylphenidate","Amphetamine"),("Methylphenidate","Quetiapine"),
    ("Methylphenidate","Sertraline"),("Amphetamine","Sertraline"),
    ("Methylphenidate","Nortriptyline"),("Amphetamine","Nortriptyline"),
    ("Aripiprazole","Nortriptyline"),("Aripiprazole","Imipramine"),
    ("Aripiprazole","Sertraline"),("Aripiprazole","Quetiapine"),
    ("Aripiprazole","Fluoxetine"),("Quetiapine","Nortriptyline"),
    ("Methylphenidate","Risperidone"),("Amphetamine","Risperidone"),
    ("Methylphenidate","Imipramine"),("Amphetamine","Imipramine"),
    ("Fluoxetine","Nortriptyline"),("Risperidone","Nortriptyline"),
    ("Quetiapine","Aripiprazole"),("Quetiapine","Fluoxetine"),
    ("Aripiprazole","Escitalopram"),("Methylphenidate","Fluoxetine"),
    ("Quetiapine","Sertraline"),("Methylphenidate","Escitalopram"),
    ("Risperidone","Aripiprazole"),("Imipramine","Nortriptyline"),
    ("Aripiprazole","Guanfacine"),("Quetiapine","Guanfacine"),
]

def main():
    print("CardioSafe Pediatric — Pediatric PK Sensitivity (GENUINE delta-APD90)")
    print(f"Beats: {N_BEATS} | Combinations: {len(FOCUS_COMBOS)}\n")

    base = run_simulation(None, PARAMS_PATH, n_beats=N_BEATS, verbose=False)
    base_apd = base["APD90"]
    base_qtc = base["QTc"]
    print(f"Baseline APD90: {base_apd:.1f} ms | baseline QTc: {base_qtc:.1f} ms\n")

    rows = []
    for drug_a, drug_b in FOCUS_COMBOS:
        a, b = DRUG_CODES[drug_a], DRUG_CODES[drug_b]

        # adult reference (therapeutic Cmax resolved from params CSV)
        r_ad = run_simulation({drug_a:"therapeutic", drug_b:"therapeutic"},
                              PARAMS_PATH, n_beats=N_BEATS, verbose=False)
        dapd_ad = round(r_ad["APD90"] - base_apd, 1)      # GENUINE
        dqtc_ad = round(r_ad["QTc"]   - base_qtc, 1)      # Bazett (reference only)

        # pediatric Cmax
        r_pd = run_simulation({drug_a:PEDIATRIC_CMAX[drug_a], drug_b:PEDIATRIC_CMAX[drug_b]},
                              PARAMS_PATH, n_beats=N_BEATS, verbose=False)
        dapd_pd = round(r_pd["APD90"] - base_apd, 1)      # GENUINE
        dqtc_pd = round(r_pd["QTc"]   - base_qtc, 1)      # Bazett (reference only)

        t_ad, t_pd = tier(dapd_ad), tier(dapd_pd)
        ca, cb = CONFIDENCE[drug_a], CONFIDENCE[drug_b]
        conf = "HIGH" if ca=="HIGH" and cb=="HIGH" else ("LOW" if "LOW" in (ca,cb) else "MEDIUM")

        rows.append({
            "combination": f"{a}+{b}",
            "genuine_dAPD_adult": dapd_ad, "tier_adult": t_ad,
            "genuine_dAPD_peds": dapd_pd,  "tier_peds": t_pd,
            "delta_genuine": round(dapd_pd - dapd_ad, 1),
            "escalated": RANK[t_pd] > RANK[t_ad],
            "bazett_dQTc_adult_REF": dqtc_ad, "bazett_dQTc_peds_REF": dqtc_pd,
            "conf": conf,
        })
        print(f"  {a}+{b:<14} genuine adult {dapd_ad:+5.1f} ({t_ad:10s}) "
              f"-> peds {dapd_pd:+5.1f} ({t_pd:10s})  [conf {conf}]")

    df = pd.DataFrame(rows).sort_values("genuine_dAPD_peds", ascending=False)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nSaved -> {OUT_CSV}")

    esc = df[df["escalated"]]
    print(f"\n── GENUINE-DOMAIN TIER ESCALATIONS ({len(esc)}) ──")
    if len(esc) == 0:
        print("  none — pediatric exposure does not move any combo to a higher genuine tier")
    for _, r in esc.iterrows():
        print(f"  {r['combination']:14s}: {r['tier_adult']} -> {r['tier_peds']}  "
              f"(genuine {r['genuine_dAPD_adult']:+.1f} -> {r['genuine_dAPD_peds']:+.1f} ms, conf {r['conf']})")

    print("\n── DISTRIBUTION (genuine tiers) ──")
    for label, col in [("Adult", "tier_adult"), ("Peds", "tier_peds")]:
        d = df[col].value_counts()
        print(f"  {label}: " + " | ".join(f"{t}:{d.get(t,0)}"
              for t in ["elevated","borderline","minimal","shortening"]))

    print("\nSanity check: the genuine_dAPD_adult column should approximately match")
    print("risk_grid_results.csv for these pairs (a few ms high at 200 vs 500 beats).")

if __name__ == "__main__":
    main()
