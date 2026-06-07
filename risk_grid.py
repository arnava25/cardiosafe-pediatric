"""
Full Polypharmacy Risk Grid — CardioSafe Pediatric
Runs all pairwise and clinically relevant triple combinations
from the psychiatric drug list and produces a risk stratification table.

Run after ord_model.py is in the same directory.
Usage: python3 risk_grid.py
Output: risk_grid_results.csv, risk_stratification_table.csv
"""

import numpy as np
import pandas as pd
import itertools
import time
from ord_model import run_simulation, run_polypharmacy_sweep, load_drug_params

PARAMS_PATH = "herg_master_params.csv"
N_BEATS     = 50    # increase to 200 for manuscript figures
CL          = 1000.0

# ── DRUG LIST ─────────────────────────────────────────────────────────────────
# Grouped by class for interpretability
DRUGS = {
    # Stimulants
    "MPH":   "Methylphenidate",
    "AMP":   "Amphetamine",
    # Antipsychotics
    "RIS":   "Risperidone",
    "QUE":   "Quetiapine",
    "ARI":   "Aripiprazole",
    # SSRIs
    "SER":   "Sertraline",
    "FLU":   "Fluoxetine",
    "ESC":   "Escitalopram",
    # Alpha-2
    "CLO":   "Clonidine",
    "GUA":   "Guanfacine",
    # TCAs
    "IMI":   "Imipramine",
    "NOR":   "Nortriptyline",
}

DRUG_NAMES = list(DRUGS.values())
ABBREVS    = {v: k for k, v in DRUGS.items()}

# ── CLINICALLY RELEVANT TRIPLES ───────────────────────────────────────────────
# Based on common adolescent psychiatric polypharmacy patterns
CLINICAL_TRIPLES = [
    # ADHD + anxiety + mood
    ["Methylphenidate", "Aripiprazole",  "Sertraline"],
    ["Methylphenidate", "Aripiprazole",  "Fluoxetine"],
    ["Methylphenidate", "Aripiprazole",  "Escitalopram"],
    ["Methylphenidate", "Risperidone",   "Sertraline"],
    ["Methylphenidate", "Risperidone",   "Fluoxetine"],
    ["Methylphenidate", "Quetiapine",    "Sertraline"],
    # ADHD + alpha-2 augmentation
    ["Methylphenidate", "Clonidine",     "Risperidone"],
    ["Methylphenidate", "Clonidine",     "Aripiprazole"],
    ["Methylphenidate", "Guanfacine",    "Aripiprazole"],
    ["Methylphenidate", "Guanfacine",    "Sertraline"],
    # Amphetamine combos
    ["Amphetamine",     "Aripiprazole",  "Sertraline"],
    ["Amphetamine",     "Risperidone",   "Sertraline"],
    # TCA-containing
    ["Imipramine",      "Methylphenidate","Risperidone"],
    ["Nortriptyline",   "Methylphenidate","Aripiprazole"],
    ["Imipramine",      "Methylphenidate","Sertraline"],
    # Antipsychotic + SSRI + alpha-2
    ["Risperidone",     "Sertraline",    "Clonidine"],
    ["Aripiprazole",    "Fluoxetine",    "Clonidine"],
    ["Quetiapine",      "Sertraline",    "Clonidine"],
]

def build_combo_dict(drug_list):
    return {d: "therapeutic" for d in drug_list}

def risk_flag(dQTc):
    if dQTc >= 20:  return "HIGH"
    if dQTc >= 10:  return "MODERATE"
    if dQTc >= 5:   return "LOW-MOD"
    return "LOW"

def risk_emoji(flag):
    return {"HIGH":"⚠️ HIGH", "MODERATE":"△ MOD", "LOW-MOD":"~ LOW-MOD", "LOW":"✓ LOW"}[flag]

# ── RUN BASELINE ──────────────────────────────────────────────────────────────
print("="*70)
print("CardioSafe Pediatric — Full Polypharmacy Risk Grid")
print("="*70)
print(f"Drugs: {len(DRUG_NAMES)} | Beats: {N_BEATS} | CL: {CL} ms")

t_start = time.time()
baseline = run_simulation(None, PARAMS_PATH, n_beats=N_BEATS, verbose=True)
bAPD = baseline["APD90"]
bQTc = baseline["QTc"]
print(f"\nBaseline: APD90={bAPD:.1f} ms | QTc={bQTc:.1f} ms")

# ── ALL PAIRWISE COMBINATIONS ─────────────────────────────────────────────────
print(f"\n── PAIRWISE COMBINATIONS ({len(list(itertools.combinations(DRUG_NAMES,2)))} total) ──")

pairs = list(itertools.combinations(DRUG_NAMES, 2))
pair_results = []

for drug_a, drug_b in pairs:
    combo = {drug_a: "therapeutic", drug_b: "therapeutic"}
    res = run_simulation(combo, PARAMS_PATH, n_beats=N_BEATS, verbose=False)
    dQTc = res["QTc"] - bQTc
    flag = risk_flag(dQTc)
    pair_results.append({
        "drug_A":         drug_a,
        "drug_B":         drug_b,
        "combination":    f"{ABBREVS[drug_a]}+{ABBREVS[drug_b]}",
        "n_drugs":        2,
        "APD90_ms":       round(res["APD90"], 1),
        "QTc_ms":         round(res["QTc"], 1),
        "ΔAPD90_ms":      round(res["APD90"] - bAPD, 1),
        "ΔQTc_ms":        round(dQTc, 1),
        "IKr_block_pct":  round(res["IKr_block_pct"], 2),
        "risk_flag":      flag,
        "CL_eff_ms":      round(res["CL_effective"], 1),
    })
    print(f"  {risk_emoji(flag):12s} {drug_a:20s} + {drug_b:20s} | ΔQTc={dQTc:+.1f} ms")

pair_df = pd.DataFrame(pair_results).sort_values("ΔQTc_ms", ascending=False)

# ── CLINICAL TRIPLES ──────────────────────────────────────────────────────────
print(f"\n── CLINICAL TRIPLE COMBINATIONS ({len(CLINICAL_TRIPLES)} total) ──")

triple_results = []

for triple in CLINICAL_TRIPLES:
    combo = build_combo_dict(triple)
    res = run_simulation(combo, PARAMS_PATH, n_beats=N_BEATS, verbose=False)
    dQTc = res["QTc"] - bQTc
    flag = risk_flag(dQTc)
    abbr = "+".join(ABBREVS[d] for d in triple)
    triple_results.append({
        "combination":    abbr,
        "drugs":          " + ".join(triple),
        "n_drugs":        3,
        "APD90_ms":       round(res["APD90"], 1),
        "QTc_ms":         round(res["QTc"], 1),
        "ΔAPD90_ms":      round(res["APD90"] - bAPD, 1),
        "ΔQTc_ms":        round(dQTc, 1),
        "IKr_block_pct":  round(res["IKr_block_pct"], 2),
        "risk_flag":      flag,
        "CL_eff_ms":      round(res["CL_effective"], 1),
    })
    print(f"  {risk_emoji(flag):12s} {abbr:30s} | ΔQTc={dQTc:+.1f} ms")

triple_df = pd.DataFrame(triple_results).sort_values("ΔQTc_ms", ascending=False)

# ── COMBINED RESULTS ──────────────────────────────────────────────────────────
all_results = pd.concat([pair_df, triple_df], ignore_index=True)
all_results = all_results.sort_values("ΔQTc_ms", ascending=False)

total_time = time.time() - t_start

# ── RISK STRATIFICATION SUMMARY ───────────────────────────────────────────────
print("\n" + "="*70)
print("RISK STRATIFICATION SUMMARY")
print("="*70)
print(f"Total combinations tested: {len(all_results)}")
print(f"Total runtime: {total_time/60:.1f} min\n")

for flag in ["HIGH", "MODERATE", "LOW-MOD", "LOW"]:
    subset = all_results[all_results["risk_flag"]==flag]
    print(f"{risk_emoji(flag)} ({len(subset)} combinations):")
    for _, row in subset.iterrows():
        combo = row.get("combination","")
        print(f"    {combo:35s} ΔQTc={row['ΔQTc_ms']:+.1f} ms  "
              f"IKr_block={row['IKr_block_pct']:.2f}%")
    print()

# ── PAIRWISE HEATMAP DATA (for visualization) ─────────────────────────────────
print("── PAIRWISE ΔQTc MATRIX ──")
matrix_data = {}
for _, row in pair_df.iterrows():
    a, b = row["drug_A"], row["drug_B"]
    abbr_a, abbr_b = ABBREVS[a], ABBREVS[b]
    dq = row["ΔQTc_ms"]
    if abbr_a not in matrix_data:
        matrix_data[abbr_a] = {}
    if abbr_b not in matrix_data:
        matrix_data[abbr_b] = {}
    matrix_data[abbr_a][abbr_b] = dq
    matrix_data[abbr_b][abbr_a] = dq

abbr_list = list(ABBREVS.values())
matrix_df = pd.DataFrame(matrix_data, index=abbr_list, columns=abbr_list).fillna(0)

print(matrix_df.round(1).to_string())

# ── KEY CLINICAL INSIGHTS ──────────────────────────────────────────────────────
print("\n" + "="*70)
print("KEY MECHANISTIC INSIGHTS")
print("="*70)

high_herg   = all_results[all_results["IKr_block_pct"] > 5]
high_dqtc   = all_results[all_results["ΔQTc_ms"] >= 10]
mech_driven = high_dqtc[high_dqtc["IKr_block_pct"] < 5]

print(f"\n1. Combinations with ΔQTc≥10ms but IKr block <5% (sympathomimetic-driven):")
for _, row in mech_driven.iterrows():
    print(f"   {row['combination']:35s} ΔQTc={row['ΔQTc_ms']:+.1f}ms  IKr={row['IKr_block_pct']:.2f}%")

print(f"\n2. Combinations with highest IKr block:")
for _, row in all_results.nlargest(5,"IKr_block_pct").iterrows():
    print(f"   {row['combination']:35s} IKr={row['IKr_block_pct']:.2f}%  ΔQTc={row['ΔQTc_ms']:+.1f}ms")

print(f"\n3. Stimulant-free combinations with ΔQTc≥5ms:")
no_stim = all_results[
    ~all_results["combination"].str.contains("MPH|AMP")
    & (all_results["ΔQTc_ms"] >= 5)
]
for _, row in no_stim.iterrows():
    print(f"   {row['combination']:35s} ΔQTc={row['ΔQTc_ms']:+.1f}ms")

print(f"\n4. Adolescent-specific flags (all Cmax values adult-derived):")
print(f"   - All ΔQTc estimates carry ~20-30% uncertainty from Cmax extrapolation")
print(f"   - CYP2D6 developmental variation could shift Cmax 1.5-3x for:")
print(f"     risperidone, aripiprazole, fluoxetine, nortriptyline, imipramine")
print(f"   - Adolescent baseline QTc differs from adult (shorter by ~10-15ms)")
print(f"   - Pubertal hormonal effects on IKs not modeled (flag for future work)")

# ── SAVE ──────────────────────────────────────────────────────────────────────
all_results.to_csv("risk_grid_results.csv", index=False)
matrix_df.to_csv("risk_matrix.csv")
pair_df.to_csv("pairwise_results.csv", index=False)
triple_df.to_csv("triple_results.csv", index=False)

print(f"\nSaved:")
print(f"  risk_grid_results.csv  — all combinations ranked by ΔQTc")
print(f"  risk_matrix.csv        — pairwise ΔQTc heatmap data")
print(f"  pairwise_results.csv   — pairs only")
print(f"  triple_results.csv     — triples only")
print(f"\nTotal runtime: {total_time/60:.1f} min")
print(f"Next: run visualize_risk.py to generate manuscript figures")
