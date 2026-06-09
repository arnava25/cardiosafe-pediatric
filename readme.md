# CardioSafe Pediatric

**Computational cardiac risk stratification for psychiatric polypharmacy in children and adolescents**

[![status](https://img.shields.io/badge/status-manuscript%20in%20preparation-blue)](https://github.com/arnava25/cardiosafe-pediatric)
[![model](https://img.shields.io/badge/model-O'Hara--Rudy%202011-navy)](https://doi.org/10.1371/journal.pcbi.1002061)
[![data](https://img.shields.io/badge/validation-FDA%20FAERS%202015--2024-green)](https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html)

---

## Overview

Children and adolescents with psychiatric conditions are frequently prescribed multiple medications simultaneously: stimulants, antipsychotics, antidepressants, and alpha-2 agonists, often in combination. These drugs carry real cardiac risks, but existing safety frameworks focus almost exclusively on hERG channel block, missing the sympathomimetic and autonomic mechanisms that dominate risk in the most common pediatric polypharmacy patterns.

CardioSafe Pediatric fills that gap. It couples **O'Hara-Rudy (ORd) 2011 ventricular action potential modeling** with a three-pathway drug effect architecture and validates predictions against **FDA FAERS pharmacovigilance data** across 552,832 pediatric cases.

### Central finding

> **28 of 31 HIGH or MODERATE risk combinations have IKr block < 5%.** Sympathomimetic heart rate elevation — not hERG block — is the dominant cardiac risk mechanism for stimulant-containing polypharmacy. This mechanism is entirely invisible to standard ion channel screening.

---

## Key Results

### Polypharmacy risk grid (84 combinations, 500-beat steady state)

| Risk tier | Combinations | ΔQTc threshold |
|---|---|---|
| **HIGH** | 2 | ≥ 20 ms |
| **MODERATE** | 29 | 10–19 ms |
| **LOW-MOD** | 14 | 5–9 ms |
| **LOW / protective** | 39 | < 5 ms |

Selected findings:

| Combination | ΔQTc | IKr block | Mechanism |
|---|---|---|---|
| MPH + AMP | +19.1 ms | 0% | Sympathomimetic (pure) |
| MPH + ARI | +18.1 ms | 4.22% | Sympathomimetic + hERG |
| MPH + NOR | +14.8 ms | 2.65% | Sympathomimetic + hERG |
| MPH + RIS | +10.0 ms | 0.36% | Sympathomimetic |
| MPH + FLU | +10.0 ms | 0.33% | Sympathomimetic |
| RIS + SER | +1.5 ms | 0.55% | Minimal |
| ESC + CLO | −4.2 ms | 0.03% | Autonomic (bradycardic) |
| CLO + GUA | −1.1 ms | 0% | Autonomic (bradycardic) |

### FAERS pharmacovigilance validation

FDA FAERS 2015Q1–2024Q4 · 16.1M cases · 552,832 pediatric cases (age < 18)

**11/12 drugs** showed cardiac pharmacovigilance signal in pediatric population. Per-drug ROR highlights:

| Drug | n | ROR | 95% CI |
|---|---|---|---|
| Nortriptyline | 200 | **10.44** | 6.47–16.86 |
| Guanfacine | 2,936 | 4.63 | 3.86–5.55 |
| Quetiapine | 4,936 | 4.21 | 3.63–4.88 |
| Fluoxetine | 5,894 | 4.08 | 3.55–4.68 |
| Methylphenidate | 10,523 | 1.70 | 1.46–1.99 |
| Clonidine | 4,306 | 1.22 | 0.92–1.61 *(no signal — consistent with model)* |

Key combination signals:

| Combination | n | FAERS ROR | Model ΔQTc | Concordant |
|---|---|---|---|---|
| MPH + SER | 578 | **12.79** [9.84–16.62] | +9.6 ms | ✓ |
| MPH + ARI | — | 8.15 [6.10–10.90] | +18.1 ms | ✓ |
| MPH + QUE | — | 5.25 [3.04–9.08] | +11.0 ms | ✓ |
| ARI + GUA | 359 | 10.25 [7.14–14.71] | +2.7 ms | ✗ (conduction gap) |
| QUE + FLU | — | 9.18 [6.75–12.49] | +3.0 ms | ✗ (CYP2D6 gap) |

**Mechanism-stratified concordance:** sensitivity 0.80 for stimulant-containing combinations (Fisher p=0.04).

### ECG calibration

Steady-state APD90 at 60 bpm (500 beats): **331.7 ms**. Systematic offset vs. adolescent reference (Rijnbeek et al. 2014, ages 12–16): **67 ± 9 ms**. Offset cancels exactly in ΔQTc calculations — confirmed numerically to floating-point precision. Fridericia correction produces identical ΔQTc at fixed CL=1000 ms.

---

## Drug Architecture

### Drug list

| Code | Drug | Class | IC50 quality | Mechanism |
|---|---|---|---|---|
| MPH | Methylphenidate | Stimulant | C* | Sympathomimetic |
| AMP | Amphetamine | Stimulant | C* | Sympathomimetic |
| RIS | Risperidone | Antipsychotic | A | hERG block |
| QUE | Quetiapine | Antipsychotic | A | hERG block |
| ARI | Aripiprazole | Antipsychotic | A | hERG block |
| SER | Sertraline | SSRI | A | Mild hERG block |
| FLU | Fluoxetine | SSRI | A | Mild hERG block + CYP2D6 inhibition |
| ESC | Escitalopram | SSRI | B | Mild hERG block |
| CLO | Clonidine | Alpha-2 agonist | B | Autonomic modulation |
| GUA | Guanfacine | Alpha-2 agonist | C* | Autonomic modulation |
| IMI | Imipramine | TCA | A | hERG block |
| NOR | Nortriptyline | TCA | A | hERG block |

Quality A = patch clamp confirmed · B = mixed/binding assay · C = estimate

### Three-pathway drug effect model

**Pathway 1 — hERG/IKr block**
Competitive binding: `block_i = C_free / (C_free + IC50_i)`. Combined: `1 - Π(1 - block_i)`. Applied as GKr scaling.

**Pathway 2 — Sympathomimetic** (stimulants)
10% CL reduction + 15% GCaL upregulation. Consistent with published pediatric methylphenidate HR effects (3–8 bpm). ΔQTc predictions reflect Bazett-apparent change, not true repolarization prolongation.

**Pathway 3 — Autonomic modulation** (alpha-2 agonists)
10% CL increase (bradycardia) + 5% GNa reduction. Produces negative ΔQTc via Bazett correction.

---

## Methods

| Component | Details |
|---|---|
| AP model | O'Hara-Rudy 2011, 41-variable ODE system |
| Integration | `scipy.integrate.odeint`, mxstep=5000 |
| Simulation | 500 beats, CL=1000 ms (60 bpm) |
| APD90 | Duration from upstroke peak to 90% repolarization |
| QTc | Bazett correction: APD90 / √(RR_s) |
| IC50 sources | ChEMBL v37 + Witchel 2002, Redfern 2003, Kongsamut 2002, Polak 2009, Kramer 2013, Perrin 2008 |
| Cmax values | Adult PK data + protein binding fractions |
| FAERS | 40 quarters 2015Q1–2024Q4, pediatric filter age < 18 |
| ROR | 0.5 continuity correction, signal = CI lower bound > 1.0 |
| Concordance | Fisher exact + 10,000-iteration permutation test |

---

## Known Limitations

- All free Cmax values are **adult-derived**. CYP2D6 developmental variation in adolescents may shift concentrations of risperidone, aripiprazole, fluoxetine, nortriptyline, and imipramine by 1.5–3×
- ORd model is **adult ventricular** — no validated adolescent-specific AP model exists in the literature
- **No PR interval / AV conduction pathway** — accounts for systematic guanfacine discordances in FAERS
- **No CYP2D6 PK interaction terms** — accounts for fluoxetine/quetiapine and fluoxetine/aripiprazole discordances
- Sympathomimetic parameters (10% CL, 15% GCaL) are phenomenological — full PKA-mediated beta-adrenergic cascade not implemented
- Active metabolites (norfluoxetine, 9-OH-risperidone) not modeled
- Sex differences in IKs (hormonal modulation) not modeled

---

## Repository Structure

```
cardiosafe-pediatric/
├── src/
│   ├── ord_model.py              # ORd AP model + drug block + simulation runner
│   ├── risk_grid.py              # Full 84-combination polypharmacy sweep
│   ├── faers.py                  # FAERS download + parse + ROR pipeline
│   ├── ecg_calibration.py        # APD90 vs pediatric ECG reference
│   ├── concordance_stats.py      # Model-FAERS classification metrics
│   ├── concordance_stats_stratified.py  # Mechanism-stratified concordance
│   ├── pull.py / pull2.py        # ChEMBL hERG IC50 pull
│   ├── lit.py                    # Literature IC50 curation
│   └── master.py                 # Consolidate ChEMBL + literature params
├── data/
│   ├── herg_master_params.csv    # Master drug parameterization table
│   ├── herg_literature.csv       # Literature IC50 values
│   └── faers_cache/              # FAERS parquet cache (gitignored — ~2GB)
├── results/
│   ├── risk_grid_results.csv     # All 84 combinations ranked by ΔQTc
│   ├── pairwise_results.csv      # Pairs only
│   ├── triple_results.csv        # Clinical triples
│   ├── risk_matrix.csv           # 12×12 pairwise ΔQTc matrix
│   ├── ecg_calibration.csv       # Calibration data
│   └── faers/
│       ├── faers_drug_ror.csv    # Per-drug ROR
│       ├── faers_combo_ror.csv   # Pairwise combination ROR
│       ├── faers_model_alignment.csv  # Model vs FAERS comparison
│       └── concordance_statistics.csv
└── docs/
    ├── report.html               # Interactive clinical risk dashboard (v1.1)
    ├── faers_analysis_memo.md    # FAERS results interpretation
    └── ecg_calibration_methods.md  # Manuscript methods paragraph
```

---

## Installation

```bash
pip install numpy scipy pandas matplotlib tqdm pyarrow requests chembl-webresource-client
```

## Usage

```python
from src.ord_model import run_simulation

# Baseline
baseline = run_simulation(None, n_beats=500)

# Single combination
result = run_simulation(
    {"Methylphenidate": "therapeutic", "Aripiprazole": "therapeutic"},
    n_beats=500
)
print(f"ΔQTc: {result['QTc'] - baseline['QTc']:+.1f} ms")
```

```bash
# Full polypharmacy sweep
python3 src/risk_grid.py

# FAERS pipeline
python3 src/faers.py --all

# ECG calibration
python3 src/ecg_calibration.py

# Concordance statistics
python3 src/concordance_stats.py
python3 src/concordance_stats_stratified.py
```

---

## Reference

O'Hara T, Virag L, Varro A, Rudy Y. Simulation of the undiseased human cardiac ventricular action potential: model formulation and experimental validation. *PLoS Comput Biol.* 2011;7(5):e1002061.

---

## Status

Manuscript in preparation. Target journals: *npj Digital Medicine* · *Clinical Pharmacology & Therapeutics* · *JACAP*

## Author

Arnav Amit · Independent researcher · arnav.amit1@gmail.com
