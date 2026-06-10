# CardioSafe Pediatric

**Computational cardiac risk stratification for psychiatric polypharmacy in children and adolescents**

[![status](https://img.shields.io/badge/status-preprint%20in%20preparation-blue)](https://github.com/arnava25/cardiosafe-pediatric)
[![model](https://img.shields.io/badge/model-O'Hara--Rudy%202011-navy)](https://doi.org/10.1371/journal.pcbi.1002061)
[![data](https://img.shields.io/badge/validation-FDA%20FAERS%202015--2024-green)](https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html)

---

## Overview

Children and adolescents with psychiatric conditions are frequently prescribed multiple medications simultaneously: stimulants, antipsychotics, antidepressants, and alpha-2 agonists, often in combination. These drugs carry real cardiac risks, but existing safety frameworks focus almost exclusively on hERG channel block, missing the sympathomimetic and autonomic mechanisms that dominate risk in the most common pediatric polypharmacy patterns.

CardioSafe Pediatric fills that gap. It couples **O'Hara-Rudy (ORd) 2011 ventricular action potential modeling** with a three-pathway drug effect architecture and validates predictions against **FDA FAERS pharmacovigilance data** across 552,832 pediatric cases spanning 10 years (2015–2024).

> **This is a hypothesis-generating computational framework, not a validated clinical guideline.** Findings identify specific drug combinations warranting prospective ECG monitoring study.

---

## Central Findings

**Finding 1 — Sympathomimetic dominance:** 28 of 31 MODERATE-risk combinations have IKr block below 5%. Sympathomimetic heart rate elevation, not hERG block, is the dominant cardiac risk mechanism for stimulant-containing polypharmacy. This mechanism is entirely invisible to standard ion channel screening frameworks including CiPA and CredibleMeds.

**Finding 2 — Age stratification:** Stimulant cardiac adverse event signals in FDA FAERS are 23 to 50 times higher in children aged 6–12 than in adolescents aged 13–17 (MPH+SER: ROR 43.71 vs 1.86; MPH+ARI: ROR 23.71 vs 0.47). CYP2D6-mediated combination signals show the opposite pattern, consistent with hepatic enzyme ontogeny.

**Finding 3 — CredibleMeds gap:** Both stimulants (MPH, AMP) have no QTDrugs designation in CredibleMeds for the general population. MPH+ARI, the highest composite-risk combination (score 72.2), generates zero cardiac safety alert from CredibleMeds.

---

## Key Results

### Polypharmacy risk grid (84 combinations, 500-beat steady state)

| Risk tier | Combinations | ΔQTc threshold |
|---|---|---|
| HIGH | 0 | ≥ 20 ms |
| **MODERATE** | **31** | 10–19 ms |
| LOW-MOD | 18 | 5–9 ms |
| LOW / protective | 35 | < 5 ms |

No pairwise combinations reach HIGH tier at 500-beat steady state. MPH+AMP at +19.1 ms is the highest pairwise delta-QTc.

Selected pairwise results:

| Combination | ΔQTc | IKr block | Mechanism |
|---|---|---|---|
| MPH + AMP | +19.1 ms | 0.00% | Sympathomimetic (pure) |
| MPH + ARI | +18.1 ms | 4.22% | Sympathomimetic + hERG |
| ARI + NOR | +15.0 ms | 6.76% | hERG block |
| MPH + NOR | +14.8 ms | 2.65% | Sympathomimetic + hERG |
| MPH + QUE | +11.0 ms | 0.93% | Sympathomimetic |
| MPH + SER | +9.6 ms | 0.19% | Sympathomimetic |
| RIS + SER | +1.5 ms | 0.55% | Minimal |
| CLO + GUA | −1.5 ms | 0.00% | Autonomic (bradycardic) |

Selected clinical triples:

| Combination | ΔQTc |
|---|---|
| MPH + ARI + FLU | +19.0 ms |
| MPH + ARI + SER | +18.6 ms |
| MPH + ARI + ESC | +18.1 ms |

### FAERS pharmacovigilance validation

FDA FAERS 2015Q1–2024Q4 · 16.1M cases · 552,832 pediatric cases (age < 18)

11 of 12 drugs showed cardiac pharmacovigilance signal. Clonidine (ROR 1.22, no signal) is consistent with its modeled cardioprotective effect.

Key combination signals:

| Combination | n | FAERS ROR | Model ΔQTc | Notes |
|---|---|---|---|---|
| MPH + SER | 578 | 12.79 [9.84–16.62] | +9.6 ms | Strongest pediatric combo signal |
| ARI + GUA | 359 | 10.25 [7.14–14.71] | +7.3 ms | Conduction gap |
| MPH + ARI | 677 | 8.15 [6.10–10.90] | +18.1 ms | Concordant |
| MPH + QUE | 274 | 5.25 [3.04–9.08] | +11.0 ms | Concordant |

**Composite score AUC-ROC: 0.771 overall, 0.812 excluding named mechanistic-gap pairs.** Binary sensitivity 0.235 reflects threshold placement; no-gap sensitivity 0.444, non-stimulant non-gap kappa 0.344.

### Age-stratified pharmacovigilance

| Combination | Children ROR (6–12y) | Adolescents ROR (13–17y) | Direction |
|---|---|---|---|
| MPH + SER | 43.71 | 1.86 | Higher in children (23x) |
| MPH + ARI | 23.71 | 0.47 | Higher in children (50x) |
| QUE + GUA | 32.15 | 3.91 | Higher in children |
| QUE + FLU | 3.75 | 10.88 | Higher in adolescents (CYP2D6) |

### Developmental PK sensitivity

| Combination | Base ΔQTc | 1.5x ARI Cmax | Tier change |
|---|---|---|---|
| MPH + ARI | +18.1 ms | +21.4 ms | MODERATE → HIGH |
| ARI + SER | +9.5 ms | +13.0 ms | LOW-MOD → MODERATE |
| ARI + FLU | +10.0 ms | +13.0 ms | LOW-MOD → MODERATE |
| MPH + RIS | +10.0 ms | +10.0 ms | LOW-MOD → MODERATE |

### CYP2D6 Ito static model (fluoxetine co-administration)

Hepatic inlet [I]total = 212.5 nM, Ki = 170 nM. Predicted AUC ratios: QUE 1.68x, ARI 1.29x, RIS 1.75x, NOR 2.00x, IMI 1.44x. FLU+NOR escalates from LOW-MOD to MODERATE at adjusted Cmax.

### Composite risk score (0–100)

Integrates delta-QTc (0.50), IKr block (0.20), FAERS ROR (0.30), plus CYP2D6 flag (+15 pts) and conduction flag (+10 pts).

| Combination | Score | Notes |
|---|---|---|
| MPH + ARI | 72.2 | Top scorer |
| ARI + GUA | 57.2 | Conduction flag lifts score |
| MPH + SER | 46.9 | FAERS component 75 despite LOW-MOD ΔQTc |
| MPH + AMP | 47.8 | Zero FAERS signal despite highest raw ΔQTc |
| QUE + FLU | 44.6 | CYP2D6 flag lifts score |

### Sympathomimetic parameter sensitivity

5×5 grid (CL reduction 5–15%, GCaL upregulation 10–20%): MPH+ARI remained MODERATE or above in all 25 parameter combinations (+14.0 to +20.8 ms).

### CredibleMeds comparison

| Drug | CredibleMeds | CardioSafe gap |
|---|---|---|
| MPH | Special Risk (congenital LQTS only) | Sympathomimetic polypharmacy risk not captured |
| AMP | Special Risk (congenital LQTS only) | Sympathomimetic polypharmacy risk not captured |
| ARI | Possible Risk (weakest category) | CardioSafe flags higher risk; FAERS ROR 8.15 with MPH |
| CLO | Not Classified | Conduction risk in polypharmacy not captured |
| GUA | Not Classified | Conduction risk; FAERS ROR 4.63 |

---

## Drug Architecture

### Three-pathway drug effect model

**Pathway 1 — hERG/IKr block**
Competitive binding: `block_i = C_free / (C_free + IC50_i)`. Combined: `1 - Π(1 - block_i)` (independent binding sites). Applied as GKr conductance scaling.

**Pathway 2 — Sympathomimetic** (MPH, AMP)
10% CL reduction (≈ +7 bpm, midpoint of published 3–8 bpm range) + 15% GCaL upregulation (PKA-mediated, Soltis & Saucerman 2010). ΔQTc predictions reflect Bazett-apparent change from HR elevation, not true repolarization prolongation.

**Pathway 3 — Autonomic modulation** (CLO, GUA)
10% CL increase (bradycardia) + 5% GNa reduction. Produces negative ΔQTc via Bazett correction. PR interval and AV conduction effects not captured.

---

## Methods Summary

| Component | Details |
|---|---|
| AP model | O'Hara-Rudy 2011, 41-variable ODE system |
| Integration | scipy.integrate.odeint, mxstep=5000, rtol=1e-6, atol=1e-8 |
| Simulation | 500 beats, CL=1000 ms (60 bpm); convergence confirmed at n>300 |
| APD90 | Duration from upstroke to 90% repolarization |
| QTc | Bazett: APD90 / sqrt(RR_s) |
| IC50 sources | ChEMBL v37 + Witchel 2002, Redfern 2003, Kongsamut 2002, Polak 2009, Kramer 2013, Perrin 2008 |
| Cmax values | Adult PK + protein binding; pediatric references in Supplementary S2 |
| FAERS | 40 quarters 2015Q1–2024Q4, pediatric filter age < 18, n=552,832 |
| ROR | 0.5 continuity correction, signal = CI lower bound > 1.0 |
| Concordance | AUC-ROC, Fisher exact, 10,000-iteration permutation test |
| PK sensitivity | 1x / 1.5x / 3x Cmax for 5 CYP2D6 substrates |
| CYP2D6 DDI | Ito static model, hepatic inlet [I]h |
| Composite score | 0.50×ΔQTc + 0.20×IKr + 0.30×FAERS ROR + flags, capped at 100 |

---

## ECG Calibration

Steady-state APD90 at 60 bpm (500 beats): **331.7 ms**. Systematic offset vs adolescent reference (Rijnbeek et al. 2014): **67 ± 9 ms**. Offset cancels exactly in ΔQTc calculations — confirmed numerically to floating-point precision.

---

## Known Limitations

- All Cmax values are adult-derived; stimulant values are conservative (children get higher weight-adjusted doses); CYP2D6 substrate uncertainty quantified in PK sensitivity analysis
- ORd model is adult ventricular; no validated adolescent-specific AP model exists
- Sympathomimetic parameters (10% CL, 15% GCaL) are phenomenological; full PKA cascade not implemented; parameter sensitivity confirms robustness
- No PR interval / AV conduction pathway; accounts for guanfacine FAERS discordances
- No base-case CYP2D6 PK interaction terms; partially addressed by Ito model analysis
- Active metabolites (norfluoxetine, 9-OH-risperidone) not modeled
- Sex differences in IKs not modeled

---

## Repository Structure

```
cardiosafe-pediatric/
├── src/
│   ├── ord_model.py                    # ORd AP model + drug effect architecture
│   ├── risk_grid.py                    # Full 84-combination polypharmacy sweep
│   ├── faers.py                        # FAERS download + parse + ROR pipeline
│   ├── faers_secondary.py              # Temporal trend + age stratification
│   ├── ecg_calibration.py              # APD90 vs pediatric ECG reference
│   ├── concordance_stats.py            # Model-FAERS classification metrics
│   ├── concordance_stats_stratified.py # Mechanism-stratified concordance
│   ├── concordance_nogap.py            # No-gap concordance + AUC-ROC
│   ├── pk_sensitivity.py               # Developmental PK sensitivity sweep
│   ├── cyp2d6_ito.py                   # CYP2D6 Ito static DDI model
│   ├── composite_score.py              # Composite cardiac risk score
│   ├── sympathomimetic_sensitivity.py  # Sympathomimetic parameter sensitivity
│   ├── generate_sim_data.py            # Auto-generate clinical simulator data
│   └── figures.py                      # Manuscript figures (Figures 1–5)
├── data/
│   ├── herg_master_params.csv
│   └── faers_cache/                    # gitignored, ~2GB
├── results/
│   ├── risk_grid_results.csv
│   ├── sympathomimetic_sensitivity.csv
│   ├── pk_sensitivity.csv
│   ├── cyp2d6_ito_results.csv
│   ├── composite_scores.csv
│   └── faers/
│       ├── faers_drug_ror.csv
│       ├── faers_combo_ror.csv
│       ├── concordance_statistics.csv
│       ├── temporal_trend.csv
│       └── age_stratification.csv
└── docs/
    ├── clinical_sim.html               # Patient-specific clinical risk simulator
    └── figures/                        # Figures 1–5 + Supplementary S1
```

---

## Installation

```bash
pip install numpy scipy pandas matplotlib seaborn scikit-learn tqdm pyarrow requests chembl-webresource-client
```

## Quick Start

```python
from src.ord_model import run_simulation

baseline = run_simulation(None, n_beats=500)
result = run_simulation(
    {"Methylphenidate": "therapeutic", "Aripiprazole": "therapeutic"},
    n_beats=500
)
print(f"ΔQTc: {result['QTc'] - baseline['QTc']:+.1f} ms")
```

---

## References

O'Hara T et al. *PLoS Comput Biol.* 2011;7(5):e1002061.
Dutta S et al. *Front Physiol.* 2017;8:616.
Rijnbeek PR et al. *J Electrocardiol.* 2014;47(6):914–921.
Redfern WS et al. *Cardiovasc Res.* 2003;58(1):32–45.
Soltis AR, Saucerman JJ. *Biophys J.* 2010;99(7):2038–2047.
Aman MG et al. *Clin Ther.* 2007;29:1476–86.
Findling RL et al. *J Clin Psychopharmacol.* 2008;28(4):441–446.
Templeton I et al. *Drug Metab Dispos.* 2016;44(1):57–65.

---

## Status

All analyses complete. Manuscript in preparation. Preprint submission pending.
Target journals: *npj Digital Medicine* · *Clinical Pharmacology & Therapeutics* · *JACAP*

## Author

Arnav Amit · Independent researcher · arnav.amit1@gmail.com