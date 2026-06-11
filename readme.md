# CardioSafe Pediatric

**Computational cardiac risk stratification for psychiatric polypharmacy in children and adolescents**

[![status](https://img.shields.io/badge/status-preprint%20in%20preparation-blue)](https://github.com/arnava25/cardiosafe-pediatric)
[![model](https://img.shields.io/badge/model-O'Hara--Rudy%202011-navy)](https://doi.org/10.1371/journal.pcbi.1002061)
[![validation](https://img.shields.io/badge/validation-FDA%20FAERS%20pediatric-green)](https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html)

---

## Overview

Children and adolescents with psychiatric conditions are frequently prescribed several medications at once: stimulants, antipsychotics, antidepressants, and alpha-2 agonists, often in combination. These drugs carry real cardiac liabilities, and standard safety screening (CiPA, CredibleMeds) classifies them almost entirely by hERG channel block.

CardioSafe Pediatric couples a validated **O'Hara-Rudy (ORd) 2011 ventricular action potential model** with a three-pathway psychiatric-drug effect module (hERG/IKr block, sympathomimetic heart-rate elevation, alpha-2 autonomic modulation) across 12 drugs and 84 polypharmacy combinations, and cross-references the predictions against **FDA FAERS pharmacovigilance data** (2015 to 2024, roughly 16 million reports, about 552,000 pediatric cases).

> This is a hypothesis-generating computational framework, not a validated clinical guideline. It is intended to identify combinations and questions worth prospective study, nothing more.

---

## Central Finding

**Most of the apparent QTc prolongation in stimulant-containing polypharmacy is a heart-rate-correction artifact, not genuine action potential prolongation.**

When the model output is separated into the genuine change in action potential duration (ΔAPD90) versus the Bazett rate-corrected ΔQTc, the two diverge sharply and are negatively correlated across the 66 drug pairs (r approximately −0.34). Stimulants raise heart rate, which inflates Bazett-corrected QTc even when the underlying action potential barely changes or shortens.

| Combination | Bazett ΔQTc | Genuine ΔAPD90 | IKr block | Reading |
|---|---|---|---|---|
| MPH + AMP | +20.0 ms | **−8.4 ms** | 0.00% | Pure rate artifact; net genuine shortening |
| ARI + MPH | +18.3 ms | +3.8 ms | 4.22% | Mostly rate artifact, small genuine change |
| MPH + NOR | +15.4 ms | +1.1 ms | 2.65% | Mostly rate artifact |
| MPH + QUE | +12.6 ms | −1.6 ms | 0.93% | Rate artifact; net shortening |
| MPH + SER | +11.6 ms | −2.5 ms | 0.19% | Rate artifact; net shortening |
| **ARI + NOR** | **+11.5 ms** | **+11.5 ms** | 6.76% | Genuine, hERG-mediated prolongation |
| CLO + GUA | −14.4 ms | **+10.6 ms** | 0.00% | Bradycardia masks genuine prolongation |
| RIS + SER | +0.5 ms | +0.5 ms | 0.55% | Minimal either way |

The implication runs in both directions. Bazett over-calls stimulant pairs, and it under-calls (even reverses the sign of) the alpha-2 combinations, where slowing the heart rate drives the corrected QTc negative while the genuine action potential lengthens.

### Over-flagging

Tiering on Bazett ΔQTc flags **29 of 84 combinations** at 10 ms or above. Tiering on genuine ΔAPD90 flags only **5**. Standard QTc-based screening would over-call this drug set by roughly 6x, almost all of it stimulant-driven rate artifact.

### Where genuine prolongation does occur

It is small and concentrated in hERG-blocking antipsychotic and tricyclic combinations (ARI + NOR, +11.5 ms genuine) and, through bradycardia-driven restitution, in some alpha-2 combinations. No pairwise or triple combination reaches a genuine ΔAPD90 of 20 ms.

Triples show the same split:

| Combination | Bazett ΔQTc | Genuine ΔAPD90 |
|---|---|---|
| NOR + MPH + ARI | +22.5 ms | +7.9 ms |
| MPH + RIS + SER | +12.1 ms | −2.1 ms |
| ARI + FLU + CLO | −1.8 ms | +11.0 ms |

ARI + FLU + CLO is the cleanest inversion: Bazett reads it as a non-event, the genuine action potential change is the largest of any triple.

---

## FAERS Pharmacovigilance

The FAERS analysis is independent of the action potential model and stands on its own.

FDA FAERS 2015Q1 to 2024Q4, roughly 16 million reports, about 552,000 pediatric cases (age under 18). Reporting odds ratios computed with a 0.5 continuity correction; a signal is the lower 95% confidence bound exceeding 1.0.

| Combination | n | FAERS ROR | Genuine ΔAPD90 |
|---|---|---|---|
| ESC + IMI | — | 28.01 | +1.0 ms |
| IMI + SER | — | 16.22 | +1.5 ms |
| MPH + SER | 578 | 12.79 [9.84 to 16.62] | −2.5 ms |
| ARI + GUA | 359 | 10.25 [7.14 to 14.71] | +10.4 ms |
| MPH + ARI | 677 | 8.15 [6.10 to 10.90] | +3.8 ms |
| MPH + QUE | 274 | 5.25 [3.04 to 9.08] | −1.6 ms |

### Age stratification (most robust empirical result)

Stimulant-containing cardiac signals are far stronger in children 6 to 12 than in adolescents 13 to 17. CYP2D6-substrate combinations show the opposite direction, consistent with hepatic enzyme ontogeny.

| Combination | Children 6 to 12 | Adolescents 13 to 17 | Direction |
|---|---|---|---|
| MPH + SER | 43.71 | 1.86 | ~23x higher in children |
| MPH + ARI | 23.71 | 0.47 | ~50x higher in children |
| QUE + GUA | 32.15 | 3.91 | Higher in children |
| QUE + FLU | 3.75 | 10.88 | Higher in adolescents |

### Concordance, stated honestly

Model-FAERS concordance is **weak overall** (AUC approximately 0.55). The model has reasonable sensitivity for stimulant-containing combinations and near-zero sensitivity for non-stimulant ones. Even the stimulant concordance is confounded: both the model's Bazett ΔQTc and FAERS reporting frequency track "a stimulant is present," so agreement there does not establish a shared mechanism. The previously reported composite score that combined ΔQTc, IKr, and FAERS ROR into a single index was circular (FAERS appeared on both sides of the comparison) and has been removed from the codebase.

---

## Drug Architecture

**Pathway 1, hERG/IKr block.** Competitive binding, `block_i = C_free / (C_free + IC50_i)`, combined across drugs as `1 - product(1 - block_i)` assuming independent sites, applied as GKr conductance scaling.

**Pathway 2, sympathomimetic (MPH, AMP).** Heart-rate elevation modeled as a reduction in cycle length (about +7 bpm, midpoint of the published 3 to 8 bpm range). The genuine action potential effect at the cell level is minimal; the apparent QTc rise is overwhelmingly a Bazett rate-correction artifact.

**Pathway 3, autonomic modulation (CLO, GUA).** Heart-rate reduction (bradycardia) plus a small sodium-conductance reduction. Drives Bazett ΔQTc negative while the genuine action potential can lengthen. PR interval and AV conduction effects are not represented.

---

## Methods Summary

| Component | Details |
|---|---|
| AP model | O'Hara-Rudy 2011 endocardial cell, validated against a frozen golden baseline |
| Engine | `src/ord_core.py`, canonical ORd; drug module in `src/ord_model.py` |
| Integration | scipy.integrate.odeint, mxstep=5000, rtol=1e-6, atol=1e-8 |
| Simulation | 500 beats, CL = 1000 ms (60 bpm); ionic steady state reached around beat 100 |
| Baseline | APD90 = 263.6 ms at 60 bpm (Bazett QTc identical at RR = 1 s) |
| APD90 | Absolute time from upstroke to 90% repolarization |
| Rate correction | Bazett ΔQTc reported alongside genuine ΔAPD90 for every combination |
| IC50 sources | ChEMBL where indexed, manually curated literature hERG values otherwise |
| FAERS | 40 quarters, 2015Q1 to 2024Q4, pediatric filter age under 18 |
| ROR | 0.5 continuity correction, signal = lower 95% bound above 1.0 |

### ECG calibration

The model baseline carries a fixed additive offset relative to the Rijnbeek 2014 adolescent reference. Because the offset is a constant, it cancels exactly in any ΔQTc or ΔAPD90 comparison; the comparative findings above are offset-invariant. The absolute calibration figure is reported in `results/`.

---

## Known Limitations

- The ORd model is an adult ventricular cell; no validated adolescent-specific action potential model exists.
- Cmax values are largely adult-derived. Stimulant values are conservative, since children receive higher weight-adjusted doses.
- The sympathomimetic and autonomic pathways are phenomenological cycle-length and conductance adjustments, not a full autonomic or PKA cascade.
- No PR interval or AV conduction pathway, which is why several guanfacine and clonidine FAERS signals are not reproduced by the model.
- Active metabolites (norfluoxetine, 9-OH-risperidone) and sex differences in IKs are not modeled.
- Dose-scaling, pediatric-Cmax, and CYP2D6 adjustments in the interactive simulator are heuristic extrapolations of the base grid values, not separately simulated points.

---

## Repository Structure

```
cardiosafe-pediatric/
├── src/
│   ├── ord_core.py             # Validated canonical ORd engine (golden baseline)
│   ├── ord_model.py            # Three-pathway drug effect module
│   ├── risk_grid.py            # 84-combination polypharmacy sweep
│   ├── rate_correction.py      # Genuine ΔAPD90 vs Bazett ΔQTc decomposition
│   ├── concordance_metrics.py  # Model vs FAERS classification metrics
│   ├── supratherapeutic_sweep.py
│   ├── faers.py                # FAERS download, parse, ROR pipeline
│   ├── faers_secondary.py      # Temporal trend + age stratification
│   ├── generate_sim_data.py    # Regenerates clinical_sim.html data from the grid
│   └── figures.py              # Manuscript figures
├── data/
│   ├── herg_master_params.csv
│   └── faers_cache/            # gitignored
├── results/
│   ├── risk_grid_results.csv
│   ├── pairwise_results.csv
│   ├── triple_results.csv
│   ├── rate_correction_comparison.csv
│   └── faers/
│       ├── faers_drug_ror.csv
│       ├── faers_combo_ror.csv
│       ├── temporal_trend.csv
│       └── age_stratification.csv
├── docs/
│   ├── clinical_sim.html       # Personal exploration tool (not a validated instrument)
│   └── figures/
└── archive/                    # Superseded artifacts from the pre-validation engine
```

---

## Installation

```bash
pip install numpy scipy pandas matplotlib seaborn scikit-learn tqdm pyarrow requests chembl-webresource-client
```

## Reproducing the results

```bash
python src/risk_grid.py            # regenerates results/risk_grid_results.csv
python src/rate_correction.py      # genuine vs Bazett decomposition
python src/faers.py --analyze      # FAERS ROR pipeline
python src/concordance_metrics.py  # model vs FAERS concordance
```

---

## References

O'Hara T et al. *PLoS Comput Biol.* 2011;7(5):e1002061.
Dutta S et al. *Front Physiol.* 2017;8:616.
Rijnbeek PR et al. *J Electrocardiol.* 2014;47(6):914 to 921.
Redfern WS et al. *Cardiovasc Res.* 2003;58(1):32 to 45.
Aman MG et al. *Clin Ther.* 2007;29:1476 to 1486.
Templeton I et al. *Drug Metab Dispos.* 2016;44(1):57 to 65.

---

## Status

Rebuilt on the validated ORd engine. Core grid, rate-correction decomposition, and FAERS pipeline complete; manuscript being updated to the rate-correction framing; supratherapeutic sweep in progress.

## Author

Arnav Amit, independent researcher, arnav.amit1@gmail.com