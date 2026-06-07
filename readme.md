# CardioSafe Pediatric

**Computational cardiac risk stratification for psychiatric polypharmacy in children and adolescents**

---

## Overview

Children and adolescents with psychiatric conditions are frequently prescribed multiple medications simultaneously. Stimulants, antipsychotics, antidepressants, and alpha-2 agonists are all commonly used, often in combination. These drugs carry real cardiac risks, but no rigorous computational framework exists for characterizing those risks in pediatric populations specifically.

CardioSafe Pediatric fills that gap. It couples O'Hara-Rudy (ORd) ventricular action potential modeling with hERG channel pharmacology and sympathomimetic pathway modeling to simulate QTc prolongation across clinically relevant polypharmacy combinations. The framework produces a mechanistic, drug-specific risk stratification output rather than a simple lookup table or empirical regression.

The central finding: the majority of high-risk combinations in this drug set are driven by sympathomimetic heart rate elevation rather than direct hERG block. This mechanism is entirely invisible to standard ion channel screening, which is the basis of most existing cardiac safety frameworks for psychiatric drugs.

---

## Background

### The clinical problem

Psychiatric polypharmacy in adolescents is common and growing. Stimulants (methylphenidate, amphetamine) are prescribed to 10-15% of school-age children. Atypical antipsychotics (risperidone, aripiprazole, quetiapine) are increasingly used off-label for irritability, aggression, and mood dysregulation. SSRIs are routinely used for pediatric anxiety and depression. These drugs are regularly co-prescribed.

Each drug class carries cardiac risk through different mechanisms:
- Antipsychotics and TCAs: hERG channel block, QTc prolongation
- SSRIs: mild hERG block, variable QTc effects
- Stimulants: sympathomimetic drive, heart rate and blood pressure elevation
- Alpha-2 agonists: bradycardia, PR prolongation, conduction slowing

No existing computational model characterizes the combined cardiac effects of these drug classes in combination, and no framework has been specifically built for adolescent pharmacokinetics.

### Why O'Hara-Rudy

The ORd 2011 model is the standard mechanistic model of the undiseased human ventricular action potential. It captures all major ionic currents (INa, ICaL, IKr, IKs, IK1, INaCa, INaK, and others) and their interactions. Because QTc prolongation is fundamentally an action potential duration (APD) phenomenon, mechanistic AP modeling is the appropriate tool for this problem. It allows drug effects to be applied at the level of specific ion channels and current pathways rather than inferred from population statistics.

---

## Drug block architecture

Three parallel drug effect pathways are modeled:

**1. hERG / IKr block** (antipsychotics, SSRIs, TCAs)

Competitive binding model:

```
block_i = C_free / (C_free + IC50_i)
```

Combined block across drugs (independent binding sites):

```
total_block = 1 - product(1 - block_i for each drug i)
```

Applied as IKr conductance scaling:

```
GKr_effective = GKr_max * (1 - total_block)
```

**2. Sympathomimetic pathway** (methylphenidate, amphetamine)

Modeled as heart rate increase (10% CL reduction) and beta-adrenergic ICaL upregulation (15% GCaL increase). Each stimulant in the combination contributes independently.

**3. Autonomic modulation** (clonidine, guanfacine)

Modeled as heart rate decrease (10% CL increase) and mild conduction slowing (5% GNa reduction). Produces negative delta-QTc via Bazett correction, consistent with the observed cardioprotective pattern in the results.

---

## Repository structure

```
cardiosafe/
├── ord_model.py              # O'Hara-Rudy AP model + drug block module + simulation runner
├── herg_master_params.csv    # Master drug parameterization table (IC50, Cmax_free, mechanism)
├── risk_grid.py              # Full polypharmacy combination sweep (all pairs + clinical triples)
├── herg_pull.py              # ChEMBL API pull for hERG IC50 data
├── herg_pull_v2.py           # Extended pull with alternate ChEMBL IDs
├── herg_literature.py        # Manually curated literature IC50 table with safety indices
├── herg_master_params.py     # Consolidation script: ChEMBL + literature -> master table
├── results/
│   ├── risk_grid_results.csv # All 84 combinations ranked by delta-QTc
│   ├── pairwise_results.csv  # Pairwise combinations only
│   ├── triple_results.csv    # Clinical triple combinations
│   └── risk_matrix.csv       # 12x12 pairwise delta-QTc matrix
└── report/
    └── cardiosafe_report.html # Interactive risk stratification report
```

---

## Results summary

84 drug combinations simulated (66 pairs + 18 clinical triples) across 12 drugs at therapeutic free plasma concentrations.

| Risk tier | Combinations | Delta-QTc threshold |
|-----------|-------------|---------------------|
| HIGH | 2 | >= 20 ms |
| MODERATE | 29 | 10-19 ms |
| LOW-MOD | 14 | 5-9 ms |
| LOW / protective | 39 | < 5 ms |

**Selected findings:**

- MPH + AMP: +23.8 ms, 0% IKr block (purely sympathomimetic)
- NOR + MPH + ARI: +21.7 ms, 6.76% IKr block (highest risk triple)
- MPH + ARI: +17.9 ms, 4.22% IKr block
- RIS + SER: +1.0 ms (low risk despite both being hERG blockers)
- ESC + CLO: -4.2 ms (alpha-2 protective effect)
- CLO + GUA: -6.1 ms (dual alpha-2 bradycardic effect)

28 of 31 HIGH or MODERATE combinations have IKr block < 5%, indicating sympathomimetic or autonomic mechanisms dominate over hERG block at therapeutic concentrations.

---

## Drug list

| Abbreviation | Drug | Class |
|---|---|---|
| MPH | Methylphenidate | Stimulant |
| AMP | Amphetamine | Stimulant |
| RIS | Risperidone | Antipsychotic |
| QUE | Quetiapine | Antipsychotic |
| ARI | Aripiprazole | Antipsychotic |
| SER | Sertraline | SSRI |
| FLU | Fluoxetine | SSRI |
| ESC | Escitalopram | SSRI |
| CLO | Clonidine | Alpha-2 agonist |
| GUA | Guanfacine | Alpha-2 agonist |
| IMI | Imipramine | TCA |
| NOR | Nortriptyline | TCA |

---

## Data sources

**hERG IC50 values:** ChEMBL v37 (target CHEMBL240), supplemented by published electrophysiology literature (Witchel 2002, Redfern 2003, Kongsamut 2002, Polak 2009, Kramer 2013, Perrin 2008). Patch clamp electrophysiology values preferred over radioligand binding assays.

**Free plasma Cmax:** Derived from published adult pharmacokinetic data and protein binding fractions.

**Data quality flags:**
- Quality A (patch clamp confirmed): IMI, NOR, RIS, QUE, ARI, SER, FLU
- Quality B (mixed/binding assay primary): ESC, CLO
- Quality C (estimate): GUA, MPH, AMP

---

## Limitations

- All free Cmax values are adult-derived. Adolescent pharmacokinetics differ due to developmental CYP2D6 variability, which metabolizes risperidone, aripiprazole, fluoxetine, nortriptyline, and imipramine. Pediatric Cmax may shift 1.5-3x.
- The ORd model represents the undiseased adult ventricular action potential. Adolescent baseline QTc is shorter by approximately 10-15 ms and autonomic tone differs developmentally. The model has not been validated against pediatric ECG data.
- Sympathomimetic and autonomic modifiers are phenomenological approximations. Full mechanistic PKA-mediated ICaL/IKs signaling cascades are not implemented.
- Active metabolite contributions are not included. Norfluoxetine and 9-OH-risperidone may contribute additional hERG block in vivo.
- Independent binding site assumption for combined IKr block may overestimate combined effect if drugs compete for the same binding site.
- 50-beat simulations may not represent true electrophysiological steady state. Full convergence requires 200+ beats; delta-QTc estimates carry approximately +/- 3-5 ms uncertainty from this source.
- Sex differences in adolescent QTc (hormonal modulation of IKs) are not modeled.

---

## Installation

```bash
pip install numpy scipy pandas matplotlib pypdf chembl-webresource-client
```

---

## Usage

Run baseline and a single drug combination:

```python
from ord_model import run_simulation

baseline = run_simulation(None, n_beats=50)

result = run_simulation(
    {"Methylphenidate": "therapeutic", "Aripiprazole": "therapeutic"},
    n_beats=50
)
print(f"Delta-QTc: {result['QTc'] - baseline['QTc']:+.1f} ms")
```

Run the full polypharmacy sweep:

```bash
python3 risk_grid.py
```

Pull fresh hERG data from ChEMBL:

```bash
python3 herg_pull_v2.py
```

---

## Reference

O'Hara T, Virag L, Varro A, Rudy Y. Simulation of the undiseased human cardiac ventricular action potential: model formulation and experimental validation. PLoS Comput Biol. 2011;7(5):e1002061.

---

## Status

In preparation. Manuscript in progress.

---

## Author

Arnav Amit  
Independent researcher  
arnav.amit1@gmail.com  
github.com/arnava25
