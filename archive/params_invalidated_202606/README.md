# Invalidated hERG parameter table and its build chain

**Withdrawn June 2026. Retained for provenance. Do not use.**

These files produced the hERG / pharmacokinetic parameter table that was
withdrawn in June 2026 following an internal citation audit. Every quantitative
result downstream of this table is superseded.

## Contents

| File | Role |
|---|---|
| `pull.py` | ChEMBL pull, built `herg_full_data.csv` and `herg_summary.csv` |
| `pull2.py` | ChEMBL pull v2, built `herg_full_data_v2.csv` and `herg_summary_v2.csv` |
| `lit.py` | Manual literature table, built `herg_literature.csv` |
| `master.py` | Consolidated the above into `herg_master_params.csv` |
| `herg_full_data.csv`, `herg_full_data_v2.csv` | Raw ChEMBL pulls |
| `herg_literature.csv` | Hand-entered literature values |
| `herg_summary.csv`, `herg_summary_v2.csv` | Per-drug summaries of the pulls |
| `herg_master_params.csv` | The withdrawn parameter table itself |

## Why it was withdrawn

The audit found three distinct error modes in the citation layer. They are not
the same failure and are recorded separately.

**1. Fabricated citations.** Aripiprazole was attributed to Kramer 2013 and
Perrin 2008. Neither paper contains an aripiprazole hERG measurement. The
citations point at real papers that do not contain the drug.

**2. Citation slippage.** Imipramine's value was essentially correct but
attributed to the wrong paper by the same author group — listed against Witchel
2002, which is a citalopram and fluoxetine paper containing no tricyclics. The
number traces to Teschemacher et al. 1999, on which Harry Witchel is fourth
author. Right number, right research group, wrong paper. This is the most
insidious of the three because the value checks out on inspection.

**3. Wrong value plus wrong citation.** Nortriptyline was listed at 1100 nM
attributed to Witchel 2002. The real value is 2200 nM from Jeon et al. 2011.

Separately, the **free Cmax column was hand typed rather than computed**, so the
protein binding fraction was never applied to any row.

## Replacement

Replaced by `params/build_params.py`, which computes free Cmax in code as
`total x (1 - fraction bound)` and reads the rebuilt table `params/herg_params_v2.csv`.
See `params/rebuild_record.md` for the full rebuild record, source-by-source
verification, and the open questions carried forward.
