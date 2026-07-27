# results/ — stale as of the July 2026 parameter rebuild

**Every file in this directory except `results/faers/` predates the July 2026
parameter rebuild and must be regenerated before use.**

These outputs were computed against `data/herg_master_params.csv`, the hERG
parameter table withdrawn in the June 2026 citation audit. That table has been
moved to `archive/params_invalidated_202606/`; see the README there for the four
error modes found. Nothing here has been deleted, because the files are needed
for provenance and for before/after comparison once the rerun is done — but no
number in them should be quoted, plotted, or carried into a manuscript.

Affected files:

- `pairwise_results.csv`
- `risk_grid_results.csv`
- `risk_matrix.csv`
- `triple_results.csv`
- `rate_correction_comparison.csv`
- `supratherapeutic_sweep.csv`
- `pk_pediatric_genuine.csv`
- `risk_grid_500_newengine.log`
- `supratherapeutic.log`

Also superseded, and archived earlier under `archive/results/`: the CYP2D6 Ito
analysis (`cyp2d6_ito_results.csv`, `cyp2d6_ito_memo.md`). **This analysis is
parameter dependent, not parameter independent as previously assumed.**
`src/cyp2d6_ito.py` carries five free Cmax values hardcoded in its own source,
copied from the withdrawn table and annotated `# from herg_master_params.csv`.
Archiving the CSV did not remove them. The script is now guarded with a
`NotImplementedError` at import time and must be rewritten to read
`params/herg_params_v2.csv` before it can run.

## Exception: `results/faers/`

`results/faers/` is **parameter independent and remains valid.** Those outputs
are derived from FDA FAERS disproportionality analysis, not from the hERG
parameter table or the ORd simulation, so the audit does not touch them. They do
not need to be regenerated.

## Before regenerating

The rebuild record (`params/rebuild_record.md`) notes that the combination sweep
should not be rerun at point estimates. Measured IC50 values span up to roughly
fourfold across laboratories for a single drug, and five of nine source papers
recorded at room temperature rather than 37 C. The rerun should be a sensitivity
analysis across each drug's plausible IC50 range, reporting which conclusions
survive.
