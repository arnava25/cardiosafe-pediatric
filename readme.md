# CardioSafe Pediatric

A mechanistic model of cardiac electrophysiological risk from psychiatric polypharmacy in children and adolescents. Couples an O'Hara-Rudy 2011 ventricular action potential model with a sourced hERG and pharmacokinetic parameter set, and compares model output against pediatric adverse event reporting in FDA FAERS.

Solo authored. Arnav Amit, arnavamit1@g.ucla.edu.

---

## Status, July 2026

**The parameter table has been rebuilt and the simulation has not yet been rerun on it.**

An internal audit in June 2026 found citation and arithmetic errors in the original hERG parameter table that invalidated every quantitative model result. The submission to the Journal of Electrocardiology was withdrawn before submission, and the interactive risk simulator was taken down. The parameter table has since been reconstructed from primary sources with every value traced to a paper read in full. The combination sweep, the rate correction decomposition and the figures all still need regenerating on the new table.

The FAERS analysis is unaffected. It is parameter independent and remains valid.

Full account: [`params/rebuild_record.md`](params/rebuild_record.md).

---

## What the model does

Twelve drugs across the classes commonly co prescribed in pediatric psychiatry: two stimulants, three atypical antipsychotics, three SSRIs, two alpha-2 agonists, two tricyclics.

Three parallel drug effect pathways feed the action potential model. hERG block scales the rapid delayed rectifier conductance through the Hill equation. A sympathomimetic pathway shortens cycle length and raises L type calcium conductance. An alpha-2 pathway slows rate. Combinations are simulated to steady state at 500 beats, and action potential duration at 90 percent repolarization is extracted from the final beat. The frozen validated drug free baseline is 263.6 ms.

Two active metabolites are modelled as first class members rather than footnotes, because neither is prescribed and both arrive obligately with their parent: 9-hydroxyrisperidone with risperidone, and norfluoxetine with fluoxetine.

---

## Current finding

Six of twelve drugs have a complete sourced parameter set and can contribute a computed hERG block: escitalopram, fluoxetine, nortriptyline, quetiapine, risperidone, sertraline. Two more are excluded by bounding argument rather than measurement, described in section 4.5 of the rebuild record. Four do not use the hERG pathway.

Fractional hERG block at pediatric therapeutic exposure is under 6 percent for every drug except escitalopram.

**Escitalopram delivers two to three times more free drug to the channel in CYP2C19 poor metabolizers than in normal metabolizers.** That relative effect is robust across every published clearance estimate. The absolute block estimate is not: it spans 4.8 to 12.5 percent depending on which clearance value is used, because published escitalopram apparent clearance in this population ranges from 14.2 to 40 L/h. The mechanism is protein binding rather than potency. Escitalopram is 56 percent protein bound against sertraline's 98 percent, so at an identical measured IC50 of 700 nM it presents roughly twenty two times more free drug at the channel.

Model rerun pending. No tier assignment or combination result in this repository should be quoted until it lands.

---

## Repository status

| Path | Status |
|---|---|
| `params/` | **Current.** Rebuilt table, the gate that validates it, and the full record. |
| `src/ord_core.py` | **Current.** Validated ORd engine, parameter independent. |
| `src/ord_model.py` | **Current.** Drug module, rewired to `params/` July 2026. |
| `src/faers.py`, `src/faers_secondary.py` | **Current.** Parameter independent. |
| `data/faers_cache/` | **Current.** All 40 quarters, 2015q1 to 2024q4, verified. See `PROVENANCE.md`. |
| `results/faers/` | **Valid.** Parameter independent, does not need regenerating. |
| Other `src/*.py` | Code is sound. Consumes parameters through `ord_model.py`. Needs running, not rewriting. |
| Other `results/*` | **Stale.** Computed on the withdrawn table. See `results/README.md`. |
| `docs/figures/` | Four of five are stale. The age stratification figure is FAERS derived and valid. |
| `archive/` | Audit trail. The withdrawn table, the scripts that built it, the pulled submission, the withdrawn simulator. Each with a README explaining why. |

---

## The audit

Four distinct error modes were found in the original parameter table. They are worth separating, because only one of them is fabrication and the others are easier to make and harder to catch.

1. **Fabricated attribution.** Aripiprazole cited to two papers, neither of which contains an aripiprazole hERG measurement.
2. **Citation slippage within a correct neighbourhood.** Imipramine's value was correct but attributed to the wrong paper by the same author group. The number checks out, which is what makes it dangerous.
3. **Wrong value and wrong citation.** Nortriptyline listed at half its real value, attributed to a paper containing no tricyclics.
4. **Value duplicated across rows.** Two drugs with different molecular weights, protein binding and concentrations listed with an identical free Cmax. The only one of the four visible by eye.

Underneath all of them, a structural cause: free Cmax was stored as a hand typed constant rather than computed, so the protein binding fraction was never applied. And the parameter values existed in three independent places that never read from each other, which is why fixing one of them did not propagate.

The rebuilt table stores total concentration and bound fraction as separately sourced columns and computes free Cmax in code. `params/build_params.py` refuses any row carrying a value without a citation. `src/ord_model.py` raises rather than warning when the parameter table is missing or a drug is incompletely specified, because the previous version returned an empty dictionary and produced clean baseline numbers for drug combinations.

---

## Reproducing

```bash
python3 params/build_params.py     # validate the table, show computed block per drug
cd src && python3 ord_model.py     # confirm the parameter wiring and baseline
```

The FAERS pipeline needs the quarterly ASCII releases, which are not tracked. `data/faers_cache/PROVENANCE.md` records what was used and where to get it.

**Do not rerun the combination sweep at point estimates.** Measured hERG IC50 values vary 1.9 fold within a single laboratory and 4.4 fold across laboratories for the same drug, and five of nine source papers recorded at room temperature rather than 37 C. `ord_model.py` exposes `ic50_scale` and `cmax_scale` so the sweep can be run as a sensitivity analysis. Report which conclusions survive the range rather than a single tier.

Combinations containing a stimulant or an alpha-2 agonist return `uses_unrebuilt_pathway: True`. Those pathways are still hardcoded and not concentration dependent, which makes any result involving them a mathematical demonstration rather than a pharmacological finding. Rebuilding them is outstanding.

---

## Not currently available

The interactive risk simulator previously served from this repository has been withdrawn. It carried a full second copy of the invalidated parameter set as JavaScript constants, and its displayed risk tier was driven by a composite score that was removed from the manuscript as circular. Archived at `archive/simulator_withdrawn_202607/`.

---

## Known limitations

Beyond the pending rerun:

- **Independent binding site assumption, untested.** Combined hERG block is computed as one minus the product of individual survivals. Nearly all hERG blockers bind the same aromatic residues in the S6 pore cavity. If two drugs in a combination compete for that site, this overestimates combined block. This underlies every polypharmacy number in the project and has never been checked.
- **Trafficking disruption not modelled.** Fluoxetine and norfluoxetine reduce hERG current both by pore block and by disrupting channel trafficking to the membrane, at similar concentrations and through a different binding site. Trafficking develops over hours and is not a Hill equation phenomenon.
- **No IKs block pathway.** Three drugs in the set block IKs at measured concentrations. None reach threshold at therapeutic free exposure, so this is completeness rather than missing risk.
- **Fixed cycle length decomposition not implemented.** Genuine change in action potential duration still contains rate dependent restitution for rate changing drugs, so it does not yet isolate channel pharmacology.
- **Room temperature source data.** Five of nine hERG IC50 values were recorded at 22 to 24 C while the model runs at 37 C. The effect is drug specific and cannot be corrected with a factor: escitalopram differs 3.7 fold between temperatures, fluoxetine not at all.

---

This is a research model. It is not a validated clinical predictor and should not be used to make prescribing decisions.