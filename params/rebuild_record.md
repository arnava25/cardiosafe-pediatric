# CardioSafe Pediatric: parameter table rebuild

**Record of work, 26 July 2026.** Written so the project state survives a gap.

This documents the rebuild of the hERG and pharmacokinetic parameter table that was invalidated by the June 2026 audit. It covers what is now sourced, what was wrong before and how, what is still open, and what has to be decided before the model is rerun.

Companion files: `herg_params_v2.csv` (the table), `build_params.py` (the gate).

---

## 1. Bottom line

The parameter foundation is rebuilt. Nine of twelve drugs now have a verified hERG IC50 traced to a primary patch clamp source with the assay method and recording temperature confirmed. Seven rows are model ready. Two rows are closed by bounding argument rather than measurement. Three drugs do not use this pathway.

**The result is not the null I expected.** Escitalopram sits at 12.5 percent hERG block in CYP2C19 normal metabolizers at the approved maximum pediatric dose, and 26.5 percent in poor metabolizers. Every other drug is under 6 percent. The mechanism is protein binding, not potency.

---

## 2. The table

Free Cmax computed in code as `total x (1 minus fraction bound)`. Fractional block computed as `free / (free + IC50)`.

| Drug | IC50 nM | Assay temp | f bound | Free Cmax nM | Block |
|---|---|---|---|---|---|
| Escitalopram | 700 | 37 C | 0.560 | 100.0 | **12.50%** |
| Quetiapine | 8300 | 22 to 24 C | 0.830 | 493.6 | 5.61% |
| Fluoxetine | 700 | 37 C | 0.945 | 22.0 | 3.05% |
| Risperidone | 160 | 22 to 24 C | 0.880 | 4.65 | 2.82% |
| Nortriptyline | 2200 | 23 C | 0.920 | 45.6 | 2.03% |
| Paliperidone | 570 | 22 to 24 C | 0.770 | 9.98 | 1.72% |
| Sertraline | 700 | 37 C | 0.980 | 10.8 | 1.52% |
| Norfluoxetine | 2500 | 37 C | not sourced | see 5.2 | excluded |
| Aripiprazole | not found | n/a | 0.990 | 5.13 | excluded |
| Imipramine | 3400 | ambient | not sourced | see 5.3 | excluded |
| Methylphenidate | n/a | n/a | n/a | n/a | sympathomimetic pathway |
| Amphetamine | n/a | n/a | n/a | n/a | sympathomimetic pathway |
| Clonidine | none published | n/a | n/a | n/a | autonomic pathway |
| Guanfacine | none published | n/a | n/a | n/a | autonomic pathway |

### 2.1 Escitalopram by CYP2C19 phenotype

Poweleit et al. 2023, measured pediatric population PK, normalized to 20 mg per day and BSA 1.73 square metres.

| Phenotype | Cmax ng/mL | Free nM | Block |
|---|---|---|---|
| Poor | 186.0 | 252.3 | 26.49% |
| Intermediate | 95.5 | 129.5 | 15.62% |
| Normal | 73.7 | 100.0 | 12.50% |
| Rapid | 70.1 | 95.1 | 11.96% |
| Ultrarapid | 65.0 | 88.2 | 11.19% |

Half life in poor metabolizers is 60.99 hours against 21.81 in normal metabolizers. Normalized trough is 151.1 ng/mL against 38.0, a 3.98 fold difference.

### 2.2 Risperidone by CYP2D6 phenotype

FDA Clinical Pharmacology Review NDA 20272/20588, allometrically scaled apparent clearance in L/h.

| | IM/PM | EM | ratio |
|---|---|---|---|
| Child, 39 kg, age 11 | 5.47 | 20.8 | 3.8x |
| Adolescent, 60 kg, age 15 | 7.56 | 28.7 | 3.8x |
| Adult, 70 kg | 8.48 | 32.2 | 3.8x |

Estimated proportion of the IM/PM subpopulation is 0.346.

---

## 3. What the old table got wrong, and how

Three distinct error modes. They need to be described separately because they are not the same failure.

**Fabricated attribution.** Aripiprazole was cited to Kramer 2013 and Perrin 2008. Neither paper contains an aripiprazole hERG measurement. The citation points at real papers that do not contain the drug.

**Citation slippage within a correct neighbourhood.** Imipramine was listed at 3388 nM and attributed to Witchel 2002. The value is essentially correct: Teschemacher et al. 1999 measured 3.4 micromolar. But Witchel 2002 is a citalopram and fluoxetine paper. Harry Witchel is the fourth author on Teschemacher 1999. Right number, right research group, wrong paper. This is the most insidious of the three because the number checks out.

**Wrong value and wrong citation.** Nortriptyline was listed at 1100 nM attributed to Witchel 2002. The real value is 2200 nM from Jeon et al. 2011, and Witchel 2002 contains no tricyclics.

Separately, quetiapine was listed at 1070 nM against its own cited source, Kongsamut 2002, which reports 5765 nM. And the entire free Cmax column was hand typed rather than computed, so the free fraction was never applied.

---

## 4. Methodological findings

These are contributions in their own right, independent of any model output.

### 4.1 The source literature is recorded at room temperature

Confirmed directly from Methods sections:

| Paper | Drug | Temperature |
|---|---|---|
| Lee 2017 | risperidone, paliperidone | room, 22 to 24 C |
| Lee 2018 | quetiapine, norquetiapine | room |
| Chae 2013 | escitalopram | room, 22 to 24 C |
| Jeon 2011 | nortriptyline | room |
| Kongsamut 2002 | antipsychotic series | room |
| Teschemacher 1999 | imipramine | ambient |
| Rajamani 2006 | fluoxetine, norfluoxetine | 23 C and 37 C |
| Lee 2012 (sertraline) | sertraline | 37 plus or minus 1 C |
| Zhang 2014 | escitalopram | physiological |

Five of nine at room temperature. The ORd model runs at 37 C.

**Temperature effect is drug specific and cannot be corrected with a factor.** Escitalopram measured at room temperature is 2.6 micromolar (Chae) and at 37 C is 0.70 micromolar (Zhang), a 3.7 fold difference in the same drug. Fluoxetine measured by Rajamani at both 23 C and 37 C gave 0.7 micromolar both times. So it matters enormously for one drug and not at all for another.

Where a 37 C value exists it should be used. Where it does not, the row carries an unquantified uncertainty of up to roughly fourfold in the direction of underestimating block.

### 4.2 Heterologous expression may also underestimate

Teschemacher note that imipramine completely blocks native cardiac IKr at 1 micromolar in guinea pig myocytes against their own heterologous value of 3.4, and that reduced sensitivity of heterologously expressed HERG relative to native IKr has been reported previously. This bias runs the same direction as temperature.

### 4.3 Point estimates misrepresent the real variability

Donovan et al. 2011 ran imipramine four times on one platform in one laboratory: 3.5, 6.5, 4.0 and 5.5 micromolar. A 1.9 fold spread within a single lab. They treat this as passing a three fold acceptance criterion, which appears to be the industry norm.

Across laboratories, fluoxetine spans 0.7 (Rajamani, 37 C), 1.5 (Witchel 2002), 2.1 (Donovan) and 3.1 (Thomas 2002). A 4.4 fold spread for one drug against one channel.

**Consequence for the rerun.** The 84 combination sweep should not be run at point estimates. It should be run as a sensitivity analysis across each drug's plausible IC50 range, reporting which conclusions survive. For most drugs this changes nothing because they sit far below threshold. For escitalopram and risperidone the band will straddle a threshold, and saying so is the correct result.

### 4.4 Total versus free concentration is handled inconsistently in the literature

Four papers, four degrees of engagement with the same arithmetic:

- **Rajamani 2006** applies it correctly and concludes therapeutic levels have minimal effect.
- **Teschemacher 1999** flags it explicitly as a caveat: interpretation must be qualified as the unbound concentration may not reflect drug levels in particular organ systems.
- **Kongsamut 2002** deliberately uses the ratio of total plasma concentration to IC50 and reports that it predicted QT prolongation.
- **Lee 2012 (sertraline)** compares total Cmax of 0.6 micromolar to their IC50 of 0.7, invokes the Redfern 30 fold margin, and concludes potentially harmful effects should be considered. At 98 percent protein binding the free margin is roughly 58 fold, comfortably inside Redfern.

The inconsistency, not the error, is the argument for a standardized sourced parameter table.

### 4.5 Two methods for closing rows where the data does not exist

**Required IC50 inversion.** Where free Cmax is known but IC50 is not, invert the Hill equation to find the IC50 the drug would need to reach a given block. Aripiprazole at the top of its pediatric therapeutic range would need an IC50 below 46 nM to reach 10 percent block, more potent than any drug in the set and approaching cisapride territory. Excluded without ever finding its IC50.

**Total concentration ceiling.** Where protein binding is unknown, check whether the total concentration alone can reach the required free concentration. Norfluoxetine at the measured pediatric mean has a total of 490 nM against a 625 nM free requirement for 20 percent block. Unreachable at zero binding, therefore unreachable. Excluded without ever finding its binding fraction.

Both belong in the Methods. Absent data is a recurring condition in pediatric psychopharmacology, not an accident.

---

## 5. Open items

### 5.1 Escitalopram Cmax, two sources disagree by twofold

Poweleit 2023 measured 38.0 ng/mL trough and 73.73 Cmax in pediatric normal metabolizers at 20 mg per day. Fekete 2020 gives a dose to concentration factor of 0.86 ng/mL per mg per day, implying roughly 17 to 20 ng/mL at the same dose.

| Anchor | Free Cmax | Block |
|---|---|---|
| Fekete implied, 20 ng/mL | 27 nM | 3.7% |
| Poweleit measured Cmax, 73.73 | 100 nM | 12.5% |

This determines whether escitalopram is the finding or is unremarkable. Poweleit is a direct population PK measurement in genotyped pediatric patients and should probably win, but the discrepancy needs resolving by reading what Fekete's factor is a factor for. **Highest priority open item.**

### 5.2 Norfluoxetine protein binding

No human value located. The only figure found is 85 to 90 percent for both fluoxetine and norfluoxetine in rat (Caccia 1990). Not used. Row closed by ceiling argument instead.

### 5.3 Imipramine protein binding and pediatric concentration

Neither sourced. Excluded regardless: reaching 20 percent block would require roughly 2383 ng/mL at 90 percent binding, ten to thirty times therapeutic.

### 5.4 Aripiprazole hERG IC50

Does not appear to exist. Multiple search passes. Row closed by inversion argument.

### 5.5 Awaiting reply

Email drafted to Sang June Hahn (sjhahn@catholic.ac.kr) asking for the recording temperature in Lee 2017, Lee 2018 and Chae 2013, and his view on the escitalopram discrepancy between his 2.6 micromolar and Zhang's 0.70. Temperature is now confirmed from the papers themselves, so the remaining value is his explanation of the discrepancy.

### 5.6 Not yet retrieved

Warrings W, Taurines R, Egberts K, Keicher F, Romanos M, Fekete S. Correlation Between Escitalopram, Sertraline, and Fluoxetine Serum Levels and QTc Interval Prolongation in Children and Adolescents. Ther Drug Monit 2026;48(3):366-372, PMID 41004670. 431 patients. This is a direct clinical test of what the table predicts, in the right population, with the right three drugs. Likely paywalled.

---

## 6. External corroboration

### 6.1 Faraj 2023, eBioMedicine

Human iPSC derived cardiomyocytes measuring APD90, APD30 and triangulation, cross referenced against a cohort of 19,742 patients on citalopram or escitalopram.

Three points of contact with this project:

**They use the same free fraction.** They tested 60 nM free and state this corresponds to 136 nM total at 56 percent protein binding. Free fraction 0.44, identical to the value used here, derived independently.

**They recommend a concentration ceiling of 100 nM total serum.** Pediatric exposure at the approved maximum dose exceeds it:

| | Total serum | vs ceiling |
|---|---|---|
| Faraj adult TDM average | 60 nM | 0.60x |
| Pediatric normal metabolizer trough, 20 mg | 117 nM | 1.17x |
| Pediatric normal metabolizer Cmax, 20 mg | 227 nM | 2.27x |
| Pediatric intermediate Cmax | 294 nM | 2.94x |
| Pediatric poor metabolizer Cmax | 573 nM | 5.73x |

Their concern population was patients over 65, where 20 percent on 10 mg exceed the threshold. The same comparison does not appear to have been made in children.

**They independently reject standard QT correction.** Their Methods state rate correction was based on an approximately linear relationship between spontaneous cycle length and APD90, and that this was significantly better at correcting for rate dependent changes than any of the existing QT correction factors. This is the June thesis, arrived at independently in human cardiomyocytes.

### 6.2 FDA quetiapine exposure response

FDA Clinical Pharmacology Review NDA 20639 SE5-045/046 modelled pediatric delta delta QTcF as 5.4 ms at 400 mg per day, 6.8 at 600 and 6.9 at 800, all under 10 ms. Observed mean delta QTcF across the two pivotal pediatric trials was about 2 ms, with no patient exceeding 500 ms or a 60 ms change.

This is a calibration target. Feeding 5.6 percent hERG block into the model should produce an APD90 change in that neighbourhood.

The same reviewers examined the QT to RR relationship in these specific pediatric patients and chose Fridericia over Bazett. A regulatory body, in the target population, declining to use Bazett.

### 6.3 FDA risperidone ECG findings

Across a 96 patient controlled study and a 79 patient open label extension: no significant mean changes in ECG parameters in any treatment group. Two subjects with prolonged QTc, both under 450 ms with changes under 60 ms. Mean pulse rate change was minus 2.4, plus 3.4 and plus 5.5 bpm across placebo, low and high dose.

---

## 7. Traps recorded

Near misses that would produce a plausible looking wrong number. Both are the same failure shape: right drug, real paper, real number, wrong current.

**Nortriptyline.** An IC50 of 2.86 micromolar exists for Kv currents in rabbit coronary arterial smooth muscle. Different channel, different tissue, different species. Numerically close to the real hERG value of 2.2 micromolar.

**Imipramine.** Most of the imipramine ion channel literature is hEAG, also called Eag1 or KV10.1, which is a different channel from hERG (KV11.1). Imipramine is a standard Eag1 tool compound in cancer proliferation research and dominates a naive search. A PubMed query for imipramine and ether-a-go-go returns mostly Eag1 papers.

---

## 8. Architectural gaps in the model

Not parameter problems. These are things the current code cannot represent.

**Metabolites are now first class rows.** Paliperidone and norfluoxetine have their own IC50s, binding fractions and measured pediatric concentrations. The model has no metabolite handling.

**Trafficking disruption.** Rajamani 2006 showed fluoxetine and norfluoxetine reduce hERG current by two mechanisms: direct pore block, and disruption of channel protein trafficking to the membrane. Trafficking dC50 values are 2.7 and 5.1 micromolar, close to their block IC50s. Mutating the pore binding site abolished block without affecting trafficking, so the sites are different. Trafficking develops over hours and recovers on washout. It is not a Hill equation phenomenon and cannot be modelled as fractional block.

**IKs block.** Three drugs in the set block IKs: sertraline at 12.3 to 15.2 micromolar, norfluoxetine at 5.3, and didesmethylcitalopram at 0.28. The model has an IKs upregulation term from the adrenergic work but no IKs block pathway. None reach threshold at therapeutic free concentrations, so this is completeness rather than missing risk.

**Stimulant pathway.** Still hardcoded and not concentration dependent. Unchanged since the June audit. Most of the 84 combinations involve methylphenidate or amphetamine, so this blocks the sweep.

**Fixed cycle length decomposition.** Still not implemented. Unchanged since the June audit.

**Independent binding site assumption.** Combined IKr block is computed as one minus the product of individual survivals. If two drugs share a hERG binding site this overestimates combined block. Untested, and it underlies every polypharmacy number in the project.

---

## 9. Decision pending

Framing. Deliberately deferred until the sweep is rerun, on the grounds that committing to a story before seeing results is what produced the inverted thesis in June.

What the evidence currently supports:

- A properly sourced pediatric hERG and PK parameter table for twelve drugs, which did not previously exist
- Escitalopram at 12.5 percent block in normal metabolizers and 26.5 in poor metabolizers, driven by protein binding rather than potency, in a drug carrying a pediatric indication down to age seven
- Two pharmacogenetic tails: CYP2C19 for escitalopram, CYP2D6 for risperidone, both around fourfold
- Near null for everything else at therapeutic exposure
- Four methodological findings, section 4
- Independent corroboration of the rate correction argument from Faraj

---

## 10. Files

| File | What it is |
|---|---|
| `herg_params_v2.csv` | The parameter table. Every value carries its own citation. Total Cmax and bound fraction stored separately; free Cmax never stored. |
| `build_params.py` | Gate. Refuses any row with a value but no source. Computes free Cmax and fractional block. Reports the required total concentration for 20 percent block on blocked rows. |
| `REBUILD_RECORD.md` | This file. |

Papers held in project context: kongsamut2002, lee2017, lee2018, chae2013, jeon2011 (nortriptyline), winter2008, aman2007, orsulak1988, Witchel 2002 (as 1s2_0S0014579301033208main), Rajamani 2006 (as 0706892a), Teschemacher 1999 (as 128-0702800a), donovan2011, Poweleit 2023 (as nihms1970824), fekete2020, Faraj 2023 (as PIIS2352396423003456), plus two risperidone FDA reviews.
