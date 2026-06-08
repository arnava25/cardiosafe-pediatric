# FAERS Pharmacovigilance Analysis — Results & Interpretation
**CardioSafe Pediatric | Analysis date: June 2026**

---

## Dataset

- Source: FDA FAERS ASCII quarterly files, 2015Q1–2024Q4
- Total cases parsed: 16,144,530 (after deduplication)
- Pediatric cases (age < 18): 552,832
- Age unknown/excluded: 6,980,255
- Drug rows matched to 12-drug list: 1,349,344
- Pediatric cases with primary cardiac PT: 5,331
- Pediatric cases with any cardiac PT (broad set): 10,606
- All 12 drugs had at least one pediatric report

Primary cardiac PTs: QT prolonged, Electrocardiogram QT prolonged, Torsade de
pointes, Ventricular tachycardia, Ventricular fibrillation, Cardiac arrest,
Sudden cardiac death, Long QT syndrome, Ventricular arrhythmia.

---

## Per-Drug Reporting Odds Ratios (Pediatric)

| Drug | n | Cardiac events | ROR | 95% CI | Signal |
|---|---|---|---|---|---|
| Nortriptyline | 200 | 18 | 10.44 | 6.47–16.86 | YES |
| Guanfacine | 2,936 | 124 | 4.63 | 3.86–5.55 | YES |
| Imipramine | 172 | 7 | 4.66 | 2.24–9.69 | YES |
| Quetiapine | 4,936 | 189 | 4.21 | 3.63–4.88 | YES |
| Fluoxetine | 5,894 | 218 | 4.08 | 3.55–4.68 | YES |
| Sertraline | 7,451 | 215 | 3.15 | 2.74–3.61 | YES |
| Escitalopram | 2,658 | 60 | 2.41 | 1.86–3.11 | YES |
| Aripiprazole | 6,701 | 170 | 2.74 | 2.34–3.19 | YES |
| Methylphenidate | 10,523 | 169 | 1.70 | 1.46–1.99 | YES |
| Amphetamine | 5,537 | 70 | 1.33 | 1.05–1.68 | YES |
| Risperidone | 14,000 | 159 | 1.19 | 1.01–1.39 | YES |
| Clonidine | 4,306 | 50 | 1.22 | 0.92–1.61 | no |

11/12 drugs show pharmacovigilance signal. Clonidine is the only drug without
significant signal, consistent with its modeled cardioprotective effect
(negative delta-QTc via alpha-2 agonism and bradycardia).

Nortriptyline has the highest per-drug ROR (10.44), consistent with known TCA
hERG block potency. Guanfacine ROR=4.63 is elevated despite the model
predicting protective QTc effects — this likely reflects PR prolongation and
AV conduction slowing, which the model does not capture (see discordance
analysis below).

---

## Pairwise Combination ROR — Top Signals

| Combination | n | ROR | 95% CI | Signal |
|---|---|---|---|---|
| Escitalopram + Imipramine | 6 | 28.01 | 4.60–170.42 | YES |
| Imipramine + Sertraline | 10 | 16.22 | 2.90–90.77 | YES |
| Aripiprazole + Imipramine | 11 | 14.67 | 2.65–81.21 | YES |
| Guanfacine + Quetiapine | 127 | 13.15 | 7.61–22.73 | YES |
| Methylphenidate + Sertraline | 578 | 12.79 | 9.84–16.62 | YES |
| Aripiprazole + Guanfacine | 359 | 10.25 | 7.14–14.71 | YES |
| Methylphenidate + Aripiprazole | — | 8.15 | 6.10–10.90 | YES |
| Aripiprazole + Sertraline | — | 6.78 | 5.24–8.77 | YES |
| Quetiapine + Fluoxetine | — | 9.18 | 6.75–12.49 | YES |
| Quetiapine + Sertraline | — | 5.84 | 4.32–7.91 | YES |

17/63 evaluable pairs showed pharmacovigilance signal (CI lower bound > 1).

---

## Model vs. FAERS Alignment

Overall concordance: **39/63 pairs (62%)**

### Concordant — model correctly predicted elevated risk

| Combination | Model dQTc | Tier | FAERS ROR |
|---|---|---|---|
| MPH + ARI | +17.9 ms | MODERATE | 8.15 [6.10–10.90] |
| MPH + QUE | +12.7 ms | MODERATE | 5.25 [3.04–9.08] |
| MPH + SER | +11.3 ms | MODERATE | 12.79 [9.84–16.62] |
| AMP + SER | +11.3 ms | MODERATE | 2.81 [1.56–5.05] |

MPH+SER is the strongest single validation point: n=578, tight CI, ROR=12.79.
This is one of the most commonly prescribed combinations in adolescent
psychiatry (ADHD + depression/anxiety). The model correctly flagged it as
MODERATE risk (+11.3 ms) via sympathomimetic mechanism, and FAERS confirms a
large real-world cardiac AE signal.

### Discordant — model underestimated risk (FAERS signal, model LOW/LOW-MOD)

| Combination | Model dQTc | Tier | FAERS ROR | Likely explanation |
|---|---|---|---|---|
| ARI + GUA | +2.7 ms | LOW | 10.25 [7.14–14.71] | PR/conduction pathway absent |
| QUE + GUA | -2.6 ms | LOW | 13.15 [7.61–22.73] | PR/conduction pathway absent |
| QUE + FLU | +2.0 ms | LOW | 9.18 [6.75–12.49] | CYP2D6 PK interaction absent |
| QUE + SER | +2.0 ms | LOW | 5.84 [4.32–7.91] | CYP inhibition |
| ARI + SER | +7.0 ms | LOW-MOD | 6.78 [5.24–8.77] | CYP2D6 PK interaction |
| SER + IMI | +1.5 ms | LOW | 16.22 [2.90–90.77] | TCA narrow margin + PK (n=10) |
| ESC + IMI | +1.0 ms | LOW | 28.01 [4.60–170.42] | TCA + SSRI serotonin + cardiac (n=6) |

### Discordant — model predicted elevated risk, FAERS no signal

| Combination | Model dQTc | Tier | FAERS ROR | Note |
|---|---|---|---|---|
| MPH + AMP | +23.8 ms | HIGH | 0.45 [0.13–1.54] | Dual stimulant co-Rx rare; n very low |
| MPH + RIS | +11.7 ms | MODERATE | 0.65 [0.28–1.50] | Risperidone has low per-drug ROR (1.19) |
| MPH + FLU | +11.7 ms | MODERATE | 0.12 [0.01–1.86] | Very low n; unstable estimate |

---

## Mechanistic Interpretation of Discordances

### 1. Guanfacine combinations — PR/conduction pathway gap

The model treats guanfacine as cardioprotective: alpha-2 agonism reduces heart
rate and produces negative delta-QTc via Bazett correction. This is
mechanistically correct for QTc specifically.

However, guanfacine produces significant AV node slowing and PR prolongation
through Gi-coupled signaling in sinoatrial and atrioventricular nodal tissue.
When combined with antipsychotics (which slow conduction via sodium channel
effects) or SSRIs (which have indirect autonomic effects), the combined
conduction slowing may produce clinically significant bradyarrhythmia or
heart block — events that appear in FAERS under cardiac PTs but are not
captured by QTc-based risk scoring.

The model has no PR interval or conduction velocity component. This is a
primary limitation and a clear future direction: adding a PR/AV nodal
conduction pathway to the risk framework.

### 2. Fluoxetine/quetiapine and fluoxetine/aripiprazole — CYP2D6 inhibition

Fluoxetine (and its active metabolite norfluoxetine) is a potent CYP2D6
inhibitor. Quetiapine and aripiprazole are both CYP2D6 substrates. Co-
administration of fluoxetine significantly increases quetiapine and
aripiprazole plasma concentrations — by 2–4x in some reports — well above
the therapeutic free Cmax values used in the model parameterization.

The model uses fixed adult Cmax values with no PK interaction term. If true
plasma concentrations are 2–4x higher due to enzyme inhibition, the effective
hERG block and cardiac exposure are substantially underestimated for these
combinations. This explains why QUE+FLU (model: +2.0 ms LOW) and ARI+SER
(model: +7.0 ms LOW-MOD) show strong FAERS signal despite low predicted
delta-QTc.

Incorporating CYP2D6-mediated PK interactions is the single highest-impact
modeling improvement for these drug classes.

### 3. MPH+AMP — mechanism-detection gap

The model's highest-risk prediction (MPH+AMP, +23.8 ms HIGH) shows no FAERS
signal. Two explanations:

First, dual stimulant co-prescription is clinically rare and generally
contraindicated, so FAERS n is very low and the ROR estimate (0.45) is
unstable noise rather than a true negative.

Second, the sympathomimetic QTc mechanism may not translate directly to the
MedDRA cardiac PTs used in the FAERS filter. Heart rate elevation and
adrenergic-mediated QTc prolongation may present clinically as palpitations,
chest pain, or anxiety rather than triggering specific QT/arrhythmia codes.
This is a known limitation of FAERS for non-channel cardiac mechanisms.

---

## Cardiac PT Breakdown (High-Risk Drugs)

For cases involving Methylphenidate, Amphetamine, Aripiprazole, Risperidone,
Imipramine, or Nortriptyline:

| PT | Count |
|---|---|
| Electrocardiogram QT prolonged | 297 |
| Heart rate increased | 230 |
| Cardiac arrest | 173 |
| Palpitations | 157 |
| Syncope | 148 |
| Electrocardiogram abnormal | 47 |
| Ventricular fibrillation | 25 |
| Ventricular tachycardia | 20 |
| Torsade de pointes | 20 |
| Long QT syndrome | 15 |
| Sudden cardiac death | 9 |
| Atrioventricular block | 8 |
| Ventricular arrhythmia | 7 |

"Heart rate increased" (n=230) is the second most common PT, consistent with
the sympathomimetic mechanism identified in the model. This supports the
framework's core finding that HR elevation is a primary cardiac risk pathway
for stimulant-containing combinations, not captured by hERG-only screening.

---

## Limitations of This FAERS Analysis

- Age unknown for 43% of FAERS cases, excluded conservatively. True pediatric
  n may be higher; unknown-age cases may be enriched for pediatric patients
  in some reporting contexts.
- FAERS captures reported adverse events, not confirmed diagnoses. QTc
  prolongation may be underreported relative to clinical incidence.
- Confounding by indication: drugs co-prescribed for more severe psychiatric
  presentations may carry higher baseline cardiac risk independent of
  pharmacology.
- ROR is not an incidence estimate. It measures disproportionate reporting
  relative to the full FAERS database background, not absolute risk.
- Duplicate reports are a known FAERS issue. Deduplication by primaryid
  was applied; some duplicates from re-submitted or follow-up reports may
  remain.
- TCA combinations (imipramine, nortriptyline) have low n due to declining
  use. Wide confidence intervals limit interpretability.

---

## Implications for Manuscript

The FAERS analysis provides convergent pharmacovigilance evidence supporting
the model's core predictions for stimulant-antipsychotic combinations, while
identifying two specific mechanistic gaps — conduction pathway modeling and
CYP-mediated PK interactions — that define the next iteration of the
framework.

The MPH+SER finding (ROR=12.79, n=578) is the strongest single result: a
highly prevalent clinical combination, correctly flagged as moderate-risk by
the model via sympathomimetic mechanism, with a large and precisely estimated
real-world cardiac AE signal. This combination is a direct candidate for
prospective ECG monitoring guidance.

The systematic guanfacine discordance (model predicts protective, FAERS shows
risk) is the most scientifically interesting finding and should be highlighted
as a hypothesis-generating result: guanfacine's cardiac risk profile in
polypharmacy contexts may be driven by conduction effects rather than QTc,
and warrants dedicated investigation.
