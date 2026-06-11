# Developmental PK Sensitivity Analysis — CardioSafe Pediatric
Generated: 2026-06-08 19:21
Beats per simulation: 200

## Rationale
All base-case simulations use adult-derived free Cmax values. CYP2D6 developmental
ontogeny and drug-drug pharmacokinetic interactions (fluoxetine inhibiting CYP2D6)
can increase exposure of substrate drugs by 1.5–3x in pediatric populations.
This analysis quantifies the impact on delta-QTc predictions.

## CYP2D6 Substrates Analyzed
- **Risperidone** (RIS): CYP2D6 primary metabolizer; children may have 1.5–2x adult Cmax
- **Aripiprazole** (ARI): CYP2D6 substrate; fluoxetine co-administration increases Cmax 2–4x
- **Fluoxetine** (FLU): CYP2D6 substrate AND inhibitor; active metabolite norfluoxetine adds to exposure
- **Nortriptyline** (NOR): CYP2D6 substrate; narrow therapeutic index; poor metabolizers at 3–10x Cmax
- **Imipramine** (IMI): CYP2D6/CYP2C19 substrate; variable pediatric exposure

## Key Findings

12 combination-multiplier combinations showed tier escalation:

| Combination | Substrate | Multiplier | 1x Tier | Scaled Tier | dQTc (ms) |
|---|---|---|---|---|---|
| MPH+ARI | ARI | x3.0 | MODERATE | HIGH | +33.3 |
| ARI+NOR | ARI | x3.0 | MODERATE | HIGH | +30.0 |
| QUE+ARI | ARI | x3.0 | MODERATE | HIGH | +26.5 |
| ARI+FLU | ARI | x3.0 | LOW-MOD | HIGH | +25.5 |
| ARI+SER | ARI | x3.0 | LOW-MOD | HIGH | +25.0 |
| ARI+NOR | NOR | x3.0 | MODERATE | HIGH | +24.0 |
| MPH+ARI | ARI | x1.5 | MODERATE | HIGH | +21.4 |
| ARI+SER | ARI | x1.5 | LOW-MOD | MODERATE | +13.0 |
| ARI+FLU | ARI | x1.5 | LOW-MOD | MODERATE | +13.0 |
| MPH+RIS | RIS | x3.0 | LOW-MOD | MODERATE | +11.0 |
| ARI+FLU | FLU | x3.0 | LOW-MOD | MODERATE | +10.5 |
| MPH+RIS | RIS | x1.5 | LOW-MOD | MODERATE | +10.0 |

## Manuscript Implication
The base-case delta-QTc predictions represent a conservative lower bound for
CYP2D6 substrate combinations in pediatric populations. At clinically plausible
exposure multipliers, several LOW or LOW-MOD combinations may approach or exceed
the MODERATE threshold, particularly for combinations involving fluoxetine
(CYP2D6 inhibitor) co-prescribed with quetiapine or aripiprazole (CYP2D6 substrates).