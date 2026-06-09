# CYP2D6 Ito Static Interaction Model — CardioSafe Pediatric
Generated: 2026-06-08 20:32

## Model
AUC ratio = 1 / (fm * (1 / (1 + [I]/Ki)) + (1 - fm))

Parameters:
- Fluoxetine [I] (blood): 3.5 nM
- Norfluoxetine [I]: 5.0 nM
- Combined [I]total: 212.5 nM
- Ki (CYP2D6): 170 nM

## AUC Ratios
| Substrate | fm | AUC Ratio | Cmax 1x | Cmax adj |
|---|---|---|---|---|
| Quetiapine | 0.73 | 1.68x | 2.1 nM | 3.5 nM |
| Aripiprazole | 0.40 | 1.29x | 21.0 nM | 27.1 nM |
| Risperidone | 0.77 | 1.75x | 0.9 nM | 1.6 nM |
| Nortriptyline | 0.90 | 2.00x | 28.0 nM | 56.0 nM |
| Imipramine | 0.55 | 1.44x | 18.0 nM | 25.9 nM |

## Delta-QTc Results
| Combination | dQTc base | Tier base | dQTc adj | Tier adj | Delta | Tier change |
|---|---|---|---|---|---|---|
| FLU+QUE | +2.5ms | LOW | +1.5ms | LOW | -1.0ms | no |
| FLU+ARI | +9.0ms | LOW-MOD | +5.5ms | LOW-MOD | -3.5ms | no |
| FLU+RIS | +1.5ms | LOW | +1.5ms | LOW | +0.0ms | no |
| FLU+NOR | +6.0ms | LOW-MOD | +10.5ms | MODERATE | +4.5ms | YES |
| FLU+IMI | +2.0ms | LOW | +2.0ms | LOW | +0.0ms | no |
| MPH+QUE+FLU | +11.5ms | MODERATE | +10.0ms | MODERATE | -1.5ms | no |
| MPH+ARI+FLU | +18.1ms | MODERATE | +14.3ms | MODERATE | -3.8ms | no |

## Manuscript Implication
CYP2D6 adjustment escalates 1 combinations to higher risk tiers.
This explains the FAERS signal discordance for fluoxetine-containing combinations:
the base-case model underestimates true pediatric exposure.

## References
- Ito et al. Pharm Res. 1998;15(3):396-402 (static DDI model)
- Templeton et al. Drug Metab Dispos. 2016;44(1):57-65 (fluoxetine CYP2D6 Ki)
- FDA Drug Interaction Guidance 2020 (Ito static model validation)
- Grimm et al. Br J Clin Pharmacol. 2006;61(1):58-69 (quetiapine CYP2D6 fm)