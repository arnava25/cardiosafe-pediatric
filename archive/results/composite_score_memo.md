# Composite Risk Score — CardioSafe Pediatric
Generated: 2026-06-08 20:39

## Formula
score = 0.5 * dQTc_component + 0.2 * IKr_component + 0.3 * FAERS_component
      + CYP2D6_flag(15pts) + conduction_flag(10pts)
      capped at 100

## Component Scaling
- dQTc: (dQTc / 20.0ms) * 100, floor 0
- IKr: (IKr_pct / 10.0%) * 100
- FAERS: (log(ROR) / log(30.0)) * 100, 0 if no signal

## Score Labels
75-100: HIGH | 50-74: MODERATE | 25-49: LOW-MOD | 0-24: LOW

## Distribution
- HIGH: 0 combinations
- MODERATE: 5 combinations
- LOW-MOD: 30 combinations
- LOW: 31 combinations

## Top 20 Combinations
| Combination | Score | Label | dQTc | IKr% | ROR | Flags |
|---|---|---|---|---|---|---|
| MPH+ARI | 72.2 | MODERATE | +18.1 | 4.22% | 8.15 | none |
| ARI+IMI | 61.0 | MODERATE | +11.0 | 4.92% | 14.67 | none |
| ARI+GUA | 57.2 | MODERATE | +7.3 | 4.22% | 10.25 | COND |
| AMP+ARI | 53.7 | MODERATE | +18.1 | 4.22% | 0.69 | none |
| ARI+NOR | 51.0 | MODERATE | +15.0 | 6.76% | 9.34 | none |
| ARI+SER | 49.4 | LOW-MOD | +9.5 | 4.40% | 6.78 | none |
| ARI+FLU | 49.1 | LOW-MOD | +10.0 | 4.53% | 1.31 | CYP2D6 |
| MPH+AMP | 47.8 | LOW-MOD | +19.1 | 0.00% | 0.45 | none |
| MPH+SER | 46.9 | LOW-MOD | +9.6 | 0.19% | 12.79 | none |
| QUE+FLU | 44.6 | LOW-MOD | +3.0 | 1.25% | 9.18 | CYP2D6 |
| MPH+QUE | 44.0 | LOW-MOD | +11.0 | 0.93% | 5.25 | none |
| MPH+NOR | 42.3 | LOW-MOD | +14.8 | 2.65% | 1.00 | none |
| AMP+NOR | 42.3 | LOW-MOD | +14.8 | 2.65% | 11.41 | none |
| ARI+ESC | 39.7 | LOW-MOD | +9.5 | 4.25% | 2.31 | none |
| QUE+ARI | 39.0 | LOW-MOD | +11.5 | 5.11% | 1.60 | none |
| RIS+ARI | 38.8 | LOW-MOD | +10.0 | 4.56% | 1.71 | none |
| MPH+GUA | 37.8 | LOW-MOD | +6.8 | 0.00% | 3.40 | COND |
| FLU+NOR | 37.2 | LOW-MOD | +6.5 | 2.98% | 9.34 | CYP2D6 |
| ESC+IMI | 35.9 | LOW-MOD | +2.0 | 0.76% | 28.01 | none |
| QUE+GUA | 34.6 | LOW-MOD | -0.1 | 0.93% | 13.15 | COND |

## Manuscript Note
The composite score integrates mechanistic and epidemiological evidence.
Key property: MPH+SER scores higher than its delta-QTc alone would suggest,
because the FAERS ROR 12.79 contributes 0.30 weight. This correctly reflects
the convergent evidence from two independent methods.