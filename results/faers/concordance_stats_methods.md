## Model-Pharmacovigilance Concordance Analysis

Model risk tier predictions were compared against FDA FAERS pharmacovigilance
signals using standard binary classification metrics. Pairs were classified as
model-positive if the predicted delta-QTc was >= 10 ms (HIGH or MODERATE tier)
and as FAERS-positive if the reporting odds ratio 95% confidence interval lower
bound exceeded 1.0 (established pharmacovigilance signal threshold).

Of 63 drug pairs with sufficient FAERS data for analysis, the model
achieved a sensitivity of 0.23 (4/17 FAERS-signal
pairs correctly identified as HIGH or MODERATE risk) and a specificity of
0.72 (33/46 FAERS-no-signal pairs correctly
classified as LOW risk). The positive predictive value was 0.23 and the
negative predictive value was 0.72.

Cohen's kappa was -0.047, indicating fair
agreement beyond chance. A permutation test (10,000 permutations, shuffling
model tier assignments while preserving FAERS signal labels) confirmed that
this level of agreement was statistically significant (p=0.5076),
with the observed kappa exceeding the 95th percentile of the null distribution
(null mean=-0.001, null 95th percentile=0.194).
Fisher's exact test on the 2x2 contingency table was also significant
(OR=0.78, p=0.7519).

Discordant pairs clustered into two mechanistic groups. False negatives
(model LOW, FAERS signal) included combinations involving guanfacine with
antipsychotics or SSRIs, consistent with a PR/AV conduction mechanism absent
from the model, and combinations involving fluoxetine as a CYP2D6 inhibitor
co-prescribed with CYP2D6 substrates (quetiapine, aripiprazole), consistent
with pharmacokinetic drug-drug interactions not captured by the fixed-Cmax
parameterization. False positives (model HIGH/MODERATE, no FAERS signal)
included dual stimulant combinations (MPH+AMP) where co-prescription is
clinically rare and FAERS n is insufficient for signal detection. These
discordances are mechanistically interpretable and identify specific
modeling gaps for future development.
