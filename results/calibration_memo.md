# Calibration memo

Generated 2026-07-26 by `src/calibrate.py`.

## Result

| Dose mg/day | Free Cmax nM | Block % | Model dAPD90 | FDA ddQTcF | Ratio |
|---|---|---|---|---|---|
| 400 | 230.9 | 2.71 | 4.5 | 5.4 | 0.840 |
| 600 | 453.8 | 5.18 | 8.5 | 6.8 | 1.255 |
| 800 | 493.6 | 5.61 | 9.5 | 6.9 | 1.382 |

Ratio spread 0.839 to 1.383, a factor of 1.65.

Candidate relationships, n = 3:

- multiplicative, ddQTcF = dAPD90 / 1.159, RMSE 1.20 ms
- additive, ddQTcF = dAPD90 minus 1.17, RMSE 1.49 ms
- linear, ddQTcF = 0.314 x dAPD90 + 4.00, RMSE 0.09 ms

## What this can and cannot support

With three anchor points spanning a narrow exposure range, this cannot
distinguish a multiplicative from an additive bias. It can show whether a
single correction factor is even roughly consistent, and the ratio spread
above is the answer to that.

Separately, and more fundamentally: single cell dAPD90 and surface ECG ddQTcF
are different measurements. QT reflects summed repolarization across a
heterogeneous ventricular wall plus conduction, not the duration of one
myocyte's action potential. Any factor here is an empirical correspondence over
a narrow range in one drug, not a conversion.

## Recommendation

If the ratio spread exceeds roughly 1.4, do not report adjusted values against
the ICH E14 10 ms threshold. Report dAPD90 directly, ranked, with the
mechanism, and give this comparison in the Discussion as the reason absolute
clinical thresholds are not claimed.

The risperidone anchor is a null across 175 pediatric patients and bounds the
model rather than calibrating it. A trial of that size would not reliably
detect a small change given typical QTc variance, so consistency with the null
is weak evidence.

## Observed versus predicted

The FDA numbers above are that agency's own exposure response MODEL. Observed
mean dQTcF in the two pivotal pediatric trials, roughly 500 patients, was about
2 ms. The model to observation gap is therefore
larger than the model to model gap, in both cases.
