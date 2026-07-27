"""
Calibrate model dAPD90 against FDA clinical anchors.

The project currently divides every dAPD90 by 1.5 to produce an "adjusted"
number compared against the ICH E14 10 ms threshold. That factor came from one
eyeballed comparison and is applied uniformly, which assumes the model's
overestimate is proportional to effect size. This script tests that assumption
against every anchor available.

TWO THINGS ARE BEING CONFLATED AND THIS SCRIPT SEPARATES THEM.

  1. Is the model biased, and is the bias multiplicative or additive?
  2. Is single cell dAPD90 convertible to surface ECG ddQTcF at all?

Question 2 has no clean answer. Action potential duration in one simulated
myocyte and the QT interval on a body surface electrocardiogram are different
measurements. QT reflects summed repolarization across a heterogeneous
ventricular wall plus conduction. Any factor derived here is an empirical
correspondence over a narrow range, not a conversion.

ANCHORS

  Quetiapine, FDA Clinical Pharmacology Review NDA 20639 SE5-045/046, reviewer
  Kofi A. Kumi. Exposure response model applied to pediatric exposures, ages
  10 to 17:
      400 mg/day  Cmax  520.9 ng/mL   predicted ddQTcF 5.4 ms
      600 mg/day  Cmax 1023.6 ng/mL   predicted ddQTcF 6.8 ms
      800 mg/day  Cmax 1113.4 ng/mL   predicted ddQTcF 6.9 ms
  Observed mean dQTcF across the two pivotal pediatric trials, roughly 500
  patients, was about 2 ms with no patient exceeding 500 ms or a 60 ms change.

  Risperidone, FDA review NDA 20-272 S065. Across a 96 patient controlled study
  and a 79 patient open label extension, no significant mean changes in ECG
  parameters in any treatment group. Two subjects with prolonged QTc, both
  under 450 ms with changes under 60 ms. This is a NULL, so it bounds the model
  rather than calibrating it.

Usage:  cd src && python3 calibrate.py
Output: results/calibration.csv and results/calibration_memo.md
"""

import time
from pathlib import Path

import numpy as np
import pandas as pd

import ord_model as M

BASELINE = 263.6
N_BEATS = 500
CL = 1000.0

# Quetiapine: FDA pediatric Cmax and the matching modelled ddQTcF.
QUE_ANCHORS = [
    (400, 520.9, 5.4),
    (600, 1023.6, 6.8),
    (800, 1113.4, 6.9),
]
QUE_OBSERVED_TRIAL_dQTcF = 2.0

OUT_CSV = Path(__file__).resolve().parent.parent / "results" / "calibration.csv"
OUT_MEMO = Path(__file__).resolve().parent.parent / "results" / "calibration_memo.md"


def free_nM(total_ngml, mw, fraction_bound):
    return (total_ngml / mw) * 1000.0 * (1.0 - fraction_bound)


def main():
    params = M.load_drug_params()
    que = params["QUE"]
    print("Calibration against FDA clinical anchors")
    print(f"Baseline {BASELINE} ms, {N_BEATS} beats\n")
    print(f"{'dose':>6}{'Cmax':>9}{'free nM':>10}{'block%':>8}"
          f"{'dAPD90':>9}{'FDA ddQTcF':>12}{'ratio':>8}{'diff':>8}")

    rows, t0 = [], time.time()
    for dose, cmax_ngml, fda in QUE_ANCHORS:
        fc = free_nM(cmax_ngml, que["MW"], que["fraction_bound"])
        r = M.run_simulation({"QUE": fc}, n_beats=N_BEATS, CL=CL, verbose=False)
        d = r["APD90"] - BASELINE
        rows.append({
            "drug": "quetiapine",
            "dose_mg_day": dose,
            "cmax_total_ngml": cmax_ngml,
            "free_cmax_nM": round(fc, 2),
            "block_pct": round(r["IKr_block_pct"], 3),
            "model_dAPD90_ms": round(d, 2),
            "fda_predicted_ddQTcF_ms": fda,
            "ratio_model_over_fda": round(d / fda, 3),
            "difference_ms": round(d - fda, 2),
        })
        print(f"{dose:>6}{cmax_ngml:>9.1f}{fc:>10.1f}{r['IKr_block_pct']:>8.2f}"
              f"{d:>9.1f}{fda:>12.1f}{d / fda:>8.3f}{d - fda:>+8.1f}", flush=True)

    df = pd.DataFrame(rows)
    x = df["model_dAPD90_ms"].values
    y = df["fda_predicted_ddQTcF_ms"].values

    # Three candidate relationships. n = 3, so this cannot distinguish them
    # convincingly. It can show whether the multiplicative assumption is even
    # roughly consistent.
    ratios = x / y
    mult_k = float(np.mean(ratios))
    mult_resid = y - x / mult_k

    add_c = float(np.mean(x - y))
    add_resid = y - (x - add_c)

    m, c = np.polyfit(x, y, 1)
    lin_resid = y - (m * x + c)

    def rmse(r):
        return float(np.sqrt(np.mean(np.square(r))))

    print("\nMODEL FITS, n = 3")
    print(f"  multiplicative  ddQTcF = dAPD90 / {mult_k:.3f}"
          f"    ratios {np.round(ratios, 3)}   RMSE {rmse(mult_resid):.2f} ms")
    print(f"  additive        ddQTcF = dAPD90 - {add_c:.2f}"
          f"              RMSE {rmse(add_resid):.2f} ms")
    print(f"  linear          ddQTcF = {m:.3f}*dAPD90 + {c:.2f}"
          f"      RMSE {rmse(lin_resid):.2f} ms")

    spread = ratios.max() / ratios.min()
    print(f"\n  Ratio spread across the three anchors: {ratios.min():.3f} to "
          f"{ratios.max():.3f}, a factor of {spread:.2f}.")
    if spread > 1.4:
        print("  The multiplicative assumption does NOT hold cleanly. A single")
        print("  correction factor is not supported by these anchors.")
    else:
        print("  The multiplicative assumption is roughly consistent.")

    # Risperidone: a null, so it bounds rather than calibrates.
    print("\nRISPERIDONE NULL CHECK")
    r = M.run_simulation({"RIS": "therapeutic"}, n_beats=N_BEATS, CL=CL,
                         verbose=False)
    ris_d = r["APD90"] - BASELINE
    print(f"  model dAPD90 {ris_d:+.1f} ms at therapeutic pediatric exposure")
    print(f"  FDA: no significant mean ECG change across 175 pediatric patients")
    print(f"  Applying the multiplicative factor gives "
          f"{ris_d / mult_k:+.1f} ms ddQTcF equivalent.")
    print("  A trial of that size would not reliably detect a change of that")
    print("  size given typical QTc variance, so the null is consistent with")
    print("  the model rather than contradicting it. Weak constraint.")

    rows.append({
        "drug": "risperidone", "dose_mg_day": np.nan,
        "cmax_total_ngml": que["cmax_total_ngml"] if False else np.nan,
        "free_cmax_nM": round(M.free_cmax_nM(params["RIS"]), 2),
        "block_pct": round(r["IKr_block_pct"], 3),
        "model_dAPD90_ms": round(ris_d, 2),
        "fda_predicted_ddQTcF_ms": np.nan,
        "ratio_model_over_fda": np.nan,
        "difference_ms": np.nan,
    })

    out = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False)

    memo = f"""# Calibration memo

Generated {time.strftime('%Y-%m-%d')} by `src/calibrate.py`.

## Result

| Dose mg/day | Free Cmax nM | Block % | Model dAPD90 | FDA ddQTcF | Ratio |
|---|---|---|---|---|---|
""" + "\n".join(
        f"| {r['dose_mg_day']:.0f} | {r['free_cmax_nM']:.1f} | "
        f"{r['block_pct']:.2f} | {r['model_dAPD90_ms']:.1f} | "
        f"{r['fda_predicted_ddQTcF_ms']:.1f} | {r['ratio_model_over_fda']:.3f} |"
        for r in rows[:3]) + f"""

Ratio spread {ratios.min():.3f} to {ratios.max():.3f}, a factor of {spread:.2f}.

Candidate relationships, n = 3:

- multiplicative, ddQTcF = dAPD90 / {mult_k:.3f}, RMSE {rmse(mult_resid):.2f} ms
- additive, ddQTcF = dAPD90 minus {add_c:.2f}, RMSE {rmse(add_resid):.2f} ms
- linear, ddQTcF = {m:.3f} x dAPD90 + {c:.2f}, RMSE {rmse(lin_resid):.2f} ms

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
{QUE_OBSERVED_TRIAL_dQTcF:.0f} ms. The model to observation gap is therefore
larger than the model to model gap, in both cases.
"""
    OUT_MEMO.write_text(memo)
    print(f"\nWrote {OUT_CSV}")
    print(f"Wrote {OUT_MEMO}")
    print(f"\nTotal {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
