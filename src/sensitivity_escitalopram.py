"""
Escitalopram sensitivity across the disputed clearance range.

Escitalopram is the only drug in the set whose APD90 effect approaches the ICH
E14 threshold of regulatory concern (10 ms for delta delta QTc). Its absolute
value depends on a pharmacokinetic parameter the literature does not agree on.

Published pediatric apparent clearance spans 14.2 to 40 L/h:

  14.2  Poweleit et al 2023, population PK, 20 mg/day normalized, CYP2C19 normal
        metabolizers. The measured value used in params/herg_params_v2.csv.
        The authors note this is BELOW the published range and attribute it to
        opportunistic sampling.
  16.3  A second low estimate Poweleit cite
  20-40 The range Poweleit report from six or seven prior adult and adolescent
        analyses
  33.0  Implied by the Fekete et al 2020 dose to concentration factor for minors

Concentration scales inversely with clearance, so cmax_scale = 14.2 / CL.

Also runs the CYP2C19 phenotype arm, which is the part of the finding that does
NOT depend on the absolute clearance. Poweleit Table 3 gives measured Cmax of
186.0 ng/mL in poor metabolizers against 73.73 in normal, a ratio of 2.52. That
ratio holds at any assumed population clearance.

See params/rebuild_record.md sections 5.1 and 6.1.

Usage:  cd src && python3 sensitivity_escitalopram.py
Output: results/sensitivity_escitalopram.csv
"""

import time
from pathlib import Path

import pandas as pd

import ord_model as M

BASELINE = 263.6
N_BEATS = 500
CL = 1000.0
ICH_THRESHOLD = 10.0      # ms, ICH E14 threshold of regulatory concern
CL_REFERENCE = 14.2       # L/h, the clearance the table's Cmax was measured at
# Measured Cmax ratio from Poweleit 2023 Table 3, normalized to 20 mg/day:
# poor metabolizer 186.0 ng/mL, normal metabolizer 73.73 ng/mL.
# Do NOT derive this from the 69 percent clearance difference. Cmax depends on
# absorption rate and volume of distribution as well as clearance, and the two
# disagree: 186.0/73.73 = 2.523 measured, against 1/(1-0.69) = 3.226 from
# clearance alone. An earlier version used 3.226 and overestimated the poor
# metabolizer arm by roughly 18 percent.
PM_RATIO = 186.0 / 73.73

CLEARANCES = [
    (14.2, "Poweleit 2023 measured, CYP2C19 normal metabolizers"),
    (16.3, "second low estimate cited by Poweleit"),
    (20.0, "bottom of the published adult and adolescent range"),
    (33.0, "implied by Fekete 2020 dose to concentration factor"),
    (40.0, "top of the published adult and adolescent range"),
]

OUT = Path(__file__).resolve().parent.parent / "results" / "sensitivity_escitalopram.csv"


def run(scale, label, phenotype, cl_assumed, rows, t_start):
    r = M.run_simulation({"ESC": "therapeutic"}, n_beats=N_BEATS, CL=CL,
                         verbose=False, cmax_scale=scale)
    apd = r["APD90"]
    d = apd - BASELINE
    blk = r["IKr_block_pct"]
    rows.append({
        "assumed_CL_L_per_h": cl_assumed,
        "phenotype": phenotype,
        "cmax_scale": round(scale, 4),
        "block_pct": round(blk, 3),
        "APD90_ms": round(apd, 2),
        "dAPD90_ms": round(d, 2),
        "above_ICH_10ms": bool(d > ICH_THRESHOLD),
        "source": label,
    })
    flag = "  ABOVE 10 ms" if d > ICH_THRESHOLD else ""
    print(f"{cl_assumed:>7.1f}{phenotype:>12}{scale:>11.3f}{blk:>9.2f}"
          f"{d:>+9.1f}{flag:<14}{time.time() - t_start:>7.0f}s", flush=True)


def main():
    n = len(CLEARANCES) * 2
    print("Escitalopram sensitivity across the disputed clearance range")
    print(f"Baseline {BASELINE} ms, {N_BEATS} beats, ICH E14 threshold "
          f"{ICH_THRESHOLD:.0f} ms")
    print(f"{n} runs, each {N_BEATS} beats. Slow.\n")
    print(f"{'CL':>7}{'phenotype':>12}{'cmax_sc':>11}{'block%':>9}"
          f"{'dAPD90':>9}{'':<14}{'elapsed':>7}")

    rows = []
    t_start = time.time()
    for cl, label in CLEARANCES:
        scale = CL_REFERENCE / cl
        run(scale, label, "normal", cl, rows, t_start)
        run(scale * PM_RATIO, label + ", CYP2C19 poor metabolizer",
            "poor", cl, rows, t_start)

    df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    print(f"\nWrote {OUT}")

    nm = df[df["phenotype"] == "normal"]
    pm = df[df["phenotype"] == "poor"]
    print(f"\nNormal metabolizers: dAPD90 spans {nm['dAPD90_ms'].min():+.1f} to "
          f"{nm['dAPD90_ms'].max():+.1f} ms across the clearance range")
    print(f"Poor metabolizers:   dAPD90 spans {pm['dAPD90_ms'].min():+.1f} to "
          f"{pm['dAPD90_ms'].max():+.1f} ms")
    print(f"\nAbove the 10 ms threshold: "
          f"{int(nm['above_ICH_10ms'].sum())} of {len(nm)} normal metabolizer "
          f"scenarios, {int(pm['above_ICH_10ms'].sum())} of {len(pm)} poor.")

    print("\nCALIBRATION CAVEAT. Quetiapine in this model gives +9.5 ms dAPD90. "
          "The FDA exposure response model for the same drug in the same "
          "population predicted 5.4 to 6.9 ms delta delta QTcF, and the two "
          "pivotal pediatric trials observed about 2 ms. This model appears to "
          "run roughly 1.4 to 1.8 times higher than the regulatory model. "
          "Single cell dAPD90 is also not surface ECG dQTc. Report these "
          "numbers with that comparison stated, not as QTc predictions.")


if __name__ == "__main__":
    main()
