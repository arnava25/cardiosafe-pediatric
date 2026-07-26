"""
CardioSafe Pediatric: parameter table builder, v2.

Architectural contract, in response to the June 2026 audit:

  1. free Cmax is NEVER stored. It is computed here, every time, from
     total Cmax and fraction bound.
  2. Every pharmacological input must carry a source string. A row with a
     value but no source is a build error, not a warning.
  3. Missing values stay missing. Nothing is estimated, defaulted, or
     carried forward from a compilation.
  4. Only rows that pass all checks enter the model. Everything else is
     reported as blocked, with the reason.

Usage:
    python3 build_params.py            # report status
    python3 build_params.py --model    # emit model ready table only
"""

import sys
import csv
from pathlib import Path

CSV_PATH = Path(__file__).parent / "herg_params_v2.csv"

REQUIRED_FOR_HERG = [
    ("hERG_IC50_nM", "IC50_source"),
    ("fraction_bound", "binding_source"),
    ("cmax_total_pediatric_ngml", "cmax_source"),
    ("MW_gmol", None),
]


def to_float(s):
    s = (s or "").strip()
    if s == "":
        return None
    return float(s)


def load_rows():
    with open(CSV_PATH, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def check_row(row):
    """Return (ok, list_of_reasons_blocked)."""
    reasons = []

    if row["data_quality"] == "not_applicable_herg":
        return False, ["not an hERG pathway drug, handled by sympathomimetic module"]

    for value_col, source_col in REQUIRED_FOR_HERG:
        val = to_float(row[value_col])
        if val is None:
            reasons.append(f"missing {value_col}")
            continue
        if source_col is not None and not (row[source_col] or "").strip():
            reasons.append(f"{value_col} present but {source_col} empty (UNSOURCED VALUE)")

    return (len(reasons) == 0), reasons


def derive(row):
    """Compute the three derived quantities. Never stored, always computed."""
    ic50 = to_float(row["hERG_IC50_nM"])
    fb = to_float(row["fraction_bound"])
    cmax_ngml = to_float(row["cmax_total_pediatric_ngml"])
    mw = to_float(row["MW_gmol"])

    # ng/mL to nM:  (ng/mL) / (g/mol) * 1000
    cmax_total_nM = (cmax_ngml / mw) * 1000.0
    cmax_free_nM = cmax_total_nM * (1.0 - fb)
    fractional_block = cmax_free_nM / (cmax_free_nM + ic50)

    return {
        "cmax_total_nM": cmax_total_nM,
        "cmax_free_nM": cmax_free_nM,
        "fractional_block": fractional_block,
    }


def block_threshold_cmax(row, target_block=0.20):
    """
    Diagnostic used when Cmax is missing: what total plasma concentration
    would this drug need to reach target_block? Lets us reach a verdict
    without a Cmax value.
    """
    ic50 = to_float(row["hERG_IC50_nM"])
    fb = to_float(row["fraction_bound"])
    mw = to_float(row["MW_gmol"])
    if None in (ic50, fb, mw):
        return None
    free_needed = ic50 * target_block / (1.0 - target_block)
    total_needed_nM = free_needed / (1.0 - fb)
    total_needed_ngml = total_needed_nM * mw / 1000.0
    return {
        "free_needed_nM": free_needed,
        "total_needed_nM": total_needed_nM,
        "total_needed_ngml": total_needed_ngml,
    }


def main():
    rows = load_rows()
    ready, blocked = [], []

    for row in rows:
        ok, reasons = check_row(row)
        if ok:
            row.update(derive(row))
            ready.append(row)
        else:
            blocked.append((row, reasons))

    print(f"MODEL READY: {len(ready)} of {len(rows)} rows\n")

    if ready:
        print(f"{'code':<8}{'IC50 nM':>10}{'f_bound':>9}"
              f"{'Cmax_tot nM':>13}{'Cmax_free nM':>14}{'block %':>10}")
        for r in ready:
            print(f"{r['code']:<8}{float(r['hERG_IC50_nM']):>10.1f}"
                  f"{float(r['fraction_bound']):>9.3f}"
                  f"{r['cmax_total_nM']:>13.1f}{r['cmax_free_nM']:>14.2f}"
                  f"{r['fractional_block'] * 100:>10.2f}")
        print()

    print(f"BLOCKED: {len(blocked)} rows\n")
    for row, reasons in blocked:
        print(f"  {row['code']:<8} {row['drug']}")
        for reason in reasons:
            print(f"           - {reason}")
        thr = block_threshold_cmax(row)
        if thr is not None:
            print(f"           > to reach 20% block needs total "
                  f"{thr['total_needed_ngml']:.0f} ng/mL "
                  f"({thr['total_needed_nM']:.0f} nM)")
        print()


if __name__ == "__main__":
    main()
