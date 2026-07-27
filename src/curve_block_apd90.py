"""
Characterize the engine's hERG block to APD90 response curve.

The mapping from fractional IKr block to APD90 prolongation is a property of
ord_core, not of any particular drug. Measure it once across a wide block range
and every drug projection becomes arithmetic.

The six-drug run on 26 July 2026 found this relationship linear at 1.79 ms per
percent block (R^2 = 0.9996) across 1.5 to 12.5 percent. Repolarization reserve
predicts it should bend upward once reserve is consumed. This finds where.

Uses escitalopram as the carrier drug purely because it has no metabolite, so
requesting an explicit concentration applies exactly that concentration. The
drug identity is irrelevant; only the resulting GKr multiplier reaches the cell.

NOTE ON PRECISION: run_simulation samples the final beat at 2000 points over
CL ms, so APD90 resolves to 0.5 ms at CL = 1000. Fine at these magnitudes.

Usage:  cd src && python3 curve_block_apd90.py
Output: results/curve_block_apd90.csv
"""

import time
from pathlib import Path

import numpy as np
import pandas as pd

import ord_model as M

BASELINE = 263.6          # frozen validated drug free APD90 at 500 beats, 60 bpm
IC50_ESC = 700.0          # nM, escitalopram, used only to invert to a concentration
N_BEATS = 500
CL = 1000.0

# Block targets. The first five bracket the measured drug range; the rest push
# well past it to find the bend.
TARGETS = [0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70]

OUT = Path(__file__).resolve().parent.parent / "results" / "curve_block_apd90.csv"


def free_conc_for_block(target, ic50):
    """Invert the Hill equation: what free concentration gives this block?"""
    return ic50 * target / (1.0 - target)


def main():
    print(f"Engine block to APD90 curve, {N_BEATS} beats at CL={CL:.0f} ms")
    print(f"Baseline reference: {BASELINE} ms")
    print(f"{len(TARGETS)} runs. This is slow, each run is {N_BEATS} beats.\n")
    print(f"{'target':>8}{'conc nM':>10}{'block %':>9}{'APD90':>9}"
          f"{'dAPD90':>9}{'ms per %':>10}{'elapsed':>9}")

    rows = []
    t_start = time.time()
    for target in TARGETS:
        conc = free_conc_for_block(target, IC50_ESC)
        r = M.run_simulation({"ESC": conc}, n_beats=N_BEATS, CL=CL, verbose=False)
        apd = r["APD90"]
        d = apd - BASELINE
        blk = r["IKr_block_pct"]
        rows.append({
            "target_block": target,
            "free_conc_nM": round(conc, 2),
            "block_pct": round(blk, 3),
            "APD90_ms": round(apd, 2),
            "dAPD90_ms": round(d, 2),
            "ms_per_pct_block": round(d / blk, 4) if blk else np.nan,
        })
        print(f"{target:>8.2f}{conc:>10.1f}{blk:>9.2f}{apd:>9.1f}"
              f"{d:>+9.1f}{d / blk:>10.3f}{time.time() - t_start:>8.0f}s",
              flush=True)

    df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    print(f"\nWrote {OUT}")

    # Where does linearity break? Fit the low-block region and see where the
    # measured curve departs from that fit by more than the 0.5 ms grid.
    low = df[df["block_pct"] <= 20.0]
    if len(low) >= 3:
        m, c = np.polyfit(low["block_pct"], low["dAPD90_ms"], 1)
        pred = m * df["block_pct"] + c
        resid = df["dAPD90_ms"] - pred
        print(f"\nFit over block <= 20%: dAPD90 = {m:.3f} * block% + {c:.3f}")
        print(f"{'block %':>9}{'measured':>10}{'linear pred':>13}{'residual':>10}")
        for b, meas, p, res in zip(df["block_pct"], df["dAPD90_ms"], pred, resid):
            flag = "  <-- departs" if abs(res) > 1.0 else ""
            print(f"{b:>9.2f}{meas:>10.1f}{p:>13.1f}{res:>+10.1f}{flag}")
        print("\nIf residuals stay near zero throughout, the response is linear "
              "across the whole range and reserve is not being consumed even at "
              "70 percent block. If they turn positive, that is the bend, and "
              "the block value where it starts is the number to report.")


if __name__ == "__main__":
    main()
