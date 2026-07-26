#!/usr/bin/env python3
"""
generate_sim_data.py - CardioSafe Pediatric (curiosity build)
Regenerates the model-derived data constants in docs/clinical_sim.html from the
validated results/risk_grid_results.csv:
  DQTC    - Bazett delta-QTc per pair
  DAPD    - genuine delta-APD90 per pair  (the real repolarization metric)
  IKR     - hERG block percent per pair
  TRIPLES - triples with genuine + Bazett + tier-on-genuine
Also removes the obsolete BAZETT decomposition block.
FAERS_ROR / FAERS_SIGNAL / AGE_ROR are model-independent and left untouched.

Usage:  python3 src/generate_sim_data.py [--dry-run]

DO NOT RUN until the risk grid is regenerated (July 2026). This script writes
its output directly into docs/clinical_sim.html, which is the live page at
arnava25.github.io and is currently carrying an "under revision" banner because
the parameter set behind it was withdrawn in the June 2026 audit. Running it
against the existing results/risk_grid_results.csv would refresh the live
calculator with the same invalidated numbers. See results/README.md and
params/rebuild_record.md.
"""
import sys, re, csv, json, argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GRID = ROOT / "results" / "risk_grid_results.csv"
HTML = ROOT / "docs" / "clinical_sim.html"

def key(c): return "+".join(sorted(str(c).split("+")))
def tier(g): return "HIGH" if g>=20 else "MODERATE" if g>=10 else "LOW-MOD" if g>=5 else "LOW" if g>=0 else "PROTECTIVE"

def load(path):
    dq, da, ik, trip = {}, {}, {}, []
    with open(path) as f:
        for r in csv.DictReader(f):
            combo = r.get("combination","").strip()
            if not combo: continue
            parts = combo.split("+")
            def col(*names):
                for n in names:
                    if n in r and r[n] not in ("", None): return r[n]
                return None
            try:
                bz = float(col("ΔQTc_ms","dQTc","delta_qtc"))
                ap = float(col("ΔAPD90_ms","dAPD90","delta_apd"))
                ir = float(col("IKr_block_pct","ikr_block_pct"))
            except (TypeError, ValueError):
                continue
            if len(parts) == 2:
                k = key(combo); dq[k]=round(bz,1); da[k]=round(ap,1); ik[k]=round(ir,2)
            elif len(parts) == 3:
                trip.append({"combo":combo,"genuine":round(ap,1),"bazett":round(bz,1),"tier":tier(ap)})
    return dq, da, ik, trip

def obj(d): return "{\n" + ",\n".join(f'  "{k}":{json.dumps(v)}' for k,v in sorted(d.items())) + "\n}"

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--dry-run", action="store_true"); a = ap.parse_args()
    dq, da, ik, trip = load(GRID)
    print(f"pairs: {len(dq)}  triples: {len(trip)}")
    blocks = {
        "DQTC": "const DQTC=" + obj(dq) + ";",
        "DAPD": "const DAPD=" + obj(da) + ";",
        "IKR":  "const IKR="  + obj(ik) + ";",
        "TRIPLES": "const TRIPLES=" + json.dumps(trip) + ";",
    }
    if a.dry_run:
        print("\n\n".join(blocks.values())); return
    html = HTML.read_text()
    html = re.sub(r'\n*//[^\n]*Bazett decomposition \(mean[^\n]*', '', html)
    html = re.sub(r'\nconst BAZETT\s*=\s*\{[\s\S]*?\n\};', '', html)
    html = re.sub(r'\nconst DAPD\s*=\s*\{[\s\S]*?\n\};', '', html)   # idempotent
    html = re.sub(r'const DQTC\s*=\s*\{[\s\S]*?\n\};', lambda m: blocks["DQTC"]+"\n\n"+blocks["DAPD"], html, count=1)
    html = re.sub(r'const IKR\s*=\s*\{[\s\S]*?\n\};',  lambda m: blocks["IKR"], html, count=1)
    html = re.sub(r'const TRIPLES\s*=\s*\[[\s\S]*?\n\];', lambda m: blocks["TRIPLES"], html, count=1)
    HTML.write_text(html)
    print("patched", HTML)

if __name__ == "__main__":
    main()