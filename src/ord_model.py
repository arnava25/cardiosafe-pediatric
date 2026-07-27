"""
CardioSafe Pediatric: drug module on the validated canonical O'Hara-Rudy core.

Rewritten July 2026 against the rebuilt parameter table (params/herg_params_v2.csv).
Same run_simulation interface as before, so risk_grid.py and rate_correction.py
do not change. Four things are different underneath:

  1. FREE CMAX IS COMPUTED, NEVER READ.
     free = (total_ngml / MW) * 1000 * (1 - fraction_bound)
     The old module read a cmax_free_nM column. That column was hand typed and
     the protein binding fraction was never applied to most of it. It no longer
     exists and nothing here will accept one.

  2. HARD FAILURE ON INCOMPLETE DATA.
     The old load_drug_params caught FileNotFoundError, printed a warning, and
     returned {}. Every drug then silently missed the lookup, total block came
     out 0.0, and the simulation produced clean baseline numbers for a drug
     combination. That is how a broken parameter path yields plausible wrong
     results. This module raises instead, at every point where the old one
     shrugged.

  3. METABOLITES ARE PAIRED TO THEIR PARENT, NOT PRESCRIBED SEPARATELY.
     Nobody prescribes paliperidone alongside risperidone; it arrives
     obligately. Requesting RIS now applies RIS and 9OHRIS block together, each
     at its own free concentration. Same for FLU and NORFLU. Metabolites cannot
     be requested directly.

  4. IC50 AND CONCENTRATION SCALING HOOKS.
     ic50_scale and cmax_scale multiply every drug's value, so the sweep can be
     run as a sensitivity analysis rather than at point estimates. Published
     hERG IC50s vary 1.9-fold within one lab and 4.4-fold across labs, and
     escitalopram clearance estimates span 14.2 to 40 L/h. Point estimates
     misrepresent both. See params/rebuild_record.md sections 4.3 and 5.1.

Run from src/ or with src/ on the path, as before.
"""

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import odeint

warnings.filterwarnings("ignore")

import ord_core as C  # the validated engine

# ---------------------------------------------------------------------------
# Parameter table
# ---------------------------------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent
PARAMS_PATH = ROOT / "params" / "herg_params_v2.csv"

# Mechanism is derived from class rather than stored, so the table stays a
# record of measurements and the model owns the modelling decisions.
_MECHANISM_BY_CLASS = {
    "Stimulant": "sympathomimetic",
    "Alpha-2": "autonomic_modulation",
    "Antipsychotic": "hERG_block",
    "Antipsychotic_metabolite": "hERG_block",
    "SSRI": "hERG_block",
    "SSRI_metabolite": "hERG_block",
    "TCA": "hERG_block",
}

# Rows excluded from hERG block by documented argument rather than by
# measurement. See rebuild_record.md section 4.5. Requesting one is allowed;
# it contributes zero block and says so.
_EXCLUDED_BY_ARGUMENT = {
    "ARI": "no hERG IC50 exists; would need IC50 below 46 nM to reach 10% block",
    "NORFLU": "total concentration 490 nM is below the 625 nM free requirement "
              "for 20% block, unreachable at zero protein binding",
}


def _f(x):
    """
    Blank-safe float. Returns None for missing, never 0.0 and never NaN.

    Pandas reads an empty CSV cell as float NaN, not as an empty string. An
    earlier version of this function returned that NaN through, which made
    model_ready() pass for rows with no IC50 and no binding fraction. NaN then
    propagated silently into the block arithmetic. Missing must be None so the
    gate can see it.
    """
    if x is None:
        return None
    try:
        if pd.isna(x):
            return None
    except (TypeError, ValueError):
        pass
    s = str(x).strip()
    if s == "" or s.lower() == "nan":
        return None
    return float(s)


def _s(x):
    """Blank-safe string, for annotation columns."""
    if x is None:
        return ""
    try:
        if pd.isna(x):
            return ""
    except (TypeError, ValueError):
        pass
    return str(x).strip()


def load_drug_params(csv_path=None):
    """
    Load the parameter table keyed by drug code (RIS, QUE, ESC...).

    Raises FileNotFoundError rather than warning. A missing parameter table
    must stop the run.
    """
    path = Path(csv_path) if csv_path else PARAMS_PATH
    if not path.exists():
        raise FileNotFoundError(
            f"Parameter table not found at {path}. This module requires the "
            f"rebuilt table. The old data/herg_master_params.csv was withdrawn "
            f"in June 2026 and is in archive/params_invalidated_202606/. "
            f"See params/rebuild_record.md."
        )

    df = pd.read_csv(path)
    required = {"code", "drug", "class", "is_metabolite", "parent_code",
                "hERG_IC50_nM", "fraction_bound",
                "cmax_total_pediatric_ngml", "MW_gmol", "data_quality"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"{path} is missing required columns: {sorted(missing)}. "
            f"If it has a cmax_free_nM column it is the withdrawn table."
        )
    if "cmax_free_nM" in df.columns:
        raise ValueError(
            f"{path} contains a cmax_free_nM column. Free Cmax must be computed, "
            f"not stored. This is the error that invalidated the June 2026 table."
        )

    params = {}
    for _, r in df.iterrows():
        cls = r["class"]
        if cls not in _MECHANISM_BY_CLASS:
            raise ValueError(f"Unknown class '{cls}' for {r['code']}. "
                             f"Add it to _MECHANISM_BY_CLASS explicitly.")
        params[r["code"]] = {
            "drug": r["drug"],
            "class": cls,
            "mechanism": _MECHANISM_BY_CLASS[cls],
            "is_metabolite": int(r["is_metabolite"]) == 1,
            "parent_code": _s(r["parent_code"]) or None,
            "ic50_nM": _f(r["hERG_IC50_nM"]),
            "fraction_bound": _f(r["fraction_bound"]),
            "cmax_total_ngml": _f(r["cmax_total_pediatric_ngml"]),
            "MW": _f(r["MW_gmol"]),
            "assay_temp_C": _s(r.get("assay_temp_C", "")),
            "data_quality": r["data_quality"],
        }

    # Index metabolites under their parent so pairing is a lookup, not a rule.
    for code, dp in params.items():
        dp["metabolites"] = [c for c, m in params.items()
                             if m["is_metabolite"] and m["parent_code"] == code]
    return params


def free_cmax_nM(dp):
    """
    Compute free Cmax from total concentration, molecular weight and bound
    fraction. Never stored. Returns None if any input is missing.
    """
    if dp["cmax_total_ngml"] is None or dp["MW"] is None or dp["fraction_bound"] is None:
        return None
    total_nM = (dp["cmax_total_ngml"] / dp["MW"]) * 1000.0
    return total_nM * (1.0 - dp["fraction_bound"])


def model_ready(dp):
    """
    True if this row can contribute a computed hERG block.

    Belt and braces: also rejects NaN explicitly, so a future change to _f
    cannot silently reopen the hole where a row with no IC50 passed the gate.
    """
    if dp["mechanism"] != "hERG_block":
        return False
    ic50, fc = dp["ic50_nM"], free_cmax_nM(dp)
    if ic50 is None or fc is None:
        return False
    if not (np.isfinite(ic50) and np.isfinite(fc)):
        return False
    return True


# ---------------------------------------------------------------------------
# Drug combination expansion
# ---------------------------------------------------------------------------

def expand_combination(drug_combo, drug_params):
    """
    Expand a requested combination to include obligate active metabolites.

    Input keys are drug codes. Values are "therapeutic" or an explicit free
    concentration in nM. Requesting RIS yields RIS plus 9OHRIS; requesting FLU
    yields FLU plus NORFLU. Metabolites cannot be requested directly, because
    they are not prescribed.
    """
    expanded = {}
    for code, val in drug_combo.items():
        if code not in drug_params:
            raise KeyError(
                f"'{code}' is not in the parameter table. Valid codes: "
                f"{sorted(drug_params)}"
            )
        dp = drug_params[code]
        if dp["is_metabolite"]:
            raise ValueError(
                f"'{code}' is a metabolite of {dp['parent_code']} and cannot be "
                f"requested directly. Request {dp['parent_code']}; the metabolite "
                f"is applied automatically."
            )
        expanded[code] = val
        for met in dp["metabolites"]:
            # Metabolite always at its own measured therapeutic concentration.
            # An explicit parent concentration does not scale the metabolite,
            # because the parent-to-metabolite ratio is CYP-dependent and not
            # modelled here.
            expanded[met] = "therapeutic"
    return expanded


def compute_ikr_block(drug_combo, drug_params, ic50_scale=1.0, cmax_scale=1.0,
                      strict=True):
    """
    Total fractional hERG block, combined as independent binding sites.

    WARNING, UNTESTED ASSUMPTION. Independent-site combination assumes the drugs
    do not compete for the same pocket. Nearly all hERG blockers bind the
    aromatic residues Y652 and F656 in the S6 pore cavity, and Rajamani 2006
    showed fluoxetine block is abolished by F656 mutation. If two drugs in a
    combination share that site, this overestimates combined block. This
    assumption underlies every polypharmacy number in the project and has never
    been checked. See rebuild_record.md section 8.

    ic50_scale and cmax_scale multiply every drug for sensitivity analysis.
    strict=True raises on a hERG-class drug that has no usable parameters.
    """
    combo = expand_combination(drug_combo, drug_params)
    blocks, notes = {}, {}

    for code, val in combo.items():
        dp = drug_params[code]

        if dp["mechanism"] != "hERG_block":
            blocks[code] = 0.0
            continue

        if code in _EXCLUDED_BY_ARGUMENT:
            blocks[code] = 0.0
            notes[code] = f"excluded by argument: {_EXCLUDED_BY_ARGUMENT[code]}"
            continue

        if not model_ready(dp):
            msg = (f"'{code}' ({dp['drug']}) is a hERG-class drug but is not "
                   f"model ready: ic50={dp['ic50_nM']}, "
                   f"free_cmax={free_cmax_nM(dp)}, "
                   f"data_quality={dp['data_quality']}")
            if strict:
                raise ValueError(
                    msg + ". Pass strict=False to treat it as zero block, but "
                    "understand that produces a silently incomplete result."
                )
            blocks[code] = 0.0
            notes[code] = "NOT MODEL READY, treated as zero block"
            continue

        ic50 = dp["ic50_nM"] * ic50_scale
        if val == "therapeutic":
            conc = free_cmax_nM(dp) * cmax_scale
        else:
            conc = float(val)
        blocks[code] = conc / (conc + ic50)

    total = 1.0 - np.prod([1.0 - b for b in blocks.values()]) if blocks else 0.0
    return total, blocks, notes


IKS_UPREG = 1.20  # per-sympathomimetic IKs upregulation; set 1.0 to disable


def compute_autonomic_modifiers(drug_combo, drug_params):
    """
    Rate and conductance modifiers for the sympathomimetic and alpha-2 pathways.

    UNCHANGED AND STILL BROKEN. These effects are hardcoded and not
    concentration dependent, which was flagged in the June 2026 audit and makes
    any stimulant result a mathematical demonstration rather than a
    pharmacological finding. Do not report stimulant-containing combinations
    until this is rebuilt against published concentration-response data.
    See rebuild_record.md section 8.
    """
    hr_mult = gcal_mult = gna_mult = gks_mult = 1.0
    for code in expand_combination(drug_combo, drug_params):
        mech = drug_params[code]["mechanism"]
        if mech == "sympathomimetic":
            hr_mult *= 0.90
            gcal_mult *= 1.15
            gks_mult *= IKS_UPREG
        elif mech == "autonomic_modulation":
            hr_mult *= 1.10
            gna_mult *= 0.95
    return hr_mult, gcal_mult, gna_mult, gks_mult


def contains_unrebuilt_pathway(drug_combo, drug_params):
    """True if any member uses the hardcoded, non-concentration-dependent path."""
    return any(
        drug_params[c]["mechanism"] in ("sympathomimetic", "autonomic_modulation")
        for c in expand_combination(drug_combo, drug_params)
    )


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------

def _apd90_qtc(V, t, CL_eff):
    Vmax = np.max(V)
    Vrest = np.min(V[len(V) // 2:])
    V90 = Vrest + 0.10 * (Vmax - Vrest)
    iu = int(np.argmax(V))
    post, tp = V[iu:], t[iu:] - t[iu]
    cr = np.where(post <= V90)[0]
    apd = tp[cr[0]] if len(cr) else np.nan
    qtc = apd / np.sqrt(CL_eff / 1000) if not np.isnan(apd) else np.nan
    return apd, qtc, Vmax


def run_simulation(drug_combination=None, drug_params_path=None, n_beats=500,
                   CL=1000.0, stim_amp=-80.0, stim_dur=0.5, verbose=True,
                   ic50_scale=1.0, cmax_scale=1.0, strict=True):
    """
    Run to steady state and return APD90, Bazett QTc and diagnostics.

    n_beats defaults to 500, not 200. The frozen validated baseline of
    263.6 ms is a 500-beat value; 200 beats does not reach steady state.
    """
    drug_params = load_drug_params(drug_params_path)

    if drug_combination:
        total_ikr_block, ikr_breakdown, block_notes = compute_ikr_block(
            drug_combination, drug_params,
            ic50_scale=ic50_scale, cmax_scale=cmax_scale, strict=strict)
        GKr_mult = 1.0 - total_ikr_block
        hr_mult, gcal_mult, gna_mult, gks_mult = compute_autonomic_modifiers(
            drug_combination, drug_params)
        CL_eff = CL * hr_mult
        unrebuilt = contains_unrebuilt_pathway(drug_combination, drug_params)
        expanded = list(expand_combination(drug_combination, drug_params))
    else:
        GKr_mult = gcal_mult = gna_mult = gks_mult = 1.0
        CL_eff = CL
        total_ikr_block, ikr_breakdown, block_notes = 0.0, {}, {}
        unrebuilt, expanded = False, []

    if verbose:
        label = " + ".join(expanded) if expanded else "Baseline"
        print(f"Simulation: {label}")
        print(f"  GKr={GKr_mult:.4f} (IKr block {total_ikr_block * 100:.2f}%)  "
              f"CL_eff={CL_eff:.1f}  beats={n_beats}")
        if ic50_scale != 1.0 or cmax_scale != 1.0:
            print(f"  sensitivity: ic50_scale={ic50_scale}  cmax_scale={cmax_scale}")
        for c, n in block_notes.items():
            print(f"  {c}: {n}")
        if unrebuilt:
            print("  WARNING: uses the hardcoded sympathomimetic or alpha-2 "
                  "pathway. Not concentration dependent. Do not report.")

    y = C.Y0.copy()
    last_t = last_y = None
    for beat in range(n_beats):
        t0 = beat * CL_eff
        ts = np.linspace(t0, t0 + CL_eff, 2000 if beat == n_beats - 1 else 300)
        sol = odeint(C.rhs, y, ts,
                     args=(GKr_mult, gcal_mult, gna_mult, gks_mult,
                           stim_amp, stim_dur, CL_eff),
                     rtol=1e-6, atol=1e-8, mxstep=8000)
        if not np.all(np.isfinite(sol[-1])):
            raise RuntimeError(
                f"Integration produced non-finite state at beat {beat} for "
                f"{expanded or 'baseline'}. The old module printed a warning and "
                f"returned NaN, which propagates silently into results tables."
            )
        y = sol[-1]
        if beat == n_beats - 1:
            last_t, last_y = ts - ts[0], sol

    V = last_y[:, 0]
    apd, qtc, vpeak = _apd90_qtc(V, last_t, CL_eff)
    if verbose:
        print(f"  APD90: {apd:.1f} ms  |  Bazett QTc: {qtc:.1f} ms  |  "
              f"Vpeak: {vpeak:.1f} mV")

    return {"APD90": apd, "QTc": qtc, "Vpeak": vpeak,
            "IKr_block_pct": total_ikr_block * 100, "GKr_scale": GKr_mult,
            "CL_effective": CL_eff, "drug_breakdown": ikr_breakdown,
            "block_notes": block_notes, "expanded_combination": expanded,
            "uses_unrebuilt_pathway": unrebuilt,
            "ic50_scale": ic50_scale, "cmax_scale": cmax_scale,
            "n_beats": n_beats, "t_trace": last_t, "V_trace": V}


def run_polypharmacy_sweep(combinations, drug_params_path=None, n_beats=500,
                           CL=1000.0, ic50_scale=1.0, cmax_scale=1.0,
                           strict=True):
    """
    Sweep a list of combinations. Returns a DataFrame plus the baseline result.

    Rows flag whether they used the unrebuilt stimulant or alpha-2 pathway, so
    a downstream table cannot mix defensible and undefensible results without
    showing it.
    """
    base = run_simulation(None, drug_params_path, n_beats, CL, verbose=False)
    bAPD, bQTc = base["APD90"], base["QTc"]
    print(f"Baseline APD90={bAPD:.1f} QTc={bQTc:.1f} at {n_beats} beats")
    if ic50_scale != 1.0 or cmax_scale != 1.0:
        print(f"Sensitivity: ic50_scale={ic50_scale} cmax_scale={cmax_scale}")
    print()

    rows = []
    for combo in combinations:
        r = run_simulation(combo, drug_params_path, n_beats, CL, verbose=False,
                           ic50_scale=ic50_scale, cmax_scale=cmax_scale,
                           strict=strict)
        rows.append({
            "combination": " + ".join(combo.keys()),
            "expanded": " + ".join(r["expanded_combination"]),
            "n_drugs_requested": len(combo),
            "APD90_ms": round(r["APD90"], 1),
            "dAPD90_ms": round(r["APD90"] - bAPD, 1),
            "QTc_ms": round(r["QTc"], 1),
            "dQTc_ms": round(r["QTc"] - bQTc, 1),
            "IKr_block_pct": round(r["IKr_block_pct"], 2),
            "uses_unrebuilt_pathway": r["uses_unrebuilt_pathway"],
            "ic50_scale": ic50_scale,
            "cmax_scale": cmax_scale,
        })
        flag = "  [UNREBUILT PATHWAY]" if r["uses_unrebuilt_pathway"] else ""
        print(f"  {rows[-1]['combination'][:40]:40s} "
              f"dAPD90={rows[-1]['dAPD90_ms']:+7.1f}  "
              f"dQTc={rows[-1]['dQTc_ms']:+7.1f}{flag}")

    return pd.DataFrame(rows).sort_values("dAPD90_ms", ascending=False), base


def model_ready_codes(drug_params=None):
    """Codes that can contribute a computed hERG block. The defensible set."""
    dp = drug_params or load_drug_params()
    return sorted(c for c, d in dp.items() if model_ready(d) and not d["is_metabolite"])


if __name__ == "__main__":
    p = load_drug_params()
    print(f"Loaded {len(p)} rows from {PARAMS_PATH}\n")
    print(f"{'code':<8}{'mechanism':<22}{'IC50':>8}{'free nM':>10}{'block':>9}  temp")
    for code in sorted(p):
        d = p[code]
        fc = free_cmax_nM(d)
        if model_ready(d):
            b = fc / (fc + d["ic50_nM"]) * 100
            print(f"{code:<8}{d['mechanism']:<22}{d['ic50_nM']:>8.0f}"
                  f"{fc:>10.2f}{b:>8.2f}%  {d['assay_temp_C']}")
        else:
            print(f"{code:<8}{d['mechanism']:<22}{'-':>8}{'-':>10}{'-':>9}  "
                  f"{d['data_quality']}")
    print(f"\nModel ready (excluding metabolites): {model_ready_codes(p)}")