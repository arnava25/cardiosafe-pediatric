"""
CardioSafe Pediatric — drug module on the validated canonical O'Hara-Rudy core.
Drop-in replacement for the old ord_model.py: same run_simulation interface,
same compute_ikr_block / compute_autonomic_modifiers, but the cell engine is
now ord_core (faithful INaK / two-compartment NCX / SR cycling), which converges
to a physiological steady state.
"""
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from scipy.integrate import odeint
import ord_core as C   # the validated engine

def load_drug_params(csv_path="herg_master_params.csv"):
    try:
        df = pd.read_csv(csv_path)
        return df.set_index("drug_name").to_dict("index")
    except FileNotFoundError:
        print(f"Warning: {csv_path} not found.")
        return {}

def compute_ikr_block(drug_combo, drug_params):
    blocks = {}
    for drug, val in drug_combo.items():
        if drug not in drug_params:
            continue
        dp = drug_params[drug]
        if dp["primary_mechanism"] != "hERG_block":
            blocks[drug] = 0.0; continue
        ic50 = float(dp["ic50_nM"])
        conc = float(dp["cmax_free_nM"]) if val == "therapeutic" else float(val)
        blocks[drug] = conc / (conc + ic50)
    total = 1.0 - np.prod([1.0-b for b in blocks.values()]) if blocks else 0.0
    return total, blocks

IKS_UPREG = 1.20   # per-sympathomimetic-drug IKs upregulation; set 1.0 to disable

def compute_autonomic_modifiers(drug_combo, drug_params):
    hr_mult=1.0; gcal_mult=1.0; gna_mult=1.0; gks_mult=1.0
    for drug in drug_combo:
        if drug not in drug_params: continue
        mech = drug_params[drug]["primary_mechanism"]
        if mech == "sympathomimetic":
            hr_mult *= 0.90; gcal_mult *= 1.15; gks_mult *= IKS_UPREG
        elif mech == "autonomic_modulation":
            hr_mult *= 1.10; gna_mult  *= 0.95
    return hr_mult, gcal_mult, gna_mult, gks_mult

def _apd90_qtc(V, t, CL_eff):
    Vmax=np.max(V); Vrest=np.min(V[len(V)//2:])
    V90=Vrest+0.10*(Vmax-Vrest); iu=int(np.argmax(V))
    post=V[iu:]; tp=t[iu:]-t[iu]; cr=np.where(post<=V90)[0]
    apd=(tp[cr[0]] if len(cr) else np.nan)
    qtc=apd/np.sqrt(CL_eff/1000) if not np.isnan(apd) else np.nan
    return apd, qtc, Vmax

def run_simulation(drug_combination=None, drug_params_path="herg_master_params.csv",
                   n_beats=200, CL=1000.0, stim_amp=-80.0, stim_dur=0.5, verbose=True):
    drug_params = load_drug_params(drug_params_path)
    if drug_combination:
        total_ikr_block, ikr_breakdown = compute_ikr_block(drug_combination, drug_params)
        GKr_mult = 1.0 - total_ikr_block
        hr_mult, gcal_mult, gna_mult, gks_mult = compute_autonomic_modifiers(drug_combination, drug_params)
        CL_eff = CL * hr_mult
    else:
        GKr_mult=1.0; gcal_mult=1.0; gna_mult=1.0; gks_mult=1.0; CL_eff=CL
        total_ikr_block=0.0; ikr_breakdown={}

    if verbose:
        ds=", ".join(drug_combination.keys()) if drug_combination else "Baseline"
        print(f"Simulation: {ds}  GKr={GKr_mult:.4f} (IKr block {total_ikr_block*100:.2f}%)  CL_eff={CL_eff:.1f}")

    y=C.Y0.copy(); last_t=last_y=None
    for beat in range(n_beats):
        t0=beat*CL_eff; n_pts=2000 if beat==n_beats-1 else 300
        ts=np.linspace(t0, t0+CL_eff, n_pts)
        sol=odeint(C.rhs, y, ts,
                   args=(GKr_mult, gcal_mult, gna_mult, gks_mult, stim_amp, stim_dur, CL_eff),
                   rtol=1e-6, atol=1e-8, mxstep=8000)
        if not np.all(np.isfinite(sol[-1])):
            if verbose: print(f"  NaN at beat {beat}"); 
            break
        y=sol[-1]
        if beat==n_beats-1:
            last_t=ts-ts[0]; last_y=sol
    if last_y is None:
        return {"APD90":np.nan,"QTc":np.nan,"IKr_block_pct":total_ikr_block*100}
    V=last_y[:,0]
    apd, qtc, vpeak = _apd90_qtc(V, last_t, CL_eff)
    if verbose:
        print(f"  APD90: {apd:.1f} ms  |  QTc: {qtc:.1f} ms  |  Vpeak: {vpeak:.1f} mV")
    return {"APD90":apd, "QTc":qtc, "Vpeak":vpeak, "IKr_block_pct":total_ikr_block*100,
            "GKr_scale":GKr_mult, "CL_effective":CL_eff, "drug_breakdown":ikr_breakdown,
            "t_trace":last_t, "V_trace":V}

def run_polypharmacy_sweep(combinations, drug_params_path="herg_master_params.csv",
                           n_beats=200, CL=1000.0):
    base=run_simulation(None, drug_params_path, n_beats, CL, verbose=False)
    bAPD=base["APD90"]; bQTc=base["QTc"]
    print(f"Baseline APD90={bAPD:.1f} QTc={bQTc:.1f}\n")
    rows=[]
    for combo in combinations:
        r=run_simulation(combo, drug_params_path, n_beats, CL, verbose=False)
        name=" + ".join(combo.keys()); dQTc=r["QTc"]-bQTc
        rows.append({"combination":name,"n_drugs":len(combo),
                     "APD90_ms":round(r["APD90"],1),"QTc_ms":round(r["QTc"],1),
                     "dQTc_ms":round(dQTc,1),"IKr_block_pct":round(r["IKr_block_pct"],2)})
        tag="HIGH" if dQTc>20 else ("MOD" if dQTc>10 else "ok")
        print(f"  {tag:4s} {name[:48]:48s} dQTc={dQTc:+.1f}")
    return pd.DataFrame(rows).sort_values("dQTc_ms",ascending=False), base