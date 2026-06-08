"""
O'Hara-Rudy 2011 (ORd) Action Potential Model — CardioSafe Pediatric
Numerically stabilized implementation with psychiatric drug block module.

Reference: O'Hara et al. PLoS Comput Biol 2011; 7(5):e1002061
Units: mV, ms, uA/uF, mM, nM
"""

import numpy as np
import pandas as pd
from scipy.integrate import odeint
import warnings
from pathlib import Path
warnings.filterwarnings("ignore")

# Resolve paths relative to this file so the module works from any cwd
_SRC_DIR  = Path(__file__).resolve().parent        # …/src/
_ROOT_DIR = _SRC_DIR.parent                        # …/cardiosafe-pediatric/
_DATA_DIR = _ROOT_DIR / "data"
_DEFAULT_PARAMS = str(_DATA_DIR / "herg_master_params.csv")

# ── CONSTANTS ────────────────────────────────────────────────────────────────
R = 8314.0; T = 310.0; F = 96485.0

# Cell geometry
L=0.01; rad=0.0011
vcell = 1000*np.pi*rad**2*L
Ageo  = 2*np.pi*rad**2 + 2*np.pi*rad*L
Acap  = 2*Ageo
vmyo=0.68*vcell; vnsr=0.0552*vcell; vjsr=0.0048*vcell; vss=0.02*vcell

# ── GHK helper (regularized near V=0) ────────────────────────────────────────
def ghk(V, ci, co, z, regularize=True):
    """Goldman-Hodgkin-Katz current equation, regularized near V=0."""
    vf = z * V * F / (R * T)
    if regularize and abs(vf) < 1e-7:
        # L'Hopital limit: phi -> z*F*(ci - co)
        return z * F * (ci - co)
    ev = np.exp(vf)
    return z * V * F / (R * T) * (ci * ev - co) / (ev - 1.0)

# ── DRUG BLOCK MODULE ─────────────────────────────────────────────────────────
def load_drug_params(csv_path=None):
    if csv_path is None:
        csv_path = _DEFAULT_PARAMS
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

def compute_autonomic_modifiers(drug_combo, drug_params):
    hr_mult=1.0; gcal_mult=1.0; gna_mult=1.0
    for drug in drug_combo:
        if drug not in drug_params: continue
        mech = drug_params[drug]["primary_mechanism"]
        if mech == "sympathomimetic":
            hr_mult *= 0.90; gcal_mult *= 1.15
        elif mech == "autonomic_modulation":
            hr_mult *= 1.10; gna_mult  *= 0.95
    return hr_mult, gcal_mult, gna_mult

# ── STATE INDICES ─────────────────────────────────────────────────────────────
# 0:V 1:CaMKt 2:nai 3:nass 4:ki 5:kss 6:cai 7:cass 8:cansr 9:cajsr
# 10:m 11:hf 12:hs 13:j 14:hsp 15:jp
# 16:mL 17:hL 18:hLp
# 19:a 20:iF 21:iS 22:ap 23:iFp 24:iSp
# 25:d 26:ff 27:fs 28:fcaf 29:fcas 30:jca 31:nca 32:ffp 33:fcafp
# 34:xrf 35:xrs 36:xs1 37:xs2 38:xk1
# 39:Jrelnp 40:Jrelp

def ord_rhs(y, t, GKr_scale=1.0, GCaL_scale=1.0, GNa_scale=1.0,
            stim_amp=-80.0, stim_dur=0.5, CL=1000.0):

    (V, CaMKt, nai, nass, ki, kss, cai, cass, cansr, cajsr,
     m, hf, hs, j, hsp, jp,
     mL, hL, hLp,
     a, iF, iS, ap, iFp, iSp,
     d, ff, fs, fcaf, fcas, jca, nca, ffp, fcafp,
     xrf, xrs, xs1, xs2, xk1,
     Jrelnp, Jrelp) = y

    nao=140.0; cao=1.8; ko=5.4
    zna=1.0; zca=2.0; zk=1.0

    # CaMKII
    CaMKo=0.05; KmCaM=0.0015; KmCaMK=0.15; aCaMK=0.05; bCaMK=0.00068
    CaMKb = CaMKo*(1-CaMKt)/(1+KmCaM/max(cass,1e-12))
    CaMKa = CaMKb + CaMKt
    dCaMKt = aCaMK*CaMKb*(CaMKb+CaMKt) - bCaMK*CaMKt
    fCaMKp = 1.0/(1+KmCaMK/max(CaMKa,1e-12))

    # Reversals
    ENa = R*T/F * np.log(nao/max(nai,1e-6))
    EK  = R*T/F * np.log(ko /max(ki ,1e-6))
    EKs = R*T/F * np.log((ko+0.01833*nao)/max(ki+0.01833*nai,1e-6))

    # ── INa ──────────────────────────────────────────────────────────────────
    GNa = GNa_scale * 75.0
    mss = 1/(1+np.exp(-(V+39.57)/9.871))
    tm  = 1/(6.765*np.exp((V+11.64)/34.77)+8.552*np.exp(-(V+77.42)/5.955))
    dm  = (mss-m)/tm
    hss = 1/(1+np.exp((V+82.9)/6.086))
    thf = 1/(1.416e-5*np.exp(-(V+14.78)/10.52)+2.045e-7*np.exp((V+13.89)/6.085))
    ths = 1/(0.009794*np.exp(-(V+17.95)/28.05)+0.3343*np.exp((V+5.73)/56.66))
    Ahf=0.99; Ahs=0.01
    dhf=(hss-hf)/thf; dhs=(hss-hs)/ths
    h = Ahf*hf+Ahs*hs
    jss=hss
    tj=2.038+1/(0.02136*np.exp(-(V+100.6)/8.281)+0.3052*np.exp((V+0.9941)/38.45))
    dj=(jss-j)/tj
    hssp=1/(1+np.exp((V+89.1)/6.086))
    thsp=3*ths; dhsp=(hssp-hsp)/thsp; hp=Ahf*hf+Ahs*hsp
    tjp=1.46*tj; djp=(jss-jp)/tjp
    fINap=fCaMKp
    INa = GNa*(V-ENa)*m**3*((1-fINap)*h*j + fINap*hp*jp)

    # ── INaL ─────────────────────────────────────────────────────────────────
    GNaL=GNa_scale*0.0075
    mLss=1/(1+np.exp(-(V+42.85)/5.264)); tmL=tm; dmL=(mLss-mL)/tmL
    hLss=1/(1+np.exp((V+87.61)/7.488)); thL=200.0; dhL=(hLss-hL)/thL
    hLssp=1/(1+np.exp((V+93.81)/7.488)); thLp=3*thL; dhLp=(hLssp-hLp)/thLp
    fINaLp=fCaMKp
    INaL=GNaL*(V-ENa)*mL*((1-fINaLp)*hL+fINaLp*hLp)

    # ── Ito ───────────────────────────────────────────────────────────────────
    Gto=0.02
    ass=1/(1+np.exp(-(V-14.34)/14.82))
    ta=1.0515/(1/(1.2089*(1+np.exp(-(V-18.41)/29.38)))+3.5/(1+np.exp((V+100)/29.38)))
    da=(ass-a)/ta
    iFss=1/(1+np.exp((V+43.94)/5.711))
    tiF=4.562+1/(0.3933*np.exp(-(V+100)/100)+0.08004*np.exp((V+50)/16.59))
    iSss=iFss
    tiS=23.62+1/(0.001416*np.exp(-(V+96.52)/59.05)+1.780e-8*np.exp((V+114.1)/8.079))
    AiF=1/(1+np.exp((V-213.6)/151.2)); AiS=1-AiF
    diF=(iFss-iF)/tiF; diS=(iSss-iS)/tiS
    i_=AiF*iF+AiS*iS
    assp=1/(1+np.exp(-(V-24.34)/14.82)); dap=(assp-ap)/ta
    diFp=(iFss-iFp)/(2.5*tiF); diSp=(iSss-iSp)/(2.5*tiS)
    ip=AiF*iFp+AiS*iSp
    fItop=fCaMKp
    Ito=Gto*(V-EK)*((1-fItop)*i_*a+fItop*ip*ap)

    # ── ICaL (GHK, regularized) ───────────────────────────────────────────────
    PCa=GCaL_scale*0.0001; PCap=1.1*PCa
    PCaNa=0.00125*PCa; PCaK=3.574e-4*PCa
    dss=1/(1+np.exp(-(V+3.94)/4.23))
    td=0.6+1/(np.exp(-0.05*(V+6))+np.exp(0.09*(V+14)))
    dd=(dss-d)/td
    fss=1/(1+np.exp((V+19.58)/3.696))
    tff=7+1/(0.0045*np.exp(-(V+20)/10)+0.0045*np.exp((V+20)/10))
    tfs=1000+1/(3.5e-5*np.exp(-(V+5)/4)+3.5e-5*np.exp((V+5)/6))
    Aff=0.6; Afs=0.4
    dff=(fss-ff)/tff; dfs=(fss-fs)/tfs
    f=Aff*ff+Afs*fs
    fcass=fss
    tfcaf=7+1/(0.04*np.exp(-(V-4)/7)+0.04*np.exp((V-4)/7))
    tfcas=100+1/(0.00012*np.exp(-V/3)+0.00012*np.exp(V/7))
    Afcaf=0.3+0.6/(1+np.exp((V-10)/10)); Afcas=1-Afcaf
    dfcaf=(fcass-fcaf)/tfcaf; dfcas=(fcass-fcas)/tfcas
    fca=Afcaf*fcaf+Afcas*fcas
    tjca=75.0; djca=(fcass-jca)/tjca
    Km=0.0006
    dnca=(10*cass**2/(10*cass**2+Km**2)-nca)/tjca
    dffp=(fss-ffp)/(2.5*tff); fp=Aff*ffp+Afs*fs
    dfcafp=(fcass-fcafp)/(2.5*tfcaf); fcap=Afcaf*fcafp+Afcas*fcas
    fICaLp=fCaMKp

    # GHK terms (regularized)
    vf2  = 2*V*F/(R*T)
    ev2  = np.exp(vf2)
    denom2 = ev2-1.0
    if abs(denom2) < 1e-10: denom2 = 1e-10
    PhiCaL = 4*V*F**2/(R*T)*(cass*ev2-0.341*cao)/denom2

    vf1  = V*F/(R*T)
    ev1  = np.exp(vf1)
    denom1 = ev1-1.0
    if abs(denom1) < 1e-10: denom1 = 1e-10
    PhiCaNa = V*F**2/(R*T)*(0.75*nass*ev1-0.75*nao)/denom1
    PhiCaK  = V*F**2/(R*T)*(0.75*kss*ev1-0.75*ko)/denom1

    gd = d*((1-fICaLp)*(f*(1-nca)+jca*fca*nca) + fICaLp*(fp*(1-nca)+jca*fcap*nca))
    ICaL  = PCa  * PhiCaL  * gd
    ICaNa = PCaNa* PhiCaNa * gd
    ICaK  = PCaK * PhiCaK  * gd

    # ── IKr ──────────────────────────────────────────────────────────────────
    GKr = GKr_scale * 0.046
    xrss=1/(1+np.exp(-(V+8.337)/6.789))
    txrf=12.98+1/(0.3652*np.exp((V-31.66)/3.869)+4.123e-5*np.exp(-(V-47.78)/20.38))
    txrs=1.865+1/(0.06629*np.exp((V-34.70)/7.355)+1.128e-5*np.exp(-(V-29.74)/25.94))
    Axrf=1/(1+np.exp((V+54.81)/38.21)); Axrs=1-Axrf
    dxrf=(xrss-xrf)/txrf; dxrs=(xrss-xrs)/txrs
    Xr=Axrf*xrf+Axrs*xrs
    rkr=1/(1+np.exp((V+55)/75))*1/(1+np.exp((V-10)/30))
    IKr=GKr*np.sqrt(ko/5.4)*Xr*rkr*(V-EK)

    # ── IKs ──────────────────────────────────────────────────────────────────
    GKs=0.0034
    xs1ss=1/(1+np.exp(-(V+11.60)/8.932))
    txs1=817.3+1/(2.326e-4*np.exp((V+48.28)/17.80)+0.001292*np.exp(-(V+210)/230))
    dxs1=(xs1ss-xs1)/txs1
    xs2ss=xs1ss
    txs2=1/(0.01*np.exp((V-50)/20)+0.0193*np.exp(-(V+66.54)/31))
    dxs2=(xs2ss-xs2)/txs2
    KsCa=1+0.6/(1+(3.8e-5/max(cai,1e-9))**1.4)
    IKs=GKs*KsCa*xs1*xs2*(V-EKs)

    # ── IK1 ──────────────────────────────────────────────────────────────────
    GK1=0.1908
    xk1ss=1/(1+np.exp(-(V+2.5538*ko+144.59)/(1.5692*ko+3.8115)))
    txk1=122.2/(np.exp(-(V+127.2)/20.36)+np.exp((V+236.8)/69.33))
    dxk1=(xk1ss-xk1)/txk1
    rk1=1/(1+np.exp((V+105.8-2.6*ko)/9.493))
    IK1=GK1*np.sqrt(ko)*rk1*xk1*(V-EK)

    # ── INaCa (simplified Luo-Rudy style for stability) ──────────────────────
    Gncx=0.0008; KmCaN=0.0025; KmCaAct=150e-6
    alpha_ncx=2.5; gamma_ncx=0.35
    num_ncx = (np.exp(gamma_ncx*V*F/(R*T))*nai**3*cao -
               np.exp((gamma_ncx-1)*V*F/(R*T))*nao**3*cai*2.5)
    denom_ncx = (5000 + 1*(87.5**3+nao**3)*(0.0013+cao)*(1+0.1*np.exp((gamma_ncx-1)*V*F/(R*T))))
    INaCa = Gncx * num_ncx / max(denom_ncx, 1e-10)

    # ── INaK ─────────────────────────────────────────────────────────────────
    Pnak=30.0; KmNai=87.5; KmKo=0.5
    INaK = Pnak*(ko/(ko+KmKo))*(nai**1.5/(nai**1.5+KmNai**1.5)) * \
           (V+150)/(V+200) * 1/(1+0.1245*np.exp(-0.1*V*F/(R*T))+0.0365*np.exp(-V*F/(R*T)))

    # ── Background currents ───────────────────────────────────────────────────
    GKb=0.003; xkbss=1/(1+np.exp(-(V-14.48)/18.34))
    IKb=GKb*xkbss*(V-EK)

    PNab=3.75e-10
    INab=PNab*ghk(V,nai,nao,1)

    PCab=2.5e-8
    ICab=PCab*ghk(V,cai,0.341*cao,2)

    GpCa=0.0005
    IpCa=GpCa*cai/(0.0005+cai)

    # ── Ca handling ───────────────────────────────────────────────────────────
    Jupnp=0.004375*cai/(cai+0.00092)
    Jupp=2.75*0.004375*cai/(cai+0.00075)
    fJupp=fCaMKp; Jup=(1-fJupp)*Jupnp+fJupp*Jupp
    Jleak=0.0039375*cansr/15
    Jtr=(cansr-cajsr)/60

    Jrel_inf=(-ICaL)*1.0/(1+(0.3/max(cajsr,1e-9))**8)
    tau_rel=max(0.001,0.5/(1+(0.3/max(cajsr,1e-9))**8))
    Jrel_inf_p=Jrel_inf*1.7; tau_rel_p=max(0.001,tau_rel/1.7)
    dJrelnp=(Jrel_inf-Jrelnp)/tau_rel
    dJrelp=(Jrel_inf_p-Jrelp)/tau_rel_p
    Jrel=(1-fCaMKp)*Jrelnp+fCaMKp*Jrelp

    Jdiff=(cass-cai)/0.2; JdiffNa=(nass-nai)/2; JdiffK=(kss-ki)/2

    cmdnmax=0.05; kmcmdn=0.00238
    Bcai=1/(1+cmdnmax*kmcmdn/(cai+kmcmdn)**2)
    BSRmax=0.047; KmBSR=0.00087; BSLmax=1.124; KmBSL=0.0087
    Bcass=1/(1+BSRmax*KmBSR/(cass+KmBSR)**2+BSLmax*KmBSL/(cass+KmBSL)**2)
    csqnmax=10.0; kmcsqn=0.8
    Bcajsr=1/(1+csqnmax*kmcsqn/(cajsr+kmcsqn)**2)

    # ── Concentration ODEs ────────────────────────────────────────────────────
    dnai  = -(INa+INaL+3*INaCa+3*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo
    dnass = -(ICaNa)*Acap/(F*vss)-JdiffNa
    dki   = -(Ito+IKr+IKs+IK1+IKb-2*INaK+ICaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo
    dkss  = -(ICaK)*Acap/(F*vss)-JdiffK
    dcai  = Bcai*(-(IpCa+ICab-2*INaCa)*Acap/(2*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo)
    dcass = Bcass*(-(ICaL)*Acap/(2*F*vss)+Jrel*vjsr/vss-Jdiff)
    dcansr= Jup-Jleak-Jtr*vjsr/vnsr
    dcajsr= Bcajsr*(Jtr-Jrel)

    # ── Stimulus ──────────────────────────────────────────────────────────────
    t_in_CL = t % CL
    Istim = stim_amp if t_in_CL < stim_dur else 0.0

    # ── Voltage ODE ───────────────────────────────────────────────────────────
    Itot = INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa+INaK+INab+ICab+IpCa+IKb
    dV   = -(Itot+Istim)

    return [dV, dCaMKt,
            dnai, dnass, dki, dkss, dcai, dcass, dcansr, dcajsr,
            dm, dhf, dhs, dj, dhsp, djp,
            dmL, dhL, dhLp,
            da, diF, diS, dap, diFp, diSp,
            dd, dff, dfs, dfcaf, dfcas, djca, dnca, dffp, dfcafp,
            dxrf, dxrs, dxs1, dxs2, dxk1,
            dJrelnp, dJrelp]


# ── INITIAL CONDITIONS ────────────────────────────────────────────────────────
Y0 = np.array([
    -87.5, 0.0, 7.0, 7.0, 145.0, 145.0, 1e-4, 1e-4, 1.5, 1.5,
    0.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    0.0, 1.0, 1.0,
    0.0, 1.0, 1.0, 0.0, 1.0, 1.0,
    0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0,
], dtype=float)

# ── SIMULATION RUNNER ─────────────────────────────────────────────────────────
def run_simulation(drug_combination=None, drug_params_path=None,
                   n_beats=100, CL=1000.0, stim_amp=-80.0, stim_dur=0.5, verbose=True):

    if drug_params_path is None:
        drug_params_path = _DEFAULT_PARAMS
    drug_params = load_drug_params(drug_params_path)

    if drug_combination:
        total_ikr_block, ikr_breakdown = compute_ikr_block(drug_combination, drug_params)
        GKr_scale = 1.0 - total_ikr_block
        hr_mult, gcal_mult, gna_mult = compute_autonomic_modifiers(drug_combination, drug_params)
        CL_eff = CL * hr_mult
    else:
        GKr_scale=1.0; gcal_mult=1.0; gna_mult=1.0; CL_eff=CL
        total_ikr_block=0.0; ikr_breakdown={}

    if verbose:
        drugs_str = ", ".join(drug_combination.keys()) if drug_combination else "Baseline"
        print(f"\n{'='*60}")
        print(f"Simulation: {drugs_str}")
        print(f"  GKr scale: {GKr_scale:.4f}  (IKr block: {total_ikr_block*100:.2f}%)")
        print(f"  CL_eff:    {CL_eff:.1f} ms")

    y = Y0.copy()
    last_t = last_y = None

    for beat in range(n_beats):
        t0 = beat * CL_eff
        t1 = t0 + CL_eff
        n_pts = 2000 if beat == n_beats-1 else 200
        t_span = np.linspace(t0, t1, n_pts)

        sol = odeint(
            ord_rhs, y, t_span,
            args=(GKr_scale, gcal_mult, gna_mult, stim_amp, stim_dur, CL_eff),
            rtol=1e-6, atol=1e-8,
            full_output=False, mxstep=5000
        )

        if not np.all(np.isfinite(sol[-1])):
            if verbose:
                print(f"  Warning: NaN at beat {beat}, stopping.")
            break

        y = sol[-1]
        if beat == n_beats-1:
            last_t = t_span - t_span[0]
            last_y = sol.T  # shape (41, n_pts)

    if last_y is None:
        return {"APD90": np.nan, "QTc": np.nan}

    V_trace = last_y[0]
    t_trace = last_t

    V_max  = np.max(V_trace)
    # Use minimum voltage in second half of beat as resting estimate
    half   = len(V_trace)//2
    V_rest = np.min(V_trace[half:])
    V_90   = V_rest + 0.10*(V_max - V_rest)
    idx_up = int(np.argmax(V_trace))
    post   = V_trace[idx_up:]
    t_post = t_trace[idx_up:]
    cross  = np.where(post <= V_90)[0]

    APD90 = t_post[cross[0]] if len(cross) else np.nan
    QTc   = APD90 / np.sqrt(CL_eff/1000) if not np.isnan(APD90) else np.nan

    if verbose:
        print(f"  APD90: {APD90:.1f} ms  |  QTc: {QTc:.1f} ms  |  Vpeak: {V_max:.1f} mV")

    return {"APD90": APD90, "QTc": QTc,
            "IKr_block_pct": total_ikr_block*100,
            "GKr_scale": GKr_scale, "CL_effective": CL_eff,
            "drug_breakdown": ikr_breakdown,
            "t_trace": t_trace, "V_trace": V_trace}


def run_polypharmacy_sweep(combinations, drug_params_path=None,
                           n_beats=100, CL=1000.0):
    if drug_params_path is None:
        drug_params_path = _DEFAULT_PARAMS
    baseline = run_simulation(None, drug_params_path, n_beats, CL, verbose=False)
    bAPD = baseline["APD90"]; bQTc = baseline["QTc"]
    print(f"Baseline  APD90={bAPD:.1f} ms  QTc={bQTc:.1f} ms\n")

    rows = []
    for combo in combinations:
        res  = run_simulation(combo, drug_params_path, n_beats, CL, verbose=False)
        name = " + ".join(combo.keys())
        dAPD = res["APD90"]-bAPD; dQTc = res["QTc"]-bQTc
        rows.append({
            "combination":    name,
            "n_drugs":        len(combo),
            "APD90_ms":       round(res["APD90"],1),
            "QTc_ms":         round(res["QTc"],1),
            "ΔAPD90_ms":      round(dAPD,1),
            "ΔQTc_ms":        round(dQTc,1),
            "IKr_block_pct":  round(res["IKr_block_pct"],2),
        })
        risk = "⚠ HIGH" if dQTc>20 else ("△ MOD" if dQTc>10 else "  OK")
        print(f"  {risk} {name[:50]:50s} APD90:{res['APD90']:.1f} QTc:{res['QTc']:.1f} ΔQTc:{dQTc:+.1f}")

    df = pd.DataFrame(rows).sort_values("ΔQTc_ms", ascending=False)
    return df, baseline


# ── MAIN ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import time
    print("CardioSafe Pediatric — O'Hara-Rudy Simulator")
    print("="*60)

    PARAMS = _DEFAULT_PARAMS
    N = 50  # beats (increase to 200+ for true steady state)

    t0 = time.time()
    baseline = run_simulation(None, PARAMS, n_beats=N, verbose=True)
    print(f"  Runtime: {time.time()-t0:.1f}s")

    print("\n── POLYPHARMACY SWEEP ──")
    combos = [
        {"Methylphenidate":"therapeutic","Aripiprazole":"therapeutic"},
        {"Methylphenidate":"therapeutic","Sertraline":"therapeutic"},
        {"Risperidone":"therapeutic","Sertraline":"therapeutic"},
        {"Risperidone":"therapeutic","Fluoxetine":"therapeutic"},
        {"Methylphenidate":"therapeutic","Risperidone":"therapeutic","Sertraline":"therapeutic"},
        {"Methylphenidate":"therapeutic","Aripiprazole":"therapeutic","Fluoxetine":"therapeutic"},
        {"Methylphenidate":"therapeutic","Clonidine":"therapeutic"},
        {"Methylphenidate":"therapeutic","Risperidone":"therapeutic","Clonidine":"therapeutic"},
        {"Imipramine":"therapeutic","Methylphenidate":"therapeutic"},
        {"Nortriptyline":"therapeutic","Risperidone":"therapeutic"},
    ]

    t1 = time.time()
    df, _ = run_polypharmacy_sweep(combos, PARAMS, n_beats=N, CL=1000.0)
    print(f"\nSweep runtime: {time.time()-t1:.1f}s")

    print("\n── RESULTS (sorted by ΔQTc) ──")
    print(df.to_string(index=False))
    df.to_csv("polypharmacy_sweep.csv", index=False)
    print("\nSaved: polypharmacy_sweep.csv")