#!/usr/bin/env python3
"""
CardioSafe Pediatric - dual-metric FAERS concordance.

Scores BOTH genuine repolarization (ΔAPD90) and rate-confounded Bazett (ΔQTc)
against the FAERS signal, by subset. Tells you whether the real-world signal
tracks true APD prolongation or the heart-rate artifact.

NOT A COMPOSITE (verified July 2026). This computes three separate
Mann-Whitney AUCs against the FAERS signal — genuine ΔAPD90, Bazett ΔQTc, and
IKr block — reported side by side and never combined into a single number. It
is a discriminator, not a composite. Retained on that basis when the composite
score and its concordance variants were archived; see archive/src/.

Usage:
    python3 src/concordance_metrics.py \
        --grid results/risk_grid_results.csv \
        --faers results/faers/faers_combo_ror.csv
"""
import argparse, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata

ABBREV = {"methylphenidate":"MPH","amphetamine":"AMP","risperidone":"RIS",
 "quetiapine":"QUE","aripiprazole":"ARI","sertraline":"SER","fluoxetine":"FLU",
 "escitalopram":"ESC","clonidine":"CLO","guanfacine":"GUA","imipramine":"IMI",
 "nortriptyline":"NOR"}
KNOWN=set(ABBREV.values()); STIM={"MPH","AMP"}

def to_abbrev(t):
    t=t.strip().lower()
    return t.upper() if t.upper() in KNOWN else ABBREV.get(t,t.strip().upper()[:3])
def norm_combo(s):
    parts=[p for p in str(s).replace("+"," ").replace(","," ").split() if p]
    return "+".join(sorted({to_abbrev(p) for p in parts}))
def find_col(df,cands,required=True,what=""):
    low={c.lower():c for c in df.columns}
    for c in cands:
        if c.lower() in low: return low[c.lower()]
    for c in df.columns:
        if any(x.lower() in c.lower() for x in cands): return c
    if required: sys.exit(f"ERROR no {what}; have {list(df.columns)}")
    return None
def auc_mw(scores,labels):
    s=np.asarray(scores,float); l=np.asarray(labels,int)
    pos=s[l==1]; neg=s[l==0]; n1,n0=len(pos),len(neg)
    if n1==0 or n0==0: return np.nan
    r=rankdata(np.concatenate([pos,neg]))
    return (r[:n1].sum()-n1*(n1+1)/2)/(n1*n0)

def load_grid(path):
    g=pd.read_csv(path)
    combo=find_col(g,["combination","combo","drugs","pair","name"],what="combo")
    dqtc=find_col(g,["ΔQTc","dqtc","delta_qtc","qtc_delta"],what="dQTc")
    dapd=find_col(g,["ΔAPD90","ΔAPD","dapd","delta_apd","apd_delta"],what="dAPD90")
    ikr=find_col(g,["ikr_block","ikr_pct","ikr"],required=False,what="ikr")
    g["key"]=g[combo].map(norm_combo)
    g=g.rename(columns={dqtc:"dqtc",dapd:"dapd"})
    if ikr: g=g.rename(columns={ikr:"ikr"})
    cols=["key","dqtc","dapd"]+(["ikr"] if ikr else [])
    return g[cols].drop_duplicates("key"), bool(ikr)

def load_faers(path):
    f=pd.read_csv(path)
    combo=find_col(f,["combination","combo","drugs","pair","name"],what="combo")
    ror=find_col(f,["ror","reporting_odds_ratio"],required=False,what="ror")
    sig=find_col(f,["signal","is_signal","ror_signal","flag"],required=False,what="signal")
    f["key"]=f[combo].map(norm_combo)
    if ror: f=f.rename(columns={ror:"ror"})
    if sig:
        f=f.rename(columns={sig:"signal"})
        f["signal"]=(f["signal"].astype(str).str.lower().isin(["1","true","yes","signal","y"])).astype(int)
    elif "ror" in f.columns:
        f["signal"]=(f["ror"]>=2.0).astype(int); print("NOTE: ROR>=2 signal proxy")
    keep=["key"]+[c for c in ["ror","signal"] if c in f.columns]
    return f[keep].drop_duplicates("key")

def block(label, df, has_ikr):
    print(f"\n{'='*60}\n{label}  (n={len(df)}, signals={int(df.signal.sum())})\n{'='*60}")
    if len(df)<3 or df.signal.nunique()<2:
        print("  too few / no signal variation"); return
    for name,col in [("ΔAPD90 (genuine repol)","dapd"),("ΔQTc Bazett (rate-conf.)","dqtc")]:
        auc=auc_mw(df[col],df.signal)
        rho=spearmanr(df[col],df.ror,nan_policy="omit")[0] if "ror" in df else np.nan
        print(f"  {name:26s} AUC={auc:.3f}   Spearman_vs_ROR={rho:+.3f}")
    if has_ikr:
        print(f"  {'IKr% block':26s} AUC={auc_mw(df.ikr,df.signal):.3f}")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--grid",required=True); ap.add_argument("--faers",required=True)
    a=ap.parse_args()
    grid,has_ikr=load_grid(a.grid); faers=load_faers(a.faers)
    m=grid.merge(faers,on="key",how="inner")
    if "signal" not in m: sys.exit("no FAERS signal column")
    m["is_stim"]=m.key.apply(lambda k: bool(set(k.split("+"))&STIM))
    print(f"\nMatched {len(m)} combos (grid {len(grid)} x FAERS {len(faers)}).")
    block("ALL", m, has_ikr)
    block("STIMULANT-CONTAINING", m[m.is_stim], has_ikr)
    block("NON-STIMULANT", m[~m.is_stim], has_ikr)
    print("\nRead it:")
    print("  - ΔAPD90 AUC high, ΔQTc AUC low  -> FAERS tracks genuine repolarization")
    print("  - ΔQTc AUC high, ΔAPD90 AUC low  -> FAERS tracks the HR/rate artifact")
    print("  - both low                       -> model doesn't predict the real-world signal")

if __name__=="__main__": main()
