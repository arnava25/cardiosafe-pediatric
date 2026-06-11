"""
Mechanism-stratified concordance analysis — addition to concordance_stats.py
Run after concordance_stats.py has produced faers_model_alignment.csv
"""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import fisher_exact

_ROOT_DIR  = Path(__file__).resolve().parent.parent
_FAERS_DIR = _ROOT_DIR / "results" / "faers"
np.random.seed(42)

def compute_ror(a,b,c,d):
    a,b,c,d = a+0.5,b+0.5,c+0.5,d+0.5
    ror = (a/b)/(c/d)
    se  = np.sqrt(1/a+1/b+1/c+1/d)
    return round(ror,3), round(np.exp(np.log(ror)-1.96*se),3), round(np.exp(np.log(ror)+1.96*se),3)

def metrics(tp,fp,fn,tn):
    total = tp+fp+fn+tn
    sens = tp/(tp+fn) if (tp+fn)>0 else np.nan
    spec = tn/(tn+fp) if (tn+fp)>0 else np.nan
    ppv  = tp/(tp+fp) if (tp+fp)>0 else np.nan
    f1   = 2*tp/(2*tp+fp+fn) if (2*tp+fp+fn)>0 else np.nan
    p_o  = (tp+tn)/total
    p_e  = (((tp+fp)/total)*((tp+fn)/total)+(( fn+tn)/total)*((fp+tn)/total))
    kappa= (p_o-p_e)/(1-p_e) if p_e<1 else np.nan
    _,pf = fisher_exact([[tp,fp],[fn,tn]], alternative='greater')
    return dict(TP=tp,FP=fp,FN=fn,TN=tn,
                sensitivity=round(sens,3),specificity=round(spec,3),
                PPV=round(ppv,3),F1=round(f1,3),kappa=round(kappa,3),
                fisher_p=round(pf,4))

def permutation_p(model_pos, faers_pos, obs_kappa, n=10000):
    hits = 0
    for _ in range(n):
        s = np.random.permutation(model_pos)
        tp=int((s&faers_pos).sum()); fp=int((s&~faers_pos).sum())
        fn=int((~s&faers_pos).sum()); tn=int((~s&~faers_pos).sum())
        tot=tp+fp+fn+tn
        po=(tp+tn)/tot
        pe=(((tp+fp)/tot)*((tp+fn)/tot)+(( fn+tn)/tot)*((fp+tn)/tot))
        k=(po-pe)/(1-pe) if pe<1 else 0
        if k>=obs_kappa: hits+=1
    return hits/n

def run():
    df = pd.read_csv(_FAERS_DIR / "faers_model_alignment.csv")
    df = df.dropna(subset=["faers_ROR"]).copy()

    STIMS = {"MPH","AMP"}
    def has_stim(combo):
        parts = combo.split("+")
        return any(p in STIMS for p in parts)

    df["has_stimulant"] = df["combination"].apply(has_stim)
    df["model_pos"]     = df["model_risk_tier"].isin(["HIGH","MODERATE"])
    df["faers_pos"]     = df["faers_signal"].astype(bool)

    strat_labels = {
        "ALL pairs":              df,
        "Stimulant-containing":   df[df["has_stimulant"]],
        "Non-stimulant only":     df[~df["has_stimulant"]],
    }

    print("="*70)
    print("MECHANISM-STRATIFIED CONCORDANCE — CardioSafe Pediatric")
    print("="*70)

    rows = []
    for label, sub in strat_labels.items():
        n = len(sub)
        mp = sub["model_pos"].values
        fp = sub["faers_pos"].values

        tp=int((mp&fp).sum()); fpos=int((mp&~fp).sum())
        fn=int((~mp&fp).sum()); tn=int((~mp&~fp).sum())
        m  = metrics(tp,fpos,fn,tn)
        pp = permutation_p(mp,fp,m["kappa"])

        n_faers_sig = int(fp.sum())
        n_model_pos = int(mp.sum())

        print(f"\n── {label}  (n={n}, FAERS signals={n_faers_sig}, model positive={n_model_pos}) ──")
        print(f"  Confusion:  TP={tp}  FP={fpos}  FN={fn}  TN={tn}")
        print(f"  Sensitivity:  {m['sensitivity']:.3f}   Specificity: {m['specificity']:.3f}")
        print(f"  PPV:          {m['PPV']:.3f}   F1:          {m['F1']:.3f}")
        print(f"  Cohen kappa:  {m['kappa']:.3f}   Fisher p:    {m['fisher_p']:.4f}")
        print(f"  Permutation p: {pp:.4f}  {'*** SIGNIFICANT' if pp<0.05 else '(not significant)'}")

        rows.append({"stratum":label,"n":n,"n_faers_signal":n_faers_sig,
                     "n_model_pos":n_model_pos,**m,"permutation_p":round(pp,4)})

    # Per-mechanism detail for non-stimulant
    print("\n── Non-stimulant pairs — full detail ──\n")
    ns = df[~df["has_stimulant"]].sort_values("model_dQTc_ms",ascending=False)
    print(f"  {'Combination':38s}  {'dQTc':>6}  {'Tier':10}  {'FAERS':5}  {'Match':5}  ROR")
    print("  "+"-"*75)
    for _,row in ns.iterrows():
        sig  = "YES" if row["faers_pos"] else "no"
        conc = "✓" if row["concordant"] else "✗"
        ror  = f"{row['faers_ROR']:.2f}" if not pd.isna(row["faers_ROR"]) else "n/a"
        print(f"  {row['combination']:38s}  {row['model_dQTc_ms']:+6.1f}  "
              f"{row['model_risk_tier']:10}  {sig:5}  {conc:5}  {ror}")

    # Save
    out = _FAERS_DIR / "concordance_stratified.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"\nSaved → {out}")

if __name__ == "__main__":
    run()
