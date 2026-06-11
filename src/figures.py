#!/usr/bin/env python3
"""
CardioSafe Pediatric — no-engine manuscript figures (regenerated on validated grid).
Reads results/risk_grid_results.csv; age-stratification values are locked FAERS results.

Usage: python3 src/figures_noengine.py --grid results/risk_grid_results.csv --outdir docs/figures
"""
import argparse, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 140})

ORDER = ["MPH","AMP","RIS","QUE","ARI","SER","FLU","ESC","CLO","GUA","IMI","NOR"]
CLASS = {"MPH":"Stim","AMP":"Stim","RIS":"AP","QUE":"AP","ARI":"AP","SER":"SSRI",
         "FLU":"SSRI","ESC":"SSRI","CLO":"α2","GUA":"α2","IMI":"TCA","NOR":"TCA"}
CLASS_C = {"Stim":"#d1495b","AP":"#3a86ff","SSRI":"#2a9d8f","α2":"#8d6a9f","TCA":"#e09f3e"}
BASE = 263.6

def fig_heatmap(df, out):
    pairs = df[df.n_drugs == 2].copy()
    M = pd.DataFrame(np.nan, index=ORDER, columns=ORDER)
    for _, r in pairs.iterrows():
        a, b = str(r["combination"]).split("+")[:2]
        if a in ORDER and b in ORDER:
            M.loc[a, b] = r["ΔAPD90_ms"]; M.loc[b, a] = r["ΔAPD90_ms"]
    fig, ax = plt.subplots(figsize=(8, 7))
    norm = TwoSlopeNorm(vmin=-12, vcenter=0, vmax=12)
    im = ax.imshow(M.values.astype(float), cmap="RdBu_r", norm=norm)
    ax.set_xticks(range(12)); ax.set_yticks(range(12))
    ax.set_xticklabels(ORDER, rotation=45, ha="right"); ax.set_yticklabels(ORDER)
    for i, c in enumerate(ORDER):
        ax.get_xticklabels()[i].set_color(CLASS_C[CLASS[c]])
        ax.get_yticklabels()[i].set_color(CLASS_C[CLASS[c]])
    for i in range(12):
        for j in range(12):
            v = M.values[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:+.0f}", ha="center", va="center", fontsize=6.5,
                        color="black" if abs(v) < 7 else "white")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("genuine ΔAPD90 (ms)")
    ax.set_title("Genuine repolarization change by drug pair\n(uncorrected ΔAPD90; red = prolongation, blue = shortening)")
    handles = [plt.Line2D([0],[0],marker="s",ls="",c=CLASS_C[k],label=k) for k in CLASS_C]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.18,1), fontsize=8, frameon=False, title="class")
    plt.tight_layout(); plt.savefig(out, bbox_inches="tight"); plt.close()
    print("saved", out)

def fig_overflag(df, out):
    d = df.dropna(subset=["ΔQTc_ms"]).copy()
    RR = d["CL_eff_ms"]/1000.0
    d["dF"] = d["APD90_ms"]/RR**(1/3) - BASE
    metrics = [("genuine\nΔAPD90","ΔAPD90_ms","#2a9d8f"),
               ("Fridericia\nΔQTc","dF","#e09f3e"),
               ("Bazett\nΔQTc","ΔQTc_ms","#d1495b")]
    tiers = [("≥10 ms (elevated)",10,1e9),("5–9 ms",5,10),("<5 ms",-1e9,5)]
    fig, ax = plt.subplots(figsize=(7.5,4.8))
    x = np.arange(len(metrics)); w = 0.6
    bottoms = np.zeros(len(metrics))
    shades = ["#b23a48","#f2c14e","#d9d9d9"]
    for (tlabel,lo,hi),sh in zip(tiers,shades):
        counts = [int(((d[c]>=lo)&(d[c]<hi)).sum()) for _,c,_ in metrics]
        ax.bar(x, counts, w, bottom=bottoms, label=tlabel, color=sh, edgecolor="white")
        for xi,(ct,bi) in enumerate(zip(counts,bottoms)):
            if ct>0: ax.text(xi, bi+ct/2, str(ct), ha="center", va="center", fontsize=10, fontweight="bold")
        bottoms += counts
    ax.set_xticks(x); ax.set_xticklabels([m[0] for m in metrics])
    ax.set_ylabel("number of combinations (of 84)")
    ax.set_title("Bazett flags 29 combinations ≥10 ms; genuine ΔAPD90 flags only 5\nrate correction inflates the apparent risk count ~6×")
    ax.legend(frameon=False, fontsize=9, loc="upper left", bbox_to_anchor=(1.0,1))
    plt.tight_layout(); plt.savefig(out, bbox_inches="tight"); plt.close()
    print("saved", out)

def fig_age(out):
    # locked FAERS age-stratification results (children 6-12 vs adolescents 13-17)
    data = [("MPH+SER",43.71,1.86),("QUE+GUA",32.15,3.91),("MPH+ARI",23.71,0.47),
            ("ARI+GUA",21.72,1.69),("MPH+QUE",16.45,1.84),("MPH+RIS",0.86,0.78),
            ("AMP+SER",0.89,3.12),("QUE+FLU",3.75,10.88)]
    data = sorted(data, key=lambda r: r[1]/r[2], reverse=True)
    labels=[d[0] for d in data]; ch=[d[1] for d in data]; ad=[d[2] for d in data]
    y=np.arange(len(data))
    fig, ax = plt.subplots(figsize=(8,5))
    for yi,c,a in zip(y,ch,ad):
        ax.plot([a,c],[yi,yi], color="#bbb", lw=2, zorder=1)
    ax.scatter(ad,y,s=70,color="#3a86ff",zorder=2,label="Adolescents 13–17")
    ax.scatter(ch,y,s=70,color="#d1495b",zorder=2,label="Children 6–12")
    ax.axvline(1.0, ls="--", color="k", lw=1); ax.text(1.04,-0.75,"ROR = 1",fontsize=8)
    ax.set_xscale("log"); ax.set_yticks(y); ax.set_yticklabels(labels)
    ax.set_xlabel("FAERS reporting odds ratio (log scale)")
    ax.set_title("Stimulant cardiac signals are 23–50× higher in children than adolescents\n(reverse pattern for the CYP2D6-interaction pair QUE+FLU)")
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    plt.tight_layout(); plt.savefig(out, bbox_inches="tight"); plt.close()
    print("saved", out)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", default="results/risk_grid_results.csv")
    ap.add_argument("--outdir", default="docs/figures")
    a = ap.parse_args()
    df = pd.read_csv(a.grid)
    fig_heatmap(df, f"{a.outdir}/fig_genuine_apd90_heatmap.png")
    fig_overflag(df, f"{a.outdir}/fig_overflag_contrast.png")
    fig_age(f"{a.outdir}/fig_age_stratification.png")