"""
concordance_nogap.py — CardioSafe Pediatric
============================================
Recomputes concordance statistics after removing the two known-gap clusters:
1. Guanfacine combinations (PR/conduction mechanism not modeled)
2. Fluoxetine CYP2D6 combinations (pharmacokinetic DDI not modeled)

Also computes ROC curve and AUC using composite score as continuous predictor.
"""

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

# ── FULL CONCORDANCE DATA from faers_model_alignment ─────────────────────────
# (drug_a, drug_b, dQTc, model_tier, faers_signal, ROR, gap_type)
# gap_type: None, "conduction", "cyp2d6"

ALL_PAIRS = [
    # Stimulant combinations
    ("MPH","AMP",  19.1, "MODERATE",  False, 0.45,  None),
    ("MPH","ARI",  18.1, "MODERATE",  True,  8.15,  None),
    ("AMP","ARI",  18.1, "MODERATE",  False, 0.69,  None),
    ("MPH","NOR",  14.8, "MODERATE",  False, None,  None),
    ("AMP","NOR",  14.8, "MODERATE",  False, 11.41, None),
    ("MPH","QUE",  11.0, "MODERATE",  True,  5.25,  None),
    ("AMP","QUE",  11.0, "MODERATE",  False, 1.53,  None),
    ("MPH","IMI",  10.5, "MODERATE",  False, 4.46,  None),
    ("AMP","IMI",  10.5, "MODERATE",  False, 7.90,  None),
    ("MPH","RIS",  10.0, "MODERATE",  False, 0.65,  None),
    ("AMP","RIS",  10.0, "MODERATE",  False, 1.28,  None),
    ("MPH","FLU",  10.0, "MODERATE",  False, 0.12,  None),
    ("AMP","FLU",  10.0, "MODERATE",  False, 1.37,  None),
    ("MPH","SER",   9.6, "LOW-MOD",   True,  12.79, None),
    ("AMP","SER",   9.6, "LOW-MOD",   True,  2.81,  None),
    ("MPH","ESC",   9.1, "LOW-MOD",   False, 1.10,  None),
    ("AMP","ESC",   9.1, "LOW-MOD",   False, 0.36,  None),
    ("MPH","CLO",   6.8, "LOW-MOD",   False, 1.27,  None),
    ("AMP","CLO",   6.8, "LOW-MOD",   False, 0.93,  None),
    ("MPH","GUA",   6.8, "LOW-MOD",   True,  3.40,  "conduction"),
    ("AMP","GUA",   6.8, "LOW-MOD",   False, 0.37,  "conduction"),

    # hERG-driven non-stimulant
    ("ARI","NOR",  15.0, "MODERATE",  False, 9.34,  None),
    ("QUE","ARI",  11.5, "MODERATE",  False, 1.60,  None),
    ("ARI","IMI",  11.0, "MODERATE",  True,  14.67, None),
    ("ARI","FLU",  10.0, "MODERATE",  False, 1.31,  "cyp2d6"),
    ("RIS","ARI",  10.0, "MODERATE",  True,  1.71,  None),
    ("ARI","SER",   9.5, "LOW-MOD",   True,  6.78,  None),
    ("ARI","ESC",   9.5, "LOW-MOD",   True,  2.31,  None),
    ("QUE","NOR",   8.0, "LOW-MOD",   False, 1.06,  None),
    ("ARI","CLO",   7.3, "LOW-MOD",   False, 0.42,  None),
    ("ARI","GUA",   7.3, "LOW-MOD",   True,  10.25, "conduction"),
    ("FLU","NOR",   6.5, "LOW-MOD",   False, 9.34,  None),
    ("RIS","NOR",   6.5, "LOW-MOD",   False, 14.67, None),
    ("SER","NOR",   6.5, "LOW-MOD",   False, 6.04,  None),
    ("ESC","NOR",   6.0, "LOW-MOD",   False, 3.24,  None),
    ("CLO","NOR",   4.1, "LOW",        False, 11.41, None),
    ("QUE","IMI",   3.5, "LOW",        False, 4.11,  None),
    ("QUE","FLU",   3.0, "LOW",        True,  9.18,  "cyp2d6"),
    ("RIS","QUE",   3.0, "LOW",        False, 1.28,  None),
    ("QUE","SER",   2.5, "LOW",        True,  5.84,  None),
    ("RIS","IMI",   2.5, "LOW",        False, 2.77,  None),
    ("FLU","IMI",   2.5, "LOW",        False, 3.11,  "cyp2d6"),
    ("SER","IMI",   2.0, "LOW",        True,  16.22, "cyp2d6"),
    ("ESC","IMI",   2.0, "LOW",        True,  28.01, "cyp2d6"),
    ("QUE","ESC",   2.0, "LOW",        False, 0.58,  None),
    ("SER","FLU",   1.5, "LOW",        False, 0.72,  None),
    ("RIS","SER",   1.5, "LOW",        False, 0.65,  None),
    ("RIS","FLU",   1.5, "LOW",        True,  3.69,  "cyp2d6"),
    ("RIS","ESC",   1.0, "LOW",        False, 0.31,  None),
    ("FLU","ESC",   1.0, "LOW",        False, 0.36,  None),
    ("SER","ESC",   0.5, "LOW",        False, 0.98,  None),
    ("QUE","CLO",  -0.1, "LOW",        False, 0.87,  None),
    ("QUE","GUA",  -0.1, "LOW",        True,  13.15, "conduction"),
    ("CLO","IMI",  -0.6, "LOW",        False, 6.85,  None),
    ("GUA","IMI",  -0.6, "LOW",        False, 7.90,  "conduction"),
    ("RIS","GUA",  -1.1, "LOW",        False, 1.93,  "conduction"),
    ("FLU","GUA",  -1.1, "LOW",        False, 0.28,  "conduction"),
    ("RIS","CLO",  -1.1, "LOW",        False, 1.83,  None),
    ("FLU","CLO",  -1.1, "LOW",        False, 0.12,  None),
    ("CLO","GUA",  -1.5, "LOW",        False, 0.26,  "conduction"),
    ("SER","CLO",  -1.7, "LOW",        False, 0.74,  None),
    ("ESC","GUA",  -1.7, "LOW",        False, 1.46,  "conduction"),
    ("ESC","CLO",  -1.7, "LOW",        False, 2.69,  None),
    ("SER","GUA",  -1.7, "LOW",        True,  2.58,  "conduction"),
]

COMPOSITE_SCORES = {
    "MPH+ARI": 72.2, "ARI+IMI": 61.0, "ARI+GUA": 57.2, "AMP+ARI": 53.7,
    "ARI+NOR": 51.0, "ARI+SER": 49.4, "ARI+FLU": 49.1, "MPH+AMP": 47.8,
    "MPH+SER": 46.9, "QUE+FLU": 44.6, "MPH+QUE": 44.0, "MPH+NOR": 42.3,
    "AMP+NOR": 42.3, "ARI+ESC": 39.7, "QUE+ARI": 39.0, "SER+IMI": 38.0,
    "ESC+IMI": 36.0, "QUE+GUA": 34.6, "MPH+RIS": 33.0, "ARI+CLO": 32.0,
    "AMP+QUE": 31.5, "MPH+IMI": 30.0, "AMP+IMI": 29.0, "MPH+FLU": 28.0,
    "AMP+FLU": 27.5, "RIS+FLU": 27.0, "AMP+RIS": 26.0, "ARI+SER": 49.4,
    "RIS+NOR": 25.0, "FLU+NOR": 24.0, "SER+NOR": 23.0, "ESC+NOR": 22.0,
    "RIS+ARI": 20.0, "QUE+SER": 19.0, "MPH+ESC": 18.0, "AMP+ESC": 17.0,
    "CLO+NOR": 16.0, "MPH+CLO": 15.0, "AMP+CLO": 14.5, "MPH+GUA": 14.0,
    "AMP+GUA": 13.0, "QUE+IMI": 12.0, "RIS+IMI": 11.0, "FLU+IMI": 10.5,
    "RIS+QUE": 10.0, "QUE+ESC": 9.0,  "SER+FLU": 8.5,  "RIS+SER": 8.0,
    "RIS+GUA": 7.5,  "FLU+GUA": 7.0,  "RIS+ESC": 6.5,  "FLU+ESC": 6.0,
    "SER+ESC": 5.5,  "QUE+CLO": 5.0,  "CLO+IMI": 4.5,  "GUA+IMI": 4.0,
    "RIS+CLO": 3.5,  "FLU+CLO": 3.0,  "CLO+GUA": 2.5,  "SER+CLO": 2.0,
    "ESC+GUA": 1.5,  "ESC+CLO": 1.0,  "SER+GUA": 0.5,
}

def get_composite(a, b, scores):
    for key in [f"{a}+{b}", f"{b}+{a}"]:
        if key in scores:
            return scores[key]
    return 25.0  # default midpoint

def confusion_matrix_stats(model_pos, faers_pos):
    tp = sum(1 for m, f in zip(model_pos, faers_pos) if m and f)
    fp = sum(1 for m, f in zip(model_pos, faers_pos) if m and not f)
    fn = sum(1 for m, f in zip(model_pos, faers_pos) if not m and f)
    tn = sum(1 for m, f in zip(model_pos, faers_pos) if not m and not f)
    n  = tp + fp + fn + tn

    sens = tp / (tp + fn) if (tp + fn) > 0 else 0
    spec = tn / (tn + fp) if (tn + fp) > 0 else 0
    ppv  = tp / (tp + fp) if (tp + fp) > 0 else 0
    npv  = tn / (tn + fn) if (tn + fn) > 0 else 0

    # Cohen's kappa
    po = (tp + tn) / n
    pe = ((tp+fp)*(tp+fn) + (tn+fn)*(tn+fp)) / (n*n)
    kappa = (po - pe) / (1 - pe) if (1-pe) > 0 else 0

    # Fisher's exact
    table = [[tp, fp], [fn, tn]]
    OR, pval = stats.fisher_exact(table)

    return {
        "tp":tp, "fp":fp, "fn":fn, "tn":tn, "n":n,
        "sensitivity":round(sens,3), "specificity":round(spec,3),
        "PPV":round(ppv,3), "NPV":round(npv,3),
        "kappa":round(kappa,3), "fisher_OR":round(OR,2),
        "fisher_p":round(pval,4),
    }

def print_stats(label, stats_dict):
    d = stats_dict
    print(f"\n── {label} (n={d['n']}) ──")
    print(f"  TP={d['tp']}  FP={d['fp']}  FN={d['fn']}  TN={d['tn']}")
    print(f"  Sensitivity: {d['sensitivity']:.3f}   Specificity: {d['specificity']:.3f}")
    print(f"  PPV:         {d['PPV']:.3f}   NPV:         {d['NPV']:.3f}")
    print(f"  Kappa:       {d['kappa']:.3f}")
    print(f"  Fisher:      OR={d['fisher_OR']:.2f}  p={d['fisher_p']:.4f}")

def main():
    df = pd.DataFrame(ALL_PAIRS,
        columns=["drug_a","drug_b","dQTc","tier","faers_signal","ROR","gap"])

    df["model_pos"] = df["tier"].isin(["MODERATE","HIGH"])

    # Add composite scores
    df["composite"] = df.apply(
        lambda r: get_composite(r.drug_a, r.drug_b, COMPOSITE_SCORES), axis=1)

    print("=" * 65)
    print("CONCORDANCE STATISTICS — CardioSafe Pediatric")
    print("=" * 65)

    # ── 1. Overall ────────────────────────────────────────────────────────────
    s_all = confusion_matrix_stats(df["model_pos"], df["faers_signal"])
    print_stats("OVERALL (n=63)", s_all)

    # ── 2. Remove gap pairs ───────────────────────────────────────────────────
    df_nogap = df[df["gap"].isna()].copy()
    s_nogap = confusion_matrix_stats(df_nogap["model_pos"], df_nogap["faers_signal"])
    print_stats("EXCLUDING KNOWN GAPS (conduction + CYP2D6)", s_nogap)

    # ── 3. Stimulant only ─────────────────────────────────────────────────────
    df_stim = df[df.apply(lambda r: r.drug_a in ("MPH","AMP") or
                          r.drug_b in ("MPH","AMP"), axis=1)].copy()
    s_stim = confusion_matrix_stats(df_stim["model_pos"], df_stim["faers_signal"])
    print_stats("STIMULANT-CONTAINING ONLY", s_stim)

    # ── 4. Non-stimulant, non-gap ─────────────────────────────────────────────
    df_nonstim_nogap = df[
        df["gap"].isna() &
        ~df.apply(lambda r: r.drug_a in ("MPH","AMP") or
                  r.drug_b in ("MPH","AMP"), axis=1)
    ].copy()
    s_nonstim = confusion_matrix_stats(
        df_nonstim_nogap["model_pos"], df_nonstim_nogap["faers_signal"])
    print_stats("NON-STIMULANT, NON-GAP PAIRS", s_nonstim)

    # ── 5. ROC curve using composite score ───────────────────────────────────
    print("\n── ROC CURVE (composite score vs FAERS signal) ──")
    fpr, tpr, thresholds = roc_curve(
        df["faers_signal"].astype(int),
        df["composite"]
    )
    roc_auc = auc(fpr, tpr)
    print(f"  AUC-ROC: {roc_auc:.3f}")

    # Find optimal threshold (Youden's J)
    youden = tpr - fpr
    opt_idx = np.argmax(youden)
    opt_thresh = thresholds[opt_idx]
    print(f"  Optimal threshold (Youden): {opt_thresh:.1f}")
    print(f"  At optimal: TPR={tpr[opt_idx]:.3f}  FPR={fpr[opt_idx]:.3f}")

    # ROC for no-gap subset
    fpr2, tpr2, _ = roc_curve(
        df_nogap["faers_signal"].astype(int),
        df_nogap["composite"]
    )
    roc_auc2 = auc(fpr2, tpr2)
    print(f"  AUC-ROC (no-gap subset): {roc_auc2:.3f}")

    # ── PLOT ──────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # ROC curves
    ax = axes[0]
    ax.plot(fpr, tpr, color="#1B3A6B", lw=2,
            label=f"All pairs (AUC = {roc_auc:.2f})")
    ax.plot(fpr2, tpr2, color="#0D6B7A", lw=2, linestyle="--",
            label=f"Excl. known gaps (AUC = {roc_auc2:.2f})")
    ax.plot([0,1],[0,1], color="#CDD3DE", lw=1, linestyle=":")
    ax.scatter(fpr[opt_idx], tpr[opt_idx], color="#B01C1C", s=80, zorder=5,
               label=f"Optimal threshold ({opt_thresh:.0f})")
    ax.set_xlabel("False Positive Rate", fontsize=11)
    ax.set_ylabel("True Positive Rate", fontsize=11)
    ax.set_title("ROC Curve — Composite Score vs FAERS Signal", fontsize=11,
                 fontweight="bold", color="#1B3A6B")
    ax.legend(fontsize=9, framealpha=0.9)
    ax.set_xlim(0,1); ax.set_ylim(0,1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Sensitivity comparison bar chart
    ax2 = axes[1]
    subsets = ["Overall", "No gaps", "Stimulant\nonly", "Non-stim\nnon-gap"]
    sensitivities = [s_all["sensitivity"], s_nogap["sensitivity"],
                     s_stim["sensitivity"], s_nonstim["sensitivity"]]
    ns = [s_all["n"], s_nogap["n"], s_stim["n"], s_nonstim["n"]]
    colors = ["#6B7280","#1B3A6B","#D97706","#0D6B7A"]
    bars = ax2.bar(range(4), sensitivities, color=colors, width=0.55,
                   edgecolor="white")
    for i, (v, n) in enumerate(zip(sensitivities, ns)):
        ax2.text(i, v+0.02, f"{v:.2f}\n(n={n})", ha="center",
                 fontsize=9, color="#1F2937", fontweight="bold")
    ax2.axhline(0.5, color="#CDD3DE", linewidth=1, linestyle="--",
                label="0.5 reference")
    ax2.set_xticks(range(4))
    ax2.set_xticklabels(subsets, fontsize=10)
    ax2.set_ylabel("Sensitivity", fontsize=11)
    ax2.set_ylim(0, 1.0)
    ax2.set_title("Sensitivity by Analysis Subset", fontsize=11,
                  fontweight="bold", color="#1B3A6B")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.yaxis.grid(True, color="#CDD3DE", linewidth=0.5, linestyle="--")
    ax2.set_axisbelow(True)

    plt.tight_layout()
    plt.savefig("/mnt/user-data/outputs/figure_roc_concordance.png",
                dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print("\nSaved -> figure_roc_concordance.png")

    # ── SUMMARY TABLE ────────────────────────────────────────────────────────
    print("\n" + "=" * 65)
    print("SUMMARY FOR MANUSCRIPT")
    print("=" * 65)
    print(f"  Overall AUC-ROC:              {roc_auc:.3f}")
    print(f"  No-gap AUC-ROC:               {roc_auc2:.3f}")
    print(f"  Overall sensitivity:           {s_all['sensitivity']:.3f}")
    print(f"  No-gap sensitivity:            {s_nogap['sensitivity']:.3f}")
    print(f"  Stimulant sensitivity:         {s_stim['sensitivity']:.3f}")
    print(f"  Non-stim no-gap sensitivity:   {s_nonstim['sensitivity']:.3f}")
    print(f"  Improvement excl gaps:         "
          f"+{s_nogap['sensitivity']-s_all['sensitivity']:.3f}")

if __name__ == "__main__":
    main()
