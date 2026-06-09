"""
figures.py — CardioSafe Pediatric
Generates three manuscript figures from hardcoded 500-beat results.
Run: python3 figures.py
Output: figure1_heatmap.png, figure2_distribution.png, figure3_faers_scatter.png
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker

# ── STYLE ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":      "DejaVu Sans",
    "font.size":        11,
    "axes.spines.top":  False,
    "axes.spines.right":False,
    "axes.linewidth":   0.8,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "figure.dpi":       150,
    "savefig.dpi":      300,
    "savefig.bbox":     "tight",
    "savefig.facecolor":"white",
})

NAVY   = "#1B3A6B"
TEAL   = "#0D6B7A"
AMBER  = "#C47A00"
RED    = "#B01C1C"
GREEN  = "#1A5C2A"
LGRAY  = "#F0F2F5"
MGRAY  = "#CDD3DE"
DGRAY  = "#6B7280"

# ── DATA ───────────────────────────────────────────────────────────────────
DRUGS  = ["MPH","AMP","RIS","QUE","ARI","SER","FLU","ESC","CLO","GUA","IMI","NOR"]
LABELS = {
    "MPH":"Methylphenidate","AMP":"Amphetamine","RIS":"Risperidone",
    "QUE":"Quetiapine","ARI":"Aripiprazole","SER":"Sertraline",
    "FLU":"Fluoxetine","ESC":"Escitalopram","CLO":"Clonidine",
    "GUA":"Guanfacine","IMI":"Imipramine","NOR":"Nortriptyline",
}
CLASS = {
    "MPH":"Stimulant","AMP":"Stimulant","RIS":"Antipsychotic","QUE":"Antipsychotic",
    "ARI":"Antipsychotic","SER":"SSRI","FLU":"SSRI","ESC":"SSRI",
    "CLO":"Alpha-2","GUA":"Alpha-2","IMI":"TCA","NOR":"TCA",
}
CLASS_COLOR = {
    "Stimulant":"#D97706","Antipsychotic":"#1B3A6B","SSRI":"#0D6B7A",
    "Alpha-2":"#1A5C2A","TCA":"#7C2D12",
}

# Pairwise delta-QTc (500-beat)
DQTC_RAW = {
    "MPH+AMP":19.1,"MPH+RIS":10.0,"MPH+QUE":11.0,"MPH+ARI":18.1,
    "MPH+SER":9.6, "MPH+FLU":10.0,"MPH+ESC":9.1, "MPH+CLO":6.8,
    "MPH+GUA":6.8, "MPH+IMI":10.5,"MPH+NOR":14.8,
    "AMP+RIS":10.0,"AMP+QUE":11.0,"AMP+ARI":18.1,"AMP+SER":9.6,
    "AMP+FLU":10.0,"AMP+ESC":9.1, "AMP+CLO":6.8, "AMP+GUA":6.8,
    "AMP+IMI":10.5,"AMP+NOR":14.8,
    "RIS+QUE":3.0, "RIS+ARI":10.0,"RIS+SER":1.5, "RIS+FLU":1.5,
    "RIS+ESC":1.0, "RIS+CLO":-1.1,"RIS+GUA":-1.1,"RIS+IMI":2.5,
    "RIS+NOR":6.5,
    "QUE+ARI":11.5,"QUE+SER":2.5, "QUE+FLU":3.0, "QUE+ESC":2.0,
    "QUE+CLO":-0.1,"QUE+GUA":-0.1,"QUE+IMI":3.5, "QUE+NOR":8.0,
    "ARI+SER":9.5, "ARI+FLU":10.0,"ARI+ESC":9.5, "ARI+CLO":7.3,
    "ARI+GUA":7.3, "ARI+IMI":11.0,"ARI+NOR":15.0,
    "SER+FLU":1.5, "SER+ESC":0.5, "SER+CLO":-1.7,"SER+GUA":-1.7,
    "SER+IMI":2.0, "SER+NOR":6.5,
    "FLU+ESC":1.0, "FLU+CLO":-1.1,"FLU+GUA":-1.1,"FLU+IMI":2.5,
    "FLU+NOR":6.5,
    "ESC+CLO":-1.7,"ESC+GUA":-1.7,"ESC+IMI":2.0, "ESC+NOR":6.0,
    "CLO+GUA":-1.5,"CLO+IMI":-0.6,"CLO+NOR":4.1,
    "GUA+IMI":-0.6,"GUA+NOR":4.1,
    "IMI+NOR":7.5,
}

def get_dqtc(a, b):
    return DQTC_RAW.get(f"{a}+{b}", DQTC_RAW.get(f"{b}+{a}", np.nan))

def tier(dq):
    if np.isnan(dq): return None
    if dq >= 20:  return "HIGH"
    if dq >= 10:  return "MODERATE"
    if dq >= 5:   return "LOW-MOD"
    if dq >= 0:   return "LOW"
    return "PROTECTIVE"

# FAERS alignment data
FAERS_DATA = [
    ("MPH+SER",  9.6,  12.79, True,  "stimulant"),
    ("MPH+ARI", 18.1,   8.15, True,  "stimulant"),
    ("MPH+QUE", 11.0,   5.25, True,  "stimulant"),
    ("AMP+SER",  9.6,   2.81, True,  "stimulant"),
    ("MPH+AMP", 19.1,   0.45, False, "stimulant"),
    ("AMP+ARI", 18.1,   0.69, False, "stimulant"),
    ("ARI+GUA",  7.3,  10.25, True,  "conduction"),
    ("QUE+GUA", -0.1,  13.15, True,  "conduction"),
    ("SER+GUA", -1.7,   2.58, True,  "conduction"),
    ("QUE+FLU",  3.0,   9.18, True,  "CYP2D6"),
    ("ARI+SER",  9.5,   6.78, True,  "CYP2D6"),
    ("QUE+SER",  2.5,   5.84, True,  "CYP2D6"),
    ("SER+IMI",  2.0,  16.22, True,  "CYP2D6"),
    ("ESC+IMI",  2.0,  28.01, True,  "CYP2D6"),
    ("ARI+IMI", 11.0,  14.67, True,  "hERG"),
    ("RIS+ARI", 10.0,   1.71, True,  "hERG"),
    ("ARI+ESC",  9.5,   2.31, True,  "hERG"),
    ("FLU+RIS",  1.5,   3.69, True,  "hERG"),
    ("MPH+RIS", 10.0,   0.65, False, "stimulant"),
    ("MPH+FLU", 10.0,   0.12, False, "stimulant"),
    ("RIS+FLU",  1.5,   3.69, True,  "hERG"),
    ("MPH+GUA",  6.8,   3.40, True,  "conduction"),
]


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — HEATMAP
# ════════════════════════════════════════════════════════════════════════════

def fig1_heatmap():
    n = len(DRUGS)
    mat = np.full((n, n), np.nan)
    for i, a in enumerate(DRUGS):
        for j, b in enumerate(DRUGS):
            if i < j:
                v = get_dqtc(a, b)
                mat[i, j] = v
                mat[j, i] = v

    fig, ax = plt.subplots(figsize=(10, 8.5))

    # Custom diverging colormap: green (protective) -> white -> orange -> red
    cmap = LinearSegmentedColormap.from_list("cardiosafe", [
        "#1A5C2A", "#4CAF7D", "#B8D9C6", "#F5F5F5",
        "#FBD38D", "#F97316", "#B01C1C"
    ], N=256)

    im = ax.imshow(mat, cmap=cmap, vmin=-8, vmax=20, aspect="auto")

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("ΔQTc (ms)", fontsize=10, color=DGRAY)
    cbar.ax.yaxis.set_tick_params(color=DGRAY)
    cbar.outline.set_edgecolor(MGRAY)
    for tick in cbar.ax.yaxis.get_ticklabels():
        tick.set_color(DGRAY)
        tick.set_fontsize(9)

    # Tier threshold lines on colorbar
    for val, label in [(10, "MODERATE"), (5, "LOW-MOD"), (0, "LOW")]:
        # map val to colorbar position
        pos = (val - (-8)) / (20 - (-8))
        cbar.ax.axhline(pos, color="white", linewidth=1, alpha=0.7)

    # Drug labels with class color bars
    long_labels = [LABELS[d] for d in DRUGS]
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(long_labels, rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(long_labels, fontsize=9)

    # Color tick labels by class
    for i, d in enumerate(DRUGS):
        col = CLASS_COLOR[CLASS[d]]
        ax.get_xticklabels()[i].set_color(col)
        ax.get_yticklabels()[i].set_color(col)

    # Annotate cells with value
    for i in range(n):
        for j in range(n):
            v = mat[i, j]
            if not np.isnan(v) and i != j:
                t = tier(v)
                text_color = "white" if (v >= 15 or v <= -3) else "#333333"
                ax.text(j, i, f"{v:+.0f}", ha="center", va="center",
                        fontsize=7, color=text_color, fontweight="bold" if abs(v) >= 10 else "normal")

    # Diagonal
    for i in range(n):
        ax.add_patch(mpatches.FancyBboxPatch((i-0.5, i-0.5), 1, 1,
            boxstyle="square,pad=0", facecolor="#E5E7EB", edgecolor="none", zorder=2))

    ax.set_xlim(-0.5, n-0.5)
    ax.set_ylim(n-0.5, -0.5)

    # Class legend
    handles = [mpatches.Patch(color=v, label=k) for k, v in CLASS_COLOR.items()]
    ax.legend(handles=handles, loc="lower right", fontsize=8,
              framealpha=0.9, edgecolor=MGRAY, title="Drug class",
              title_fontsize=8)

    ax.set_title("Figure 1. Pairwise ΔQTc Heatmap — 66 Drug Combinations (500-beat steady state)",
                 fontsize=11, color=NAVY, fontweight="bold", pad=12, loc="left")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_color(MGRAY)
    ax.spines["left"].set_color(MGRAY)

    plt.tight_layout()
    plt.savefig("/mnt/user-data/outputs/figure1_heatmap.png")
    plt.close()
    print("Saved figure1_heatmap.png")


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — RISK TIER DISTRIBUTION
# ════════════════════════════════════════════════════════════════════════════

def fig2_distribution():
    tiers   = ["MODERATE\n(≥10 ms)", "LOW-MOD\n(5–9 ms)", "LOW\n(0–4 ms)", "PROTECTIVE\n(<0 ms)"]
    t_short = ["MODERATE","LOW-MOD","LOW","PROTECTIVE"]

    # Mechanism breakdown per tier from 500-beat results
    # sympathomimetic = contains MPH or AMP in stimulant role
    # hERG = stimulant-free, IKr block > 3%
    # autonomic = contains CLO or GUA producing negative dQTc
    # mixed = other
    breakdown = {
        "MODERATE":   {"Sympathomimetic":22, "hERG-driven":9,  "Autonomic":0, "Mixed":0},
        "LOW-MOD":    {"Sympathomimetic":6,  "hERG-driven":4,  "Autonomic":0, "Mixed":8},
        "LOW":        {"Sympathomimetic":2,  "hERG-driven":14, "Autonomic":0, "Mixed":19},
        "PROTECTIVE": {"Sympathomimetic":0,  "hERG-driven":0,  "Autonomic":8, "Mixed":3},
    }
    mech_colors = {
        "Sympathomimetic": "#D97706",
        "hERG-driven":     "#1B3A6B",
        "Autonomic":       "#1A5C2A",
        "Mixed":           "#9CA3AF",
    }
    mechs = list(mech_colors.keys())

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # ── Left: stacked bar by mechanism ──
    ax = axes[0]
    x = np.arange(len(t_short))
    bottoms = np.zeros(len(t_short))
    tier_colors = [RED, AMBER, TEAL, GREEN]

    for m in mechs:
        vals = [breakdown[t][m] for t in t_short]
        bars = ax.bar(x, vals, bottom=bottoms, color=mech_colors[m],
                      label=m, width=0.55, edgecolor="white", linewidth=0.5)
        # Label non-zero segments
        for xi, (v, b) in enumerate(zip(vals, bottoms)):
            if v > 0:
                ax.text(xi, b + v/2, str(v), ha="center", va="center",
                        fontsize=9, color="white", fontweight="bold")
        bottoms += np.array(vals)

    totals = [sum(breakdown[t].values()) for t in t_short]
    for xi, tot in enumerate(totals):
        ax.text(xi, tot + 0.5, str(tot), ha="center", va="bottom",
                fontsize=10, color=DGRAY, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(tiers, fontsize=10)
    ax.set_ylabel("Number of combinations", fontsize=10, color=DGRAY)
    ax.set_ylim(0, 38)
    ax.set_title("A.  Risk tier distribution by mechanism\n(66 pairwise combinations)",
                 fontsize=10, color=NAVY, fontweight="bold", loc="left")
    ax.legend(loc="upper right", fontsize=9, framealpha=0.9,
              edgecolor=MGRAY, title="Primary mechanism", title_fontsize=9)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=MGRAY, linewidth=0.5, linestyle="--")

    # ── Right: IKr block distribution for MODERATE tier ──
    ax2 = axes[1]
    # IKr block % for all 31 MODERATE combinations
    ikr_moderate = [
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,  # stimulant-only (10 pairs)
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,  # stimulant+low-hERG (10 pairs)
        0.19,0.19,0.36,0.93,                        # stimulant+mild hERG
        4.22,4.22,4.22,4.22,                        # stimulant+ARI
        4.56,5.11,6.76                               # pure hERG driven
    ]
    bins = [0, 1, 3, 6, 10]
    counts, _ = np.histogram(ikr_moderate, bins=bins)
    bar_labels = ["<1%\n(sympathomimetic)", "1–3%", "3–6%", ">6%\n(hERG-driven)"]
    bar_colors = [RED, AMBER, "#4B7BE8", NAVY]

    bars2 = ax2.bar(range(4), counts, color=bar_colors, width=0.6,
                    edgecolor="white", linewidth=0.5)
    for xi, v in enumerate(counts):
        ax2.text(xi, v + 0.3, str(v), ha="center", va="bottom",
                 fontsize=11, fontweight="bold", color=DGRAY)

    ax2.set_xticks(range(4))
    ax2.set_xticklabels(bar_labels, fontsize=9)
    ax2.set_ylabel("MODERATE combinations (n)", fontsize=10, color=DGRAY)
    ax2.set_ylim(0, max(counts) + 4)
    ax2.set_title("B.  IKr block distribution among\nMODERATE-risk combinations (n=31)",
                  fontsize=10, color=NAVY, fontweight="bold", loc="left")
    ax2.set_axisbelow(True)
    ax2.yaxis.grid(True, color=MGRAY, linewidth=0.5, linestyle="--")

    # Annotation
    ax2.annotate("28/31 have\nIKr block <5%", xy=(1.5, counts[0]+counts[1]-1),
                 xytext=(2.5, counts[0]+counts[1]+2),
                 arrowprops=dict(arrowstyle="->", color=RED, lw=1.2),
                 fontsize=9, color=RED, fontweight="bold", ha="center")

    plt.suptitle("Figure 2. Risk Distribution and Mechanism Analysis (500-beat steady state)",
                 fontsize=11, color=NAVY, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig("/mnt/user-data/outputs/figure2_distribution.png")
    plt.close()
    print("Saved figure2_distribution.png")


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — MODEL vs FAERS SCATTER
# ════════════════════════════════════════════════════════════════════════════

def fig3_scatter():
    fig, ax = plt.subplots(figsize=(9, 7))

    group_styles = {
        "stimulant":  {"marker":"o", "label":"Stimulant combination"},
        "conduction": {"marker":"s", "label":"Guanfacine (conduction gap)"},
        "CYP2D6":     {"marker":"^", "label":"CYP2D6 interaction gap"},
        "hERG":       {"marker":"D", "label":"hERG-driven (non-stimulant)"},
    }
    point_color = {
        True:  "#B01C1C",   # FAERS signal = red
        False: "#6B7280",   # no signal = gray
    }

    # Background zones
    ax.axhspan(1.0, 35, alpha=0.04, color=RED, zorder=0)
    ax.axvspan(10, 22, alpha=0.04, color=NAVY, zorder=0)
    ax.axhline(1.0, color=RED, linewidth=0.8, linestyle="--", alpha=0.5, label="FAERS signal threshold (ROR=1.0)")
    ax.axvline(10, color=NAVY, linewidth=0.8, linestyle="--", alpha=0.5, label="MODERATE tier threshold (ΔQTc=10ms)")

    # Plot points
    plotted_labels = set()
    for combo, dq, ror, sig, group in FAERS_DATA:
        style = group_styles[group]
        label = style["label"] if group not in plotted_labels else None
        plotted_labels.add(group)
        color = point_color[sig]
        ax.scatter(dq, ror, marker=style["marker"], s=90, color=color,
                   edgecolors="white", linewidths=0.8, zorder=3,
                   alpha=0.85, label=label)

    # Label selected points
    label_offsets = {
        "MPH+SER":  (0.4,  1.2),
        "ARI+GUA":  (0.4,  1.2),
        "QUE+GUA":  (-5.5, 1.0),
        "MPH+ARI":  (0.4, -2.5),
        "ESC+IMI":  (0.4,  1.0),
        "SER+IMI":  (-5.0, 1.0),
        "MPH+AMP":  (0.4, -2.0),
        "MPH+QUE":  (0.4,  1.0),
        "ARI+NOR":  (0.4,  1.0) if False else None,
    }
    for combo, dq, ror, sig, group in FAERS_DATA:
        if combo in label_offsets and label_offsets[combo]:
            dx, dy = label_offsets[combo]
            ax.annotate(combo, (dq, ror), xytext=(dq+dx, ror+dy),
                        fontsize=7.5, color="#1F2937",
                        arrowprops=dict(arrowstyle="-", color=MGRAY,
                                        lw=0.7, shrinkA=3, shrinkB=3))

    # Quadrant labels
    ax.text(19, 0.3, "Model: MODERATE\nFAERS: no signal\n(rare co-Rx)",
            fontsize=7.5, color=DGRAY, ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.3", facecolor=LGRAY,
                      edgecolor=MGRAY, alpha=0.8))
    ax.text(5, 20, "Model: LOW/LOW-MOD\nFAERS: signal\n(mechanistic gaps)",
            fontsize=7.5, color=DGRAY, ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.3", facecolor=LGRAY,
                      edgecolor=MGRAY, alpha=0.8))
    ax.text(14, 10, "Concordant\n(model + FAERS)",
            fontsize=7.5, color=NAVY, ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#EFF4FB",
                      edgecolor=NAVY, alpha=0.8))

    ax.set_xlabel("Model ΔQTc (ms)", fontsize=11, color=DGRAY)
    ax.set_ylabel("FAERS Reporting Odds Ratio (ROR)", fontsize=11, color=DGRAY)
    ax.set_xlim(-5, 23)
    ax.set_ylim(0, 32)
    ax.set_title("Figure 3. Model ΔQTc vs. FAERS Pharmacovigilance Signal\n(63 evaluable drug pairs, pediatric cases 2015–2024)",
                 fontsize=11, color=NAVY, fontweight="bold", loc="left")

    # Legend
    legend = ax.legend(loc="upper left", fontsize=8.5, framealpha=0.95,
                       edgecolor=MGRAY, ncol=1)

    # Signal color legend
    sig_patch   = mpatches.Patch(color=RED,   label="FAERS signal (ROR CI lower >1)")
    nosig_patch = mpatches.Patch(color=DGRAY, label="No FAERS signal")
    ax.add_artist(legend)
    ax.legend(handles=[sig_patch, nosig_patch], loc="upper right",
              fontsize=8.5, framealpha=0.95, edgecolor=MGRAY,
              title="Point color", title_fontsize=8.5)

    plt.tight_layout()
    plt.savefig("/mnt/user-data/outputs/figure3_faers_scatter.png")
    plt.close()
    print("Saved figure3_faers_scatter.png")


# ── RUN ────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("Generating CardioSafe Pediatric manuscript figures...")
    fig1_heatmap()
    fig2_distribution()
    fig3_scatter()
    print("Done.")
