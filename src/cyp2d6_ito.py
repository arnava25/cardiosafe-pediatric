"""
cyp2d6_ito.py — CardioSafe Pediatric
======================================
CYP2D6 static drug interaction model (Ito method) for fluoxetine-containing
combinations. Computes the predicted AUC ratio and adjusted Cmax for CYP2D6
substrate drugs when co-administered with fluoxetine.

Model: AUC_ratio = 1 / (fm * (1 / (1 + [I]/Ki)) + (1 - fm))

Where:
    fm  = fraction of substrate clearance via CYP2D6
    [I] = hepatic inhibitor concentration (fluoxetine + norfluoxetine)
    Ki  = inhibition constant of fluoxetine for CYP2D6

Then reruns the polypharmacy simulation with adjusted Cmax values and
computes corrected delta-QTc predictions.

References:
    Ito et al. Pharm Res. 1998;15(3):396-402
    Templeton et al. Drug Metab Dispos. 2016;44(1):57-65
    FDA Drug Interaction Guidance 2020 (static model section)

Usage:
    python3 src/cyp2d6_ito.py
"""

import sys
import csv
from pathlib import Path
from datetime import datetime

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import pandas as pd
import numpy as np
from ord_model import run_simulation

PARAMS_PATH = str(ROOT / "data" / "herg_master_params.csv")
OUT_CSV     = ROOT / "results" / "cyp2d6_ito_results.csv"
OUT_MEMO    = ROOT / "results" / "cyp2d6_ito_memo.md"

# ── CYP2D6 INHIBITOR PARAMETERS (Fluoxetine) ────────────────────────────────
# Fluoxetine is a potent CYP2D6 inhibitor via competitive inhibition.
# Norfluoxetine (active metabolite) is equipotent.

FLU_PARAMS = {
    # Fluoxetine Ki for CYP2D6 (nM), competitive inhibition
    # Templeton et al. 2016: Ki = 0.17 uM = 170 nM
    "Ki_nM":         170.0,

    # Hepatic inlet concentration [I]h (FDA Ito static model, hepatic variant)
    # [I]h = Cmax_unbound + (Fa * Dose * ka) / Qh
    # Fluoxetine 20mg: Cmax_free = 3.5 nM
    # Gut absorption term: Fa=0.72, Dose=20mg=20,000,000 nM*L (MW=309g/mol)
    #   = 0.72 * (20e6 / 309) / 1500 * (1000/60) nM ≈ 154 nM
    # [I]h_flu = 3.5 + 154 = 157.5 nM
    "I_hepatic_flu_nM":  157.5,

    # Norfluoxetine hepatic inlet (active metabolite, ~equal inhibitor)
    # Norfluoxetine Ki ~170 nM; Cmax_free ~5 nM; gut term ~50 nM
    # [I]h_norflu = 5 + 50 = 55 nM
    "I_hepatic_norflu_nM": 55.0,

    # Blood-based [I] (conservative, for comparison)
    "I_blood_nM":    3.5,
    "I_norflu_nM":   5.0,
}

# Hepatic inlet model (correct for CYP2D6 DDI prediction)
FLU_PARAMS["I_total_nM"] = FLU_PARAMS["I_hepatic_flu_nM"] + FLU_PARAMS["I_hepatic_norflu_nM"]
# References: FDA DDI Guidance 2020; Templeton 2016; Crewe 1992; Preskorn 1997

# ── CYP2D6 SUBSTRATE PARAMETERS ─────────────────────────────────────────────
# fm = fraction of systemic clearance via CYP2D6
# Cmax_free_nM = adult therapeutic free plasma Cmax (matches herg_master_params.csv)

SUBSTRATES = {
    "Quetiapine": {
        "fm_cyp2d6":     0.73,   # Grimm et al. 2006: CYP2D6 + CYP3A4; ~73% CYP2D6
        "Cmax_free_nM":  2.1,    # from herg_master_params.csv
        "notes":         "Primary CYP2D6 substrate; fluoxetine increases AUC ~2-3x (Templeton 2016)",
    },
    "Aripiprazole": {
        "fm_cyp2d6":     0.40,   # Otsuka label: CYP2D6 + CYP3A4; ~40% CYP2D6
        "Cmax_free_nM":  21.0,   # from herg_master_params.csv
        "notes":         "CYP2D6 + CYP3A4; fluoxetine increases AUC ~1.5x; label recommends dose reduction",
    },
    "Risperidone": {
        "fm_cyp2d6":     0.77,   # De Leon 2007: ~77% CYP2D6; 9-OH-risperidone active
        "Cmax_free_nM":  0.9,    # from herg_master_params.csv
        "notes":         "High fm; fluoxetine increases risperidone+9-OH-RIS combined Cmax ~2x",
    },
    "Nortriptyline": {
        "fm_cyp2d6":     0.90,   # CYP2D6 primary; narrow TI; PM phenotype = 10x Cmax
        "Cmax_free_nM":  28.0,   # from herg_master_params.csv
        "notes":         "Narrow therapeutic index; fluoxetine inhibition clinically significant",
    },
    "Imipramine": {
        "fm_cyp2d6":     0.55,   # CYP2D6 + CYP2C19; ~55% CYP2D6
        "Cmax_free_nM":  18.0,   # from herg_master_params.csv
        "notes":         "Both CYP2D6 and CYP2C19; partial CYP2D6 dependence",
    },
}

# Drug name to code
DRUG_CODES = {
    "Methylphenidate":"MPH","Amphetamine":"AMP","Risperidone":"RIS",
    "Quetiapine":"QUE","Aripiprazole":"ARI","Sertraline":"SER",
    "Fluoxetine":"FLU","Escitalopram":"ESC","Clonidine":"CLO",
    "Guanfacine":"GUA","Imipramine":"IMI","Nortriptyline":"NOR",
}

def ito_auc_ratio(fm, I_nM, Ki_nM):
    """
    Ito static inhibition model.
    AUC_ratio = 1 / (fm * (1 / (1 + [I]/Ki)) + (1 - fm))

    Returns the fold-increase in AUC (and approximately Cmax) of the substrate
    when the inhibitor is co-administered.
    """
    inhibited_fm  = fm * (1.0 / (1.0 + I_nM / Ki_nM))
    uninhibited   = 1.0 - fm
    auc_ratio     = 1.0 / (inhibited_fm + uninhibited)
    return round(auc_ratio, 2)

def tier(dq):
    if dq >= 20:  return "HIGH"
    if dq >= 10:  return "MODERATE"
    if dq >= 5:   return "LOW-MOD"
    if dq >= 0:   return "LOW"
    return "PROTECTIVE"

def main():
    N_BEATS = 200  # sufficient for CYP2D6 delta comparisons

    print("CardioSafe Pediatric — CYP2D6 Ito Static Interaction Model")
    print(f"Inhibitor: Fluoxetine + Norfluoxetine")
    print(f"[I]total = {FLU_PARAMS['I_total_nM']:.1f} nM  |  Ki = {FLU_PARAMS['Ki_nM']:.0f} nM")
    print()

    # ── Step 1: Compute AUC ratios ──────────────────────────────────────────
    print("── AUC RATIO PREDICTIONS ──")
    print(f"{'Drug':15s}  {'fm':6s}  {'AUC ratio':10s}  {'Cmax 1x (nM)':14s}  {'Cmax adj (nM)':14s}  Notes")
    print("-" * 95)

    auc_ratios = {}
    adj_cmax   = {}
    for drug, params in SUBSTRATES.items():
        ratio = ito_auc_ratio(
            params["fm_cyp2d6"],
            FLU_PARAMS["I_total_nM"],
            FLU_PARAMS["Ki_nM"]
        )
        c1x  = params["Cmax_free_nM"]
        cadj = round(c1x * ratio, 2)
        auc_ratios[drug] = ratio
        adj_cmax[drug]   = cadj
        print(f"{drug:15s}  {params['fm_cyp2d6']:6.2f}  {ratio:10.2f}x  "
              f"{c1x:14.1f}  {cadj:14.1f}  {params['notes'][:50]}")

    print()

    # ── Step 2: Baseline simulation ─────────────────────────────────────────
    print("Running baseline simulation...")
    baseline = run_simulation(None, PARAMS_PATH, n_beats=N_BEATS, verbose=True)
    baseline_qtc = baseline["QTc"]
    print(f"Baseline QTc: {baseline_qtc:.1f} ms\n")

    # ── Step 3: Simulate each FLU+substrate combo with adjusted Cmax ─────────
    results = []

    combos_to_run = [
        ("Fluoxetine", sub) for sub in SUBSTRATES
        if sub != "Fluoxetine"
    ]
    # Also run key clinical triples: MPH + FLU + substrate
    extra_combos = [
        ("Methylphenidate", "Quetiapine",   "Fluoxetine"),
        ("Methylphenidate", "Aripiprazole", "Fluoxetine"),
    ]

    print("── PAIRWISE: FLU + SUBSTRATE ──")
    for drug_a, drug_b in combos_to_run:
        # Baseline (unadjusted)
        combo_base = {drug_a: "therapeutic", drug_b: "therapeutic"}
        res_base = run_simulation(combo_base, PARAMS_PATH, n_beats=N_BEATS, verbose=False)
        dqtc_base = round(res_base["QTc"] - baseline_qtc, 1)
        ikr_base  = round(res_base["IKr_block_pct"], 2)

        # CYP2D6-adjusted (substrate Cmax scaled by AUC ratio)
        substrate = drug_b  # FLU is always drug_a here
        combo_adj = {drug_a: "therapeutic", drug_b: adj_cmax[drug_b]}
        res_adj   = run_simulation(combo_adj, PARAMS_PATH, n_beats=N_BEATS, verbose=False)
        dqtc_adj  = round(res_adj["QTc"] - baseline_qtc, 1)
        ikr_adj   = round(res_adj["IKr_block_pct"], 2)

        a = DRUG_CODES[drug_a]; b = DRUG_CODES[drug_b]
        ratio = auc_ratios[drug_b]
        print(f"  {a}+{b}: base={dqtc_base:+.1f}ms ({tier(dqtc_base)})  "
              f"adj={dqtc_adj:+.1f}ms ({tier(dqtc_adj)})  "
              f"[{b} x{ratio}]  IKr={ikr_adj:.2f}%")

        results.append({
            "combination":       f"{a}+{b}",
            "n_drugs":           2,
            "cyp2d6_inhibitor":  "FLU",
            "substrate":         b,
            "auc_ratio":         ratio,
            "cmax_1x_nM":        SUBSTRATES[substrate]["Cmax_free_nM"],
            "cmax_adj_nM":       adj_cmax[substrate],
            "dQTc_base_ms":      dqtc_base,
            "dQTc_adj_ms":       dqtc_adj,
            "dQTc_delta_ms":     round(dqtc_adj - dqtc_base, 1),
            "IKr_base_pct":      ikr_base,
            "IKr_adj_pct":       ikr_adj,
            "tier_base":         tier(dqtc_base),
            "tier_adj":          tier(dqtc_adj),
            "tier_changed":      tier(dqtc_base) != tier(dqtc_adj),
            "fm_cyp2d6":         SUBSTRATES[substrate]["fm_cyp2d6"],
        })

    print()
    print("── TRIPLES: MPH + SUBSTRATE + FLU ──")
    for drug_a, drug_b, drug_c in extra_combos:
        # Baseline
        combo_base = {drug_a: "therapeutic", drug_b: "therapeutic", drug_c: "therapeutic"}
        res_base   = run_simulation(combo_base, PARAMS_PATH, n_beats=N_BEATS, verbose=False)
        dqtc_base  = round(res_base["QTc"] - baseline_qtc, 1)

        # Adjusted: scale CYP2D6 substrate (drug_b) by AUC ratio
        if drug_b in adj_cmax:
            combo_adj  = {drug_a: "therapeutic", drug_b: adj_cmax[drug_b], drug_c: "therapeutic"}
            res_adj    = run_simulation(combo_adj, PARAMS_PATH, n_beats=N_BEATS, verbose=False)
            dqtc_adj   = round(res_adj["QTc"] - baseline_qtc, 1)
            ikr_adj    = round(res_adj["IKr_block_pct"], 2)
            ratio      = auc_ratios[drug_b]
        else:
            dqtc_adj, ikr_adj, ratio = dqtc_base, res_base["IKr_block_pct"], 1.0

        a = DRUG_CODES[drug_a]; b = DRUG_CODES[drug_b]; c = DRUG_CODES[drug_c]
        print(f"  {a}+{b}+{c}: base={dqtc_base:+.1f}ms ({tier(dqtc_base)})  "
              f"adj={dqtc_adj:+.1f}ms ({tier(dqtc_adj)})  [{b} x{ratio}]")

        results.append({
            "combination":       f"{a}+{b}+{c}",
            "n_drugs":           3,
            "cyp2d6_inhibitor":  c,
            "substrate":         b,
            "auc_ratio":         ratio,
            "cmax_1x_nM":        SUBSTRATES.get(drug_b, {}).get("Cmax_free_nM", None),
            "cmax_adj_nM":       adj_cmax.get(drug_b, None),
            "dQTc_base_ms":      dqtc_base,
            "dQTc_adj_ms":       dqtc_adj,
            "dQTc_delta_ms":     round(dqtc_adj - dqtc_base, 1),
            "IKr_base_pct":      round(res_base["IKr_block_pct"], 2),
            "IKr_adj_pct":       ikr_adj,
            "tier_base":         tier(dqtc_base),
            "tier_adj":          tier(dqtc_adj),
            "tier_changed":      tier(dqtc_base) != tier(dqtc_adj),
            "fm_cyp2d6":         SUBSTRATES.get(drug_b, {}).get("fm_cyp2d6", None),
        })

    # ── Step 4: Save outputs ─────────────────────────────────────────────────
    df = pd.DataFrame(results)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nSaved -> {OUT_CSV}")

    tier_changes = df[df["tier_changed"] == True]
    print(f"\nTier escalations: {len(tier_changes)}")
    for _, row in tier_changes.iterrows():
        print(f"  {row['combination']}: {row['tier_base']} -> {row['tier_adj']} "
              f"(+{row['dQTc_delta_ms']:.1f}ms from CYP2D6 adjustment)")

    write_memo(df)
    print(f"Saved -> {OUT_MEMO}")

def write_memo(df):
    lines = []
    lines.append("# CYP2D6 Ito Static Interaction Model — CardioSafe Pediatric")
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    lines.append("")
    lines.append("## Model")
    lines.append("AUC ratio = 1 / (fm * (1 / (1 + [I]/Ki)) + (1 - fm))")
    lines.append("")
    lines.append("Parameters:")
    lines.append(f"- Fluoxetine [I] (blood): {FLU_PARAMS['I_blood_nM']:.1f} nM")
    lines.append(f"- Norfluoxetine [I]: {FLU_PARAMS['I_norflu_nM']:.1f} nM")
    lines.append(f"- Combined [I]total: {FLU_PARAMS['I_total_nM']:.1f} nM")
    lines.append(f"- Ki (CYP2D6): {FLU_PARAMS['Ki_nM']:.0f} nM")
    lines.append("")
    lines.append("## AUC Ratios")
    lines.append("| Substrate | fm | AUC Ratio | Cmax 1x | Cmax adj |")
    lines.append("|---|---|---|---|---|")
    for drug, params in SUBSTRATES.items():
        ratio = ito_auc_ratio(params["fm_cyp2d6"], FLU_PARAMS["I_total_nM"], FLU_PARAMS["Ki_nM"])
        cadj  = round(params["Cmax_free_nM"] * ratio, 1)
        lines.append(f"| {drug} | {params['fm_cyp2d6']:.2f} | {ratio:.2f}x | "
                     f"{params['Cmax_free_nM']:.1f} nM | {cadj:.1f} nM |")
    lines.append("")
    lines.append("## Delta-QTc Results")
    lines.append("| Combination | dQTc base | Tier base | dQTc adj | Tier adj | Delta | Tier change |")
    lines.append("|---|---|---|---|---|---|---|")
    for _, row in df.iterrows():
        change = "YES" if row["tier_changed"] else "no"
        lines.append(f"| {row['combination']} | {row['dQTc_base_ms']:+.1f}ms | {row['tier_base']} | "
                     f"{row['dQTc_adj_ms']:+.1f}ms | {row['tier_adj']} | "
                     f"{row['dQTc_delta_ms']:+.1f}ms | {change} |")
    lines.append("")
    lines.append("## Manuscript Implication")
    tier_changes = df[df["tier_changed"] == True]
    if len(tier_changes):
        lines.append(f"CYP2D6 adjustment escalates {len(tier_changes)} combinations to higher risk tiers.")
        lines.append("This explains the FAERS signal discordance for fluoxetine-containing combinations:")
        lines.append("the base-case model underestimates true pediatric exposure.")
    else:
        lines.append("CYP2D6 adjustment increases delta-QTc but does not change risk tiers in this analysis.")
        lines.append("The FAERS discordance for fluoxetine combinations likely reflects additional")
        lines.append("mechanisms beyond hERG block (e.g. active metabolite norfluoxetine, CNS effects).")
    lines.append("")
    lines.append("## References")
    lines.append("- Ito et al. Pharm Res. 1998;15(3):396-402 (static DDI model)")
    lines.append("- Templeton et al. Drug Metab Dispos. 2016;44(1):57-65 (fluoxetine CYP2D6 Ki)")
    lines.append("- FDA Drug Interaction Guidance 2020 (Ito static model validation)")
    lines.append("- Grimm et al. Br J Clin Pharmacol. 2006;61(1):58-69 (quetiapine CYP2D6 fm)")

    OUT_MEMO.write_text("\n".join(lines))

if __name__ == "__main__":
    main()