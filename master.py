"""
Master hERG Parameterization Table - CardioSafe Pediatric Project
Consolidates ChEMBL pulls + literature values into one clean table
for O'Hara-Rudy IKr block simulation.

Drug effect on IKr modeled as:
    IKr_blocked = IKr_baseline * (1 / (1 + [C]_free / IC50_hERG))

Two cardiac risk pathways:
    1. hERG / IKr block  -> antipsychotics, SSRIs, tricyclics
    2. Sympathomimetic   -> stimulants (separate autonomic module)
    3. Mixed autonomic + mild IKr -> alpha-2 agonists
"""

import pandas as pd
import numpy as np

# ── CONSOLIDATED PARAMETERIZATION TABLE ─────────────────────────────────────
# ic50_nM:        best electrophysiology estimate
# ic50_low/high:  uncertainty bounds (from assay variability or source spread)
# cmax_free_nM:   free plasma Cmax at therapeutic dose (ADULT - flag for adolescent)
# primary_mechanism: determines which model pathway is used
# data_quality:   A = patch clamp confirmed, B = mixed/binding, C = estimate only

PARAMS = [
    # ── TRICYCLICS ────────────────────────────────────────────────────────────
    {
        "drug_name":          "Imipramine",
        "drug_class":         "TCA",
        "ic50_nM":            3388,
        "ic50_low_nM":        3000,
        "ic50_high_nM":       3400,
        "cmax_free_nM":       25.0,
        "primary_mechanism":  "hERG_block",
        "secondary_mechanism":"NE_reuptake + anticholinergic",
        "data_quality":       "A",
        "chembl_source":      "CHEMBL11 (patch clamp confirmed)",
        "lit_source":         "Witchel 2002; Redfern 2003",
        "adolescent_flag":    "Cmax adult-derived; CYP2D6 variability expected",
        "notes":              "Anchor drug; best characterized. Active metabolite desipramine also hERG active ~2000 nM",
    },
    {
        "drug_name":          "Nortriptyline",
        "drug_class":         "TCA",
        "ic50_nM":            1100,
        "ic50_low_nM":        1100,
        "ic50_high_nM":       1540,
        "cmax_free_nM":       30.0,
        "primary_mechanism":  "hERG_block",
        "secondary_mechanism":"NE_reuptake + anticholinergic",
        "data_quality":       "A",
        "chembl_source":      "sparse (binding only)",
        "lit_source":         "Witchel 2002; Perrin 2008",
        "adolescent_flag":    "Cmax adult-derived",
        "notes":              "More potent hERG block than imipramine despite same class",
    },

    # ── ANTIPSYCHOTICS ───────────────────────────────────────────────────────
    {
        "drug_name":          "Risperidone",
        "drug_class":         "Atypical antipsychotic",
        "ic50_nM":            560,
        "ic50_low_nM":        560,
        "ic50_high_nM":       700,
        "cmax_free_nM":       2.0,
        "primary_mechanism":  "hERG_block",
        "secondary_mechanism":"D2 + 5HT2A antagonism",
        "data_quality":       "A",
        "chembl_source":      "not indexed under parent ID",
        "lit_source":         "Kongsamut 2002; Kramer 2013",
        "adolescent_flag":    "Cmax adult-derived; active metabolite 9-OH-risperidone also QTc active",
        "notes":              "Lowest IC50 among antipsychotics; highest intrinsic hERG risk in class. Safety index ~280x but active metabolite adds",
    },
    {
        "drug_name":          "Quetiapine",
        "drug_class":         "Atypical antipsychotic",
        "ic50_nM":            1070,
        "ic50_low_nM":        1070,
        "ic50_high_nM":       2300,
        "cmax_free_nM":       10.0,
        "primary_mechanism":  "hERG_block",
        "secondary_mechanism":"D2 + H1 + alpha1 antagonism",
        "data_quality":       "A",
        "chembl_source":      "not indexed (known gap)",
        "lit_source":         "Kongsamut 2002; Vieweg 2004",
        "adolescent_flag":    "Cmax adult-derived; heavy off-label use in adolescents",
        "notes":              "Wide assay variability (2x spread binding vs patch clamp). Off-label pediatric use very common",
    },
    {
        "drug_name":          "Aripiprazole",
        "drug_class":         "Atypical antipsychotic",
        "ic50_nM":            1135,   # mean of 896 and 1350 from patch clamp
        "ic50_low_nM":        896,
        "ic50_high_nM":       2378,
        "cmax_free_nM":       50.0,
        "primary_mechanism":  "hERG_block",
        "secondary_mechanism":"D2 partial agonist + 5HT1A",
        "data_quality":       "A",
        "chembl_source":      "CHEMBL1112 (patch clamp confirmed)",
        "lit_source":         "Kramer 2013; Perrin 2008",
        "adolescent_flag":    "Cmax adult-derived; FDA approved pediatric use (bipolar/irritability)",
        "notes":              "Best hERG safety profile among antipsychotics. High protein binding means free Cmax low",
    },

    # ── SSRIs ────────────────────────────────────────────────────────────────
    {
        "drug_name":          "Sertraline",
        "drug_class":         "SSRI",
        "ic50_nM":            518,    # median of 446 and 590
        "ic50_low_nM":        446,
        "ic50_high_nM":       590,
        "cmax_free_nM":       1.0,
        "primary_mechanism":  "hERG_block",
        "secondary_mechanism":"SERT inhibition",
        "data_quality":       "A",
        "chembl_source":      "sparse",
        "lit_source":         "Polak 2009; Kramer 2013",
        "adolescent_flag":    "Cmax adult-derived; commonly used in adolescents",
        "notes":              "Most potent hERG block among SSRIs by IC50 but very high protein binding drives safety index up. Still worth modeling in combination",
    },
    {
        "drug_name":          "Fluoxetine",
        "drug_class":         "SSRI",
        "ic50_nM":            1513,
        "ic50_low_nM":        1513,
        "ic50_high_nM":       3162,
        "cmax_free_nM":       5.0,
        "primary_mechanism":  "hERG_block",
        "secondary_mechanism":"SERT inhibition",
        "data_quality":       "A",
        "chembl_source":      "CHEMBL41 (patch clamp confirmed)",
        "lit_source":         "Polak 2009; Tarantino 2011",
        "adolescent_flag":    "FDA approved pediatric depression age 8+; Cmax adult-derived",
        "notes":              "Active metabolite norfluoxetine IC50 ~3000 nM; long half-life means accumulation. Must model parent + metabolite",
    },
    {
        "drug_name":          "Escitalopram",
        "drug_class":         "SSRI",
        "ic50_nM":            9500,
        "ic50_low_nM":        9500,
        "ic50_high_nM":       14800,
        "cmax_free_nM":       3.0,
        "primary_mechanism":  "hERG_block",
        "secondary_mechanism":"SERT inhibition",
        "data_quality":       "B",
        "chembl_source":      "not indexed",
        "lit_source":         "Polak 2009; Tarantino 2011",
        "adolescent_flag":    "Cmax adult-derived; known QTc prolongation in overdose",
        "notes":              "Weak single-drug hERG but FDA has QTc warning. May have hERG-independent QTc mechanism. Flag for polypharmacy modeling",
    },

    # ── ALPHA-2 AGONISTS ─────────────────────────────────────────────────────
    {
        "drug_name":          "Clonidine",
        "drug_class":         "Alpha-2 agonist",
        "ic50_nM":            17378,
        "ic50_low_nM":        17378,
        "ic50_high_nM":       28000,
        "cmax_free_nM":       0.5,
        "primary_mechanism":  "autonomic_modulation",
        "secondary_mechanism":"weak_hERG_block",
        "data_quality":       "B",
        "chembl_source":      "CHEMBL600 (1 record, unspecified assay)",
        "lit_source":         "Redfern 2003",
        "adolescent_flag":    "Very commonly used in pediatric ADHD/aggression; Cmax adult-derived",
        "notes":              "Primary cardiac risk is bradycardia + conduction slowing via alpha-2, not QTc. Model via autonomic pathway. hERG block negligible at therapeutic Cmax",
    },
    {
        "drug_name":          "Guanfacine",
        "drug_class":         "Alpha-2 agonist",
        "ic50_nM":            50000,  # estimate; sparse data
        "ic50_low_nM":        30000,
        "ic50_high_nM":       70000,
        "cmax_free_nM":       0.3,
        "primary_mechanism":  "autonomic_modulation",
        "secondary_mechanism":"weak_hERG_block",
        "data_quality":       "C",
        "chembl_source":      "no records",
        "lit_source":         "Perrin 2008 (estimate)",
        "adolescent_flag":    "FDA approved pediatric ADHD; Cmax adult-derived",
        "notes":              "IC50 is estimate only - flag as high uncertainty. Same autonomic mechanism as clonidine. PR prolongation is primary ECG concern",
    },

    # ── STIMULANTS ───────────────────────────────────────────────────────────
    {
        "drug_name":          "Methylphenidate",
        "drug_class":         "Stimulant",
        "ic50_nM":            100000,
        "ic50_low_nM":        100000,
        "ic50_high_nM":       200000,
        "cmax_free_nM":       5.0,
        "primary_mechanism":  "sympathomimetic",
        "secondary_mechanism":"negligible_hERG",
        "data_quality":       "C",
        "chembl_source":      "no records",
        "lit_source":         "Darpö 2001",
        "adolescent_flag":    "Most prescribed psychiatric drug in children; Cmax adult-derived",
        "notes":              "Cardiac risk via catecholamine release: HR increase, BP increase. Model via beta-adrenergic pathway (ICaL, If modulation). hERG block negligible",
    },
    {
        "drug_name":          "Amphetamine",
        "drug_class":         "Stimulant",
        "ic50_nM":            110000,
        "ic50_low_nM":        110000,
        "ic50_high_nM":       200000,
        "cmax_free_nM":       20.0,
        "primary_mechanism":  "sympathomimetic",
        "secondary_mechanism":"negligible_hERG",
        "data_quality":       "C",
        "chembl_source":      "no records",
        "lit_source":         "Darpö 2001",
        "adolescent_flag":    "Cmax adult-derived",
        "notes":              "Same sympathomimetic pathway as methylphenidate but more pronounced catecholamine release. Model HR/BP effects not QTc",
    },
]

df = pd.DataFrame(PARAMS)

# Compute safety index
df["safety_index"]   = df["ic50_nM"] / df["cmax_free_nM"]
df["safety_index_low"]  = df["ic50_low_nM"]  / df["cmax_free_nM"]
df["safety_index_high"] = df["ic50_high_nM"] / df["cmax_free_nM"]

df["risk_hERG_single"] = df["safety_index"].apply(
    lambda x: "HIGH" if x < 10 else ("MODERATE" if x < 30 else "LOW")
)

# Fractional IKr block at therapeutic Cmax
# block = C / (C + IC50); C = cmax_free
df["IKr_block_pct"] = 100 * df["cmax_free_nM"] / (df["cmax_free_nM"] + df["ic50_nM"])
df["IKr_block_pct"] = df["IKr_block_pct"].round(2)

print("="*100)
print("MASTER hERG PARAMETERIZATION TABLE - CardioSafe Pediatric")
print("="*100)

display_cols = [
    "drug_name","drug_class","ic50_nM","ic50_low_nM","ic50_high_nM",
    "cmax_free_nM","IKr_block_pct","safety_index","risk_hERG_single",
    "primary_mechanism","data_quality"
]
print(df[display_cols].to_string(index=False))

print("\n" + "="*100)
print("IKr BLOCK % AT THERAPEUTIC Cmax (single drug)")
print("  -> This is the conductance scaling factor for O'Hara-Rudy IKr")
print("  -> GKr_effective = GKr_max * (1 - IKr_block_pct/100)")
print("="*100)
for _, row in df.sort_values("IKr_block_pct", ascending=False).iterrows():
    bar = "█" * int(row["IKr_block_pct"] / 2)
    print(f"  {row['drug_name']:20s} {row['IKr_block_pct']:5.2f}%  {bar}")

print("\n" + "="*100)
print("DATA QUALITY FLAGS")
print("  A = patch clamp electrophysiology confirmed")
print("  B = mixed sources / binding assay primary")
print("  C = estimate only - propagate as high uncertainty in model")
print("="*100)
for q in ["A","B","C"]:
    drugs = df[df["data_quality"]==q]["drug_name"].tolist()
    print(f"  Quality {q}: {', '.join(drugs)}")

print("\n" + "="*100)
print("MODEL PATHWAY ASSIGNMENT")
print("="*100)
for mech in df["primary_mechanism"].unique():
    drugs = df[df["primary_mechanism"]==mech]["drug_name"].tolist()
    print(f"  {mech}: {', '.join(drugs)}")

print("\nADOLESCENT FLAGS (all Cmax values are adult-derived - key limitation):")
for _, row in df.iterrows():
    print(f"  {row['drug_name']:20s}: {row['adolescent_flag']}")

# Save
df.to_csv("herg_master_params.csv", index=False)
print(f"\nSaved: herg_master_params.csv")
print("\nNext step: Feed ic50_nM and IKr_block_pct into O'Hara-Rudy IKr block module")
print("           Build polypharmacy combinations as additive block: ")
print("           total_block = 1 - product((1 - block_i) for each drug i)")
