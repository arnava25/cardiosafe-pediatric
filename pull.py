"""
hERG IC50 Data Pull - CardioSafe Pediatric Project
Target: CHEMBL240 (hERG / KCNH2)
Run locally: pip install chembl-webresource-client pandas
"""

from chembl_webresource_client.new_client import new_client
import pandas as pd

DRUGS = {
    "Methylphenidate":   "CHEMBL796",
    "Amphetamine":       "CHEMBL405",
    "Dextroamphetamine": "CHEMBL1201070",
    "Aripiprazole":      "CHEMBL1112",
    "Quetiapine":        "CHEMBL675",
    "Risperidone":       "CHEMBL425",
    "Fluoxetine":        "CHEMBL41",
    "Sertraline":        "CHEMBL1046",
    "Escitalopram":      "CHEMBL1508",
    "Clonidine":         "CHEMBL600",
    "Guanfacine":        "CHEMBL1057",
    "Imipramine":        "CHEMBL11",
    "Nortriptyline":     "CHEMBL407",
}

HERG_TARGET    = "CHEMBL240"
ACTIVITY_TYPES = ["IC50", "Ki"]

activity_client = new_client.activity
records = []

print(f"Pulling hERG data for {len(DRUGS)} drugs...\n")

for drug_name, chembl_id in DRUGS.items():
    print(f"  {drug_name} ({chembl_id})...", end=" ")
    try:
        results = activity_client.filter(
            molecule_chembl_id=chembl_id,
            target_chembl_id=HERG_TARGET,
        ).only([
            "molecule_chembl_id",
            "standard_type",
            "standard_value",
            "standard_units",
            "pchembl_value",
            "assay_chembl_id",
            "assay_description",
            "bao_label",
            "activity_comment",
        ])

        count = 0
        for r in results:
            if r.get("standard_type") in ACTIVITY_TYPES:
                records.append({
                    "drug_name":        drug_name,
                    "chembl_id":        chembl_id,
                    "standard_type":    r.get("standard_type"),
                    "standard_value":   r.get("standard_value"),
                    "standard_units":   r.get("standard_units"),
                    "pchembl_value":    r.get("pchembl_value"),
                    "assay_id":         r.get("assay_chembl_id"),
                    "assay_description":str(r.get("assay_description", ""))[:150],
                    "assay_format":     r.get("bao_label"),
                    "activity_comment": r.get("activity_comment"),
                })
                count += 1
        print(f"{count} records")
    except Exception as e:
        print(f"ERROR: {e}")

df = pd.DataFrame(records)

if df.empty:
    print("\nNo data returned.")
else:
    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
    df["pchembl_value"]  = pd.to_numeric(df["pchembl_value"],  errors="coerce")

    def flag_assay(desc):
        desc = str(desc).lower()
        if any(x in desc for x in ["patch clamp", "electrophysiol", "qpatch", "manual patch", "voltage clamp"]):
            return "electrophysiology"
        elif any(x in desc for x in ["binding", "radioligand", "displacement", "dofetilide"]):
            return "binding_assay"
        elif any(x in desc for x in ["fluorescence", "thallium", "flux", "fret"]):
            return "fluorescence_flux"
        else:
            return "other"

    df["assay_quality"] = df["assay_description"].apply(flag_assay)

    summary = (
        df.groupby(["drug_name", "standard_type", "assay_quality"])
        .agg(
            n               = ("standard_value", "count"),
            median_nM       = ("standard_value", "median"),
            min_nM          = ("standard_value", "min"),
            max_nM          = ("standard_value", "max"),
            median_pchembl  = ("pchembl_value",  "median"),
        )
        .reset_index()
        .sort_values(["drug_name", "standard_type"])
    )

    print("\n" + "="*80)
    print("hERG BINDING SUMMARY (IC50/Ki in nM unless noted)")
    print("="*80)
    print(summary.to_string(index=False))

    df.to_csv("herg_full_data.csv", index=False)
    summary.to_csv("herg_summary.csv", index=False)

    print(f"\nSaved: herg_full_data.csv ({len(df)} rows)")
    print(f"Saved: herg_summary.csv ({len(summary)} rows)")
    print(f"\nAssay quality breakdown:\n{df['assay_quality'].value_counts().to_string()}")
    print(f"\nDrugs with NO hERG data: "
          f"{set(DRUGS.keys()) - set(df['drug_name'].unique())}")
