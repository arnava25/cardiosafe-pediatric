"""
hERG IC50 Data Pull v2 - CardioSafe Pediatric Project
Tries multiple ChEMBL IDs per drug (parent + salt forms)
Target: CHEMBL240 (hERG / KCNH2)
"""

from chembl_webresource_client.new_client import new_client
import pandas as pd

# Multiple IDs per drug to catch salt forms / alternative entries
DRUGS = {
    "Methylphenidate":   ["CHEMBL796",   "CHEMBL1201136"],
    "Amphetamine":       ["CHEMBL405",   "CHEMBL1706",   "CHEMBL2062428"],
    "Dextroamphetamine": ["CHEMBL1201070","CHEMBL405"],
    "Aripiprazole":      ["CHEMBL1112",  "CHEMBL1201244"],
    "Quetiapine":        ["CHEMBL675",   "CHEMBL1524",   "CHEMBL1201585"],
    "Risperidone":       ["CHEMBL425",   "CHEMBL1201168"],
    "Fluoxetine":        ["CHEMBL41",    "CHEMBL1201136"],
    "Sertraline":        ["CHEMBL1046",  "CHEMBL1201149"],
    "Escitalopram":      ["CHEMBL1508",  "CHEMBL1615467","CHEMBL1201188"],
    "Citalopram":        ["CHEMBL638",   "CHEMBL1201024"],  # parent of escitalopram, often has more data
    "Clonidine":         ["CHEMBL600",   "CHEMBL1697"],
    "Guanfacine":        ["CHEMBL1057",  "CHEMBL1201454"],
    "Imipramine":        ["CHEMBL11",    "CHEMBL1692"],
    "Nortriptyline":     ["CHEMBL407",   "CHEMBL1200686"],
}

HERG_TARGET    = "CHEMBL240"
ACTIVITY_TYPES = ["IC50", "Ki"]

activity_client = new_client.activity
records = []

print(f"Pulling hERG data for {len(DRUGS)} drugs (multiple IDs per drug)...\n")

for drug_name, chembl_ids in DRUGS.items():
    drug_records = []
    ids_tried = []
    for chembl_id in chembl_ids:
        ids_tried.append(chembl_id)
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
            for r in results:
                if r.get("standard_type") in ACTIVITY_TYPES:
                    drug_records.append({
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
        except Exception as e:
            print(f"    {chembl_id} error: {e}")

    # Deduplicate by assay_id + value
    seen = set()
    unique = []
    for r in drug_records:
        key = (r["assay_id"], r["standard_value"])
        if key not in seen:
            seen.add(key)
            unique.append(r)

    print(f"  {drug_name}: {len(unique)} unique records (tried {', '.join(ids_tried)})")
    records.extend(unique)

df = pd.DataFrame(records)

if df.empty:
    print("\nNo data returned.")
else:
    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
    df["pchembl_value"]  = pd.to_numeric(df["pchembl_value"],  errors="coerce")

    def flag_assay(desc):
        desc = str(desc).lower()
        if any(x in desc for x in ["patch clamp","electrophysiol","qpatch","voltage clamp","manual patch","hek","ik,"]):
            return "electrophysiology"
        elif any(x in desc for x in ["binding","radioligand","displacement","dofetilide","filter"]):
            return "binding_assay"
        elif any(x in desc for x in ["fluorescence","thallium","flux","fret","fluo"]):
            return "fluorescence_flux"
        else:
            return "other"

    df["assay_quality"] = df["assay_description"].apply(flag_assay)

    summary = (
        df.groupby(["drug_name", "standard_type", "assay_quality"])
        .agg(
            n              = ("standard_value", "count"),
            median_nM      = ("standard_value", "median"),
            min_nM         = ("standard_value", "min"),
            max_nM         = ("standard_value", "max"),
            median_pchembl = ("pchembl_value",  "median"),
        )
        .reset_index()
        .sort_values(["drug_name", "standard_type"])
    )

    print("\n" + "="*90)
    print("hERG BINDING SUMMARY v2 (IC50/Ki in nM)")
    print("="*90)
    print(summary.to_string(index=False))

    df.to_csv("herg_full_data_v2.csv", index=False)
    summary.to_csv("herg_summary_v2.csv", index=False)

    no_data = set(DRUGS.keys()) - set(df['drug_name'].unique())
    print(f"\nDrugs still with NO hERG data: {no_data}")
    print(f"\nTotal records: {len(df)}")
    print(f"\nAssay quality breakdown:\n{df['assay_quality'].value_counts().to_string()}")
