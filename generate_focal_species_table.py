import pandas as pd


raw = pd.read_csv("Raw_Multivariate_Regression_Results.csv")

focal_prefixes = ("Streptococcus_", "Veillonella_", "Rothia_")

# 1) Keep positive, nominally significant focal taxa hits
focal_hits = raw[
    raw["Species"].str.startswith(focal_prefixes, na=False)
    & (raw["Estimate"] > 0)
    & (raw["Pvalue"] <= 0.05)
].copy()

focal_hits = focal_hits[["Cohort", "Species", "Drug"]].rename(
    columns={"Drug": "Corresponding Medication"}
)

# 2) Combine, deduplicate, and sort
focal_table = (
    focal_hits
    .drop_duplicates()
    .sort_values(["Cohort", "Species", "Corresponding Medication"])
    .reset_index(drop=True)
)

# 3) Export the final CSV
focal_table.to_csv("Focal_Species_Cohort_Medication_Table.csv", index=False)

# Validation summary
raw_set = set(
    map(
        tuple,
        focal_hits[["Cohort", "Species", "Corresponding Medication"]].itertuples(index=False, name=None),
    )
)
final_set = set(
    map(
        tuple,
        focal_table[["Cohort", "Species", "Corresponding Medication"]].itertuples(index=False, name=None),
    )
)

print("=== Validation Summary ===")
print(f"Raw focal hits (positive & p<0.05): {len(raw_set)}")
print(f"Final rows written: {len(focal_table)}")
print(f"Duplicate rows in final table: {focal_table.duplicated().sum()}")
print(f"Rows in final not from raw: {len(final_set - raw_set)}")
print("Output file: Focal_Species_Cohort_Medication_Table.csv")
