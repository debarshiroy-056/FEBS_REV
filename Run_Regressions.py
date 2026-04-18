# ==============================================================================
# SCRIPT 1 (PYTHON): Parallel Regression & Japanese-4D Integration
# ==============================================================================
import pandas as pd
import numpy as np
import statsmodels.api as sm
from joblib import Parallel, delayed
import time
import os


def load_all_abundance_df():
    """Load AllCombinedSpProfile from CSV, or fall back to RData when CSV is empty."""
    csv_path = "AllCombinedSpProfile.csv"

    if os.path.exists(csv_path):
        try:
            csv_df = pd.read_csv(csv_path, index_col=0)
            if not csv_df.empty and csv_df.shape[1] > 0:
                print(f"Loaded abundance matrix from {csv_path}")
                return csv_df.fillna(0)
            print(f"{csv_path} exists but is empty; falling back to RData...")
        except pd.errors.EmptyDataError:
            print(f"{csv_path} is empty; falling back to RData...")
    else:
        print(f"{csv_path} not found; falling back to RData...")

    try:
        import pyreadr
    except ImportError as exc:
        raise ImportError(
            "pyreadr is required to read AllCombinedSpProfile from RData. "
            "Install with: pip install pyreadr"
        ) from exc

    rdata_candidates = [
        "Required_Data_FEBS_Pipeline.RData",
        "Required_Data_FEBS_Pipepline.RData"
    ]

    for rdata_path in rdata_candidates:
        if not os.path.exists(rdata_path):
            continue

        print(f"Trying to load abundance matrix from {rdata_path}...")
        r_objects = pyreadr.read_r(rdata_path)

        if "AllCombinedSpProfile" not in r_objects:
            continue

        all_abundance_df = r_objects["AllCombinedSpProfile"]
        if not isinstance(all_abundance_df, pd.DataFrame) or all_abundance_df.empty:
            raise ValueError(
                f"AllCombinedSpProfile in {rdata_path} is missing or empty."
            )

        print(f"Loaded abundance matrix from {rdata_path} (object: AllCombinedSpProfile)")
        return all_abundance_df.fillna(0)

    raise FileNotFoundError(
        "Could not load AllCombinedSpProfile. Checked AllCombinedSpProfile.csv and "
        "RData files: Required_Data_FEBS_Pipeline.RData / Required_Data_FEBS_Pipepline.RData"
    )


def load_target_42_drugs(xlsx_path="42_medications.xlsx"):
    """Extract the 42 medication names from the Japanese-4D species-level sheet."""
    japanese_sheet_candidates = ["S9", "Our genus_medication_multi"]
    jap_raw = None
    used_sheet = None

    for sheet_name in japanese_sheet_candidates:
        try:
            jap_raw = pd.read_excel(xlsx_path, sheet_name=sheet_name, header=None)
            used_sheet = sheet_name
            break
        except ValueError:
            continue

    if jap_raw is None:
        raise ValueError(
            "Could not find a usable Japanese-4D worksheet in 42_medications.xlsx. "
            "Tried: " + ", ".join(japanese_sheet_candidates)
        )

    drugs = (
        pd.Series(jap_raw.iloc[0, 1:85:2])
        .fillna("")
        .astype(str)
        .str.strip()
    )
    drugs = [d for d in drugs if d and d.lower() != "nan"]

    if len(drugs) != 42:
        raise ValueError(f"Expected 42 medications from sheet '{used_sheet}', found {len(drugs)}")

    print(f"Loaded 42 target medications from {xlsx_path} (sheet: {used_sheet})")
    return drugs, jap_raw, used_sheet

print("==================================================")
print("PART 1: PARALLEL MULTIVARIATE REGRESSION (MetaCardis)")
print("==================================================")
start_time = time.time()

# Load and freeze the exact 42-drug panel used for Figure 4 across both cohorts.
target_42_drugs, jap_raw, used_japanese_sheet = load_target_42_drugs("42_medications.xlsx")

# --- 1. Load MetaCardis Data ---
all_abundance_df = load_all_abundance_df()
drug_matrix = pd.read_csv("drug_matrix.csv").fillna(0)

# Make column names unique
cols = pd.Series(all_abundance_df.columns)
for dup in cols[cols.duplicated()].unique(): 
    cols[cols[cols == dup].index.values.tolist()] = [dup + '.' + str(i) if i != 0 else dup for i in range(sum(cols == dup))]
all_abundance_df.columns = cols

all_species = all_abundance_df.columns.tolist()
all_abundance_df['Sample'] = all_abundance_df.index
available_meta_drugs = [col for col in drug_matrix.columns if col != 'Sample']
selected_drugs = [d for d in target_42_drugs if d in available_meta_drugs]

if len(selected_drugs) != 42:
    missing_from_meta = sorted(set(target_42_drugs) - set(selected_drugs))
    raise ValueError(
        "Figure 4 requires all 42 target medications in drug_matrix.csv. Missing: "
        + ", ".join(missing_from_meta)
    )

print(f"Using {len(selected_drugs)} harmonized medications for MetaCardis regression.")

analysis_df_all = pd.merge(all_abundance_df, drug_matrix[['Sample'] + selected_drugs], on='Sample', how='inner')

# --- 2. Parallel Regression ---
print(f"Running highly parallel OLS models for {len(all_species)} species...")
X = analysis_df_all[selected_drugs]
X = sm.add_constant(X) 
X_np = X.values 
drug_names = X.columns.tolist()

def fit_ols_fast(sp_name, y_values):
    model = sm.OLS(y_values, X_np).fit()
    return sp_name, model.params, model.pvalues

parallel_results = Parallel(n_jobs=-1)(
    delayed(fit_ols_fast)(sp, analysis_df_all[sp].values) for sp in all_species
)

results_list = []
for sp_name, params, pvals in parallel_results:
    for idx, drug in enumerate(drug_names):
        if drug == 'const': continue
        results_list.append({
            'Species': sp_name,
            'Drug': drug,
            'Estimate': params[idx],
            'Pvalue': pvals[idx],
            'Cohort': 'MetaCardis'
        })

all_results_mc = pd.DataFrame(results_list).dropna(subset=['Pvalue'])
print("-> MetaCardis models complete and NA values removed!")

print("==================================================")
print("PART 2: PARSING JAPANESE-4D EXCEL DATA")
print("==================================================")

# --- 3. Extract Japanese-4D from Excel ---
print(f"Using Japanese-4D sheet: {used_japanese_sheet} (shape={jap_raw.shape})")

# Row index 0 has the drug names. Column 0 has the species.
drug_names_j4d = jap_raw.iloc[0, 1:85:2].fillna("").astype(str).values
jap_species = (
    jap_raw.iloc[3:, 0]
    .fillna("")
    .astype(str)
    .str.strip()
    .str.replace(" ", "_", regex=False)
    .values
)

# Remove empty/non-data species labels while preserving row alignment.
valid_species_mask = np.array([
    (sp != "") and (sp != "nan") and (not sp.startswith("Name_as_per"))
    for sp in jap_species
])

j4d_list = []
for i, drug in enumerate(drug_names_j4d):
    if drug in selected_drugs:
        est_col = 1 + 2*i
        pval_col = 2 + 2*i
        est_vals = pd.to_numeric(jap_raw.iloc[3:, est_col], errors='coerce').values
        pval_vals = pd.to_numeric(jap_raw.iloc[3:, pval_col], errors='coerce').values

        species_filtered = jap_species[valid_species_mask]
        est_filtered = est_vals[valid_species_mask]
        pval_filtered = pval_vals[valid_species_mask]
        
        df_j4d = pd.DataFrame({
            'Species': species_filtered,
            'Drug': drug,
            'Estimate': est_filtered,
            'Pvalue': pval_filtered,
            'Cohort': 'Japanese-4D'
        })
        j4d_list.append(df_j4d)

all_results_j4d = pd.concat(j4d_list, ignore_index=True).dropna(subset=['Pvalue'])
print(f"-> Japanese-4D parsed successfully! Found {len(all_results_j4d)} valid NA-free models.")

# --- 4. Merge and Export ---
print("==================================================")
print("PART 3: EXPORTING CLEANED RESULTS")
print("==================================================")

# Combine both cohorts into a master file (used by the R script)
all_results_combined = pd.concat([all_results_mc, all_results_j4d], ignore_index=True)
all_results_combined.to_csv("All_Cohorts_Regression_Results.csv", index=False)

# Filter out the PI's target species (Streptococcus, Veillonella, Rothia)
target_species_all = [sp for sp in all_results_combined['Species'].unique() if sp.startswith(('Streptococcus_', 'Veillonella_', 'Rothia_'))]
pi_stats = all_results_combined[all_results_combined['Species'].isin(target_species_all)].sort_values(by=['Cohort', 'Pvalue'])

# Export the CSV for the Professor
pi_stats.to_csv("Raw_Multivariate_Regression_Results.csv", index=False)

print("-> Saved 'All_Cohorts_Regression_Results.csv' (For R Script Heatmap)")
print("-> Saved 'Raw_Multivariate_Regression_Results.csv' (For Professor Ghosh)")
print(f"Python processing complete in {round(time.time() - start_time, 2)} seconds!")