# ==============================================================================
# 100% PURE PYTHON PIPELINE: Age Constraints, Regressions & Visualizations
# ==============================================================================
import os
import time
import pyreadr
import pandas as pd
import numpy as np
import statsmodels.api as sm
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import seaborn as sns

print("==================================================")
print("PART 1: LOADING .RData & PREPARING DATA")
print("==================================================")
start_time = time.time()

# 1. Load Data
print("Loading 'Required_Data_FEBS_Pipepline.RData'...")
try:
    rdata = pyreadr.read_r("Required_Data_FEBS_Pipepline.RData")
except Exception as e:
    raise FileNotFoundError("Could not read .RData file. Please ensure it is in the same directory.")

AllCombinedMetadata = rdata["AllCombinedMetadata"]
AllCombinedSpProfile = rdata["AllCombinedSpProfile"]
SpeciesScores = rdata["SpeciesScores_20241021"]

# Ensure age is numeric
AllCombinedMetadata['age_num'] = pd.to_numeric(AllCombinedMetadata['age'], errors='coerce')

print("==================================================")
print("PART A: FINDING ELIGIBLE STUDIES (New Age Constraints)")
print("==================================================")
# Mimic R Part A: Find studies with Min 41-45 and Max 75-80 (from 41-80 cohort)
age_filtered_for_summary = AllCombinedMetadata[(AllCombinedMetadata['age_num'] > 40) & (AllCombinedMetadata['age_num'] <= 80)]

study_age_ranges = age_filtered_for_summary.groupby('study_name')['age_num'].agg(
    Overall_Min_Age='min', Overall_Max_Age='max'
).reset_index()

eligible_studies_df = study_age_ranges[
    (study_age_ranges['Overall_Min_Age'] >= 41) & (study_age_ranges['Overall_Min_Age'] <= 45) &
    (study_age_ranges['Overall_Max_Age'] >= 75) & (study_age_ranges['Overall_Max_Age'] <= 80)
]
eligible_studies_list = eligible_studies_df['study_name'].tolist()
print(f"Found {len(eligible_studies_list)} eligible studies matching constraints.")

print("==================================================")
print("PART B: PREPARING COHORTS")
print("==================================================")
# Taxonomic Renaming
rename_map = {
    "Ruminococcus_obeum": "Blautia_obeum",
    "Eubacterium_biforme": "Holdemanella_biformis",
    "Clostridium_bartlettii": "Intestinibacter_bartlettii",
    "Clostridium_hathewayi": "Hungatella_hathewayi",
    "Clostridium_ramosum": "Erysipelatoclostridium_ramosum"
}
AllCombinedSpProfile.rename(columns=rename_map, inplace=True)

# THE FIX: Drop any duplicate columns created by the renaming
AllCombinedSpProfile = AllCombinedSpProfile.loc[:, ~AllCombinedSpProfile.columns.duplicated()]

# Filter 201 Taxa
all_201_taxa_names = SpeciesScores.index.tolist()
common_taxa = [t for t in all_201_taxa_names if t in AllCombinedSpProfile.columns]
filtered_profile = AllCombinedSpProfile[common_taxa]

# Merge Meta and Profile
selected_meta = AllCombinedMetadata[['age_num', 'study_name', 'diseaseCat']]
master_data = pd.merge(selected_meta, filtered_profile, left_index=True, right_index=True)

# Filter by Eligible Studies and Age 41-80
master_data = master_data[master_data['study_name'].isin(eligible_studies_list)]
master_data = master_data[(master_data['age_num'] > 40) & (master_data['age_num'] <= 80)]

# Split Control & Disease
control_cohort = master_data[master_data['diseaseCat'] == 'control'].copy()
disease_cohort = master_data[master_data['diseaseCat'] != 'control'].copy()

# Un-nest multi-disease rows (Equivalent to separate_rows in R)
disease_cohort = disease_cohort.assign(disease_name=disease_cohort['diseaseCat'].str.split(';')).explode('disease_name')

target_disease_names = ["CRC", "CVD", "NAFLD_FLD", "Parkinsons", "Polyps", "T2D"]
disease_cohort = disease_cohort[disease_cohort['disease_name'].isin(target_disease_names)]

common_studies = list(set(control_cohort['study_name']).intersection(set(disease_cohort['study_name'])))

# Filter for >= 3 Studies
MIN_STUDIES = 3
final_diseases = {}
for disease in target_disease_names:
    d_df = disease_cohort[disease_cohort['disease_name'] == disease]
    d_studies = set(d_df['study_name']).intersection(common_studies)
    if len(d_studies) >= MIN_STUDIES:
        final_diseases[disease] = list(d_studies)

print(f"Kept {len(final_diseases)} diseases with >= 3 studies.")

print("==================================================")
print("PART C: RUNNING PARALLEL REGRESSIONS (WHITEBOARD LOGIC)")
print("==================================================")

# --- 1. Disease + Age Regression ---
print("Running 'Species ~ Status + Age' models...")
def run_disease_model(disease, study, taxon):
    d_df = disease_cohort[(disease_cohort['disease_name'] == disease) & (disease_cohort['study_name'] == study)].copy()
    d_df['Status'] = 1
    c_df = control_cohort[control_cohort['study_name'] == study].copy()
    c_df['Status'] = 0
    
    combined = pd.concat([d_df, c_df]).dropna(subset=['Status', 'age_num', taxon])
    
    if len(d_df) == 0 or len(c_df) == 0 or combined[taxon].nunique() <= 1 or combined['Status'].nunique() <= 1:
        return (disease, study, taxon, 0)
    
    X = sm.add_constant(combined[['Status', 'age_num']])
    try:
        model = sm.OLS(combined[taxon], X).fit()
        if model.pvalues.get('Status', 1.0) < 0.1:
            return (disease, study, taxon, np.sign(model.params.get('Status', 0)))
    except: pass
    return (disease, study, taxon, 0)

disease_tasks = [(d, s, t) for d in final_diseases for s in final_diseases[d] for t in common_taxa]
disease_results = Parallel(n_jobs=-1)(delayed(run_disease_model)(d, s, t) for d, s, t in disease_tasks)
disease_res_df = pd.DataFrame(disease_results, columns=['Disease', 'Study', 'Taxon', 'Direction'])

# --- 2. Age Regression (Controls) ---
print("Running 'Species ~ Age' models in Controls...")
def run_age_model(study, taxon):
    c_df = control_cohort[control_cohort['study_name'] == study].dropna(subset=['age_num', taxon])
    if len(c_df) == 0 or c_df[taxon].nunique() <= 1 or c_df['age_num'].nunique() <= 1:
        return (study, taxon, 0)
    
    X = sm.add_constant(c_df[['age_num']])
    try:
        model = sm.OLS(c_df[taxon], X).fit()
        if model.pvalues.get('age_num', 1.0) < 0.1:
            return (study, taxon, np.sign(model.params.get('age_num', 0)))
    except: pass
    return (study, taxon, 0)

control_studies_all = list(set([s for studies in final_diseases.values() for s in studies]))
age_tasks = [(s, t) for s in control_studies_all for t in common_taxa]
age_results = Parallel(n_jobs=-1)(delayed(run_age_model)(s, t) for s, t in age_tasks)
age_res_df = pd.DataFrame(age_results, columns=['Study', 'Taxon', 'Direction'])

print("==================================================")
print("PART D: ASSOCIATION SCORES & FINAL MATRIX")
print("==================================================")
def calc_as(matrix_df):
    SP = (matrix_df == 1).sum(axis=1)
    SN = (matrix_df == -1).sum(axis=1)
    T = matrix_df.shape[1]
    pmax = np.maximum(SP, SN) + 0.00001
    pmin = np.minimum(SP, SN) + 0.00001
    consensus = 1 - (pmin / pmax)
    consensus[(SP == 0) & (SN == 0)] = 1.0
    return ((SP - SN) / T) * consensus

# Age Scores
age_matrix = age_res_df.pivot(index='Taxon', columns='Study', values='Direction').fillna(0)
final_scores = {'age_association_controls': calc_as(age_matrix)}

# Disease Scores
for d in final_diseases:
    d_matrix = disease_res_df[disease_res_df['Disease'] == d].pivot(index='Taxon', columns='Study', values='Direction').fillna(0)
    final_scores[d] = calc_as(d_matrix)

df_combined_scores = pd.DataFrame(final_scores).reindex(all_201_taxa_names).fillna(0)
df_combined_scores.to_csv("Final_Association_Scores_Age_Adjusted.csv")

print("==================================================")
print("PART E: GENERATING PLOTS")
print("==================================================")
output_dir = "individual_disease_plots_NO_LEGEND_Python"
os.makedirs(output_dir, exist_ok=True)

color_map = {
    "Positive Concordance": "mediumslateblue",
    "Negative Concordance": "deeppink",
    "Discordant": "khaki",
    "On Axis": "gray"
}

for disease in final_diseases.keys():
    plt.figure(figsize=(8, 6))
    
    x = df_combined_scores['age_association_controls']
    y = df_combined_scores[disease]
    
    colors = []
    for ix, iy in zip(x, y):
        if ix > 0 and iy > 0: colors.append(color_map["Positive Concordance"])
        elif ix < 0 and iy < 0: colors.append(color_map["Negative Concordance"])
        elif np.sign(ix) != np.sign(iy) and ix != 0 and iy != 0: colors.append(color_map["Discordant"])
        else: colors.append(color_map["On Axis"])

    plt.scatter(x, y, c=colors, s=100, edgecolors='none', alpha=0.9)
    plt.axhline(0, color='black', linewidth=2)
    plt.axvline(0, color='black', linewidth=2)
    
    plt.xlabel("Age Association Score (Controls)", fontsize=14, fontweight='bold')
    plt.ylabel(f"Association Score ({disease})", fontsize=14, fontweight='bold')
    plt.title(f"Age vs {disease}", fontsize=16, fontweight='bold')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    
    file_path = os.path.join(output_dir, f"plot_age_vs_{disease}.png")
    plt.savefig(file_path, dpi=300)
    plt.close()
    print(f"Saved plot: {file_path}")

print("==================================================")
print(f"SUCCESS! Entire pipeline finished in {round(time.time() - start_time, 2)} seconds.")
