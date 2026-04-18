# ==============================================================================
# PURE PYTHON PIPELINE: Unadjusted vs Age-Adjusted Association Scores
# ==============================================================================
import os
import time
import pyreadr
import pandas as pd
import numpy as np
import statsmodels.api as sm
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
from scipy.stats import pearsonr  # <-- NEW IMPORT ADDED HERE

print("==================================================")
print("PART 1: LOADING & PREPARING DATA")
print("==================================================")
start_time = time.time()

# 1. Load Data
try:
    rdata = pyreadr.read_r("Required_Data_FEBS_Pipepline.RData")
except Exception as e:
    raise FileNotFoundError("Could not read .RData file.")

AllCombinedMetadata = rdata["AllCombinedMetadata"]
AllCombinedSpProfile = rdata["AllCombinedSpProfile"]
SpeciesScores = rdata["SpeciesScores_20241021"]

AllCombinedMetadata['age_num'] = pd.to_numeric(AllCombinedMetadata['age'], errors='coerce')

# 2. Mimic R Part A: Find eligible studies
age_filtered_for_summary = AllCombinedMetadata[(AllCombinedMetadata['age_num'] > 40) & (AllCombinedMetadata['age_num'] <= 80)]
study_age_ranges = age_filtered_for_summary.groupby('study_name')['age_num'].agg(Overall_Min_Age='min', Overall_Max_Age='max').reset_index()

eligible_studies_df = study_age_ranges[
    (study_age_ranges['Overall_Min_Age'] >= 41) & (study_age_ranges['Overall_Min_Age'] <= 45) &
    (study_age_ranges['Overall_Max_Age'] >= 75) & (study_age_ranges['Overall_Max_Age'] <= 80)
]
eligible_studies_list = eligible_studies_df['study_name'].tolist()

# 3. Rename and clean profile
rename_map = {
    "Ruminococcus_obeum": "Blautia_obeum",
    "Eubacterium_biforme": "Holdemanella_biformis",
    "Clostridium_bartlettii": "Intestinibacter_bartlettii",
    "Clostridium_hathewayi": "Hungatella_hathewayi",
    "Clostridium_ramosum": "Erysipelatoclostridium_ramosum"
}
AllCombinedSpProfile.rename(columns=rename_map, inplace=True)
AllCombinedSpProfile = AllCombinedSpProfile.loc[:, ~AllCombinedSpProfile.columns.duplicated()]

common_taxa = [t for t in SpeciesScores.index.tolist() if t in AllCombinedSpProfile.columns]
filtered_profile = AllCombinedSpProfile[common_taxa]

# 4. Merge and Split
selected_meta = AllCombinedMetadata[['age_num', 'study_name', 'diseaseCat']]
master_data = pd.merge(selected_meta, filtered_profile, left_index=True, right_index=True)
master_data = master_data[master_data['study_name'].isin(eligible_studies_list)]
master_data = master_data[(master_data['age_num'] > 40) & (master_data['age_num'] <= 80)]

control_cohort = master_data[master_data['diseaseCat'] == 'control'].copy()
disease_cohort = master_data[master_data['diseaseCat'] != 'control'].copy()
disease_cohort = disease_cohort.assign(disease_name=disease_cohort['diseaseCat'].str.split(';')).explode('disease_name')

target_disease_names = ["CRC", "CVD", "NAFLD_FLD", "Parkinsons", "Polyps", "T2D"]
disease_cohort = disease_cohort[disease_cohort['disease_name'].isin(target_disease_names)]
common_studies = list(set(control_cohort['study_name']).intersection(set(disease_cohort['study_name'])))

final_diseases = {}
for disease in target_disease_names:
    d_df = disease_cohort[disease_cohort['disease_name'] == disease]
    d_studies = set(d_df['study_name']).intersection(common_studies)
    if len(d_studies) >= 3:
        final_diseases[disease] = list(d_studies)

print("==================================================")
print("PART 2: RUNNING PARALLEL REGRESSIONS (UNADJUSTED vs ADJUSTED)")
print("==================================================")

def run_both_models(disease, study, taxon):
    d_df = disease_cohort[(disease_cohort['disease_name'] == disease) & (disease_cohort['study_name'] == study)].copy()
    d_df['Status'] = 1
    c_df = control_cohort[control_cohort['study_name'] == study].copy()
    c_df['Status'] = 0
    
    # Drop NA so BOTH models use the exact same N samples
    combined = pd.concat([d_df, c_df]).dropna(subset=['Status', 'age_num', taxon])
    
    if len(d_df) == 0 or len(c_df) == 0 or combined[taxon].nunique() <= 1 or combined['Status'].nunique() <= 1:
        return (disease, study, taxon, 0, 0)
    
    # 1. Unadjusted Model: Species ~ Status
    X_unadj = sm.add_constant(combined[['Status']])
    dir_unadj = 0
    try:
        mod_unadj = sm.OLS(combined[taxon], X_unadj).fit()
        if mod_unadj.pvalues.get('Status', 1.0) < 0.1:
            dir_unadj = np.sign(mod_unadj.params.get('Status', 0))
    except: pass
    
    # 2. Age-Adjusted Model: Species ~ Status + Age
    X_adj = sm.add_constant(combined[['Status', 'age_num']])
    dir_adj = 0
    try:
        mod_adj = sm.OLS(combined[taxon], X_adj).fit()
        if mod_adj.pvalues.get('Status', 1.0) < 0.1:
            dir_adj = np.sign(mod_adj.params.get('Status', 0))
    except: pass
    
    return (disease, study, taxon, dir_unadj, dir_adj)

tasks = [(d, s, t) for d in final_diseases for s in final_diseases[d] for t in common_taxa]
print(f"Running {len(tasks)} dual-regression tasks across all CPU cores...")

results = Parallel(n_jobs=-1)(delayed(run_both_models)(d, s, t) for d, s, t in tasks)
res_df = pd.DataFrame(results, columns=['Disease', 'Study', 'Taxon', 'Dir_Unadjusted', 'Dir_Adjusted'])

print("==================================================")
print("PART 3: CALCULATING ASSOCIATION SCORES")
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

final_scores = {}
for d in final_diseases:
    d_data = res_df[res_df['Disease'] == d]
    
    unadj_matrix = d_data.pivot(index='Taxon', columns='Study', values='Dir_Unadjusted').fillna(0)
    adj_matrix = d_data.pivot(index='Taxon', columns='Study', values='Dir_Adjusted').fillna(0)
    
    final_scores[f"{d}_Unadjusted"] = calc_as(unadj_matrix)
    final_scores[f"{d}_Adjusted"] = calc_as(adj_matrix)

df_combined_scores = pd.DataFrame(final_scores).reindex(SpeciesScores.index.tolist()).fillna(0)
df_combined_scores.to_csv("Comparison_Unadjusted_vs_AgeAdjusted_AS.csv")

print("==================================================")
print("PART 4: GENERATING PLOTS WITH STATS")
print("==================================================")

output_dir = "Adjustment_Comparison_Plots"
os.makedirs(output_dir, exist_ok=True)

for disease in final_diseases.keys():
    plt.figure(figsize=(8, 8)) # Square plot is better for X vs Y comparisons
    
    x = df_combined_scores[f"{disease}_Unadjusted"]
    y = df_combined_scores[f"{disease}_Adjusted"]
    
    # --- NEW: Calculate Correlation and P-value ---
    corr, p_val = pearsonr(x, y)
    
    # Base styling
    plt.scatter(x, y, color='royalblue', s=100, edgecolors='white', alpha=0.8)
    
    # Draw X=0 and Y=0 axes
    plt.axhline(0, color='black', linewidth=1.5)
    plt.axvline(0, color='black', linewidth=1.5)
    
    # Draw the Y=X diagonal line (Red Dashed)
    plt.plot([-1, 1], [-1, 1], color='red', linestyle='--', linewidth=2, label="y = x (No change after adjustment)")
    
    # --- NEW: Add the Stats Text Box to the Plot ---
    # Scientific notation for very small p-values looks much cleaner on plots
    stats_text = f"Pearson r = {corr:.3f}\nP-value = {p_val:.2e}"
    # Placing it in the bottom right corner so it doesn't block the legend in the top left
    plt.text(0.65, 0.05, stats_text, transform=plt.gca().transAxes, fontsize=12, fontweight='bold',
             bbox=dict(facecolor='white', alpha=0.9, edgecolor='black', boxstyle='round,pad=0.5'))
    
    # Format Plot
    plt.xlim(-1.1, 1.1)
    plt.ylim(-1.1, 1.1)
    plt.xlabel(f"{disease} Association Score (Unadjusted)", fontsize=14, fontweight='bold')
    plt.ylabel(f"{disease} Association Score (Age Adjusted)", fontsize=14, fontweight='bold')
    plt.title(f"Impact of Age Adjustment: {disease}", fontsize=16, fontweight='bold')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(loc='upper left')
    plt.tight_layout()
    
    file_path = os.path.join(output_dir, f"compare_adjustment_{disease}.png")
    plt.savefig(file_path, dpi=300)
    plt.close()
    print(f"Saved plot: {file_path} (r={corr:.3f}, p={p_val:.2e})")

print("==================================================")
print(f"SUCCESS! Finished in {round(time.time() - start_time, 2)} seconds.")