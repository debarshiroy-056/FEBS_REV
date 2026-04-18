import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap  # <-- NEW IMPORT

print("==================================================")
print("PART 1: LOADING & STRICT FILTERING")
print("==================================================")

# 1. Load the age-adjusted association scores
df = pd.read_csv("Final_Association_Scores_Age_Adjusted.csv", index_col=0)

# 2. Select only the diseases shown in the figure, in the desired order
target_diseases = [
    "CRC",
    "CVD",
    "T2D",
    "Parkinsons",
    "NAFLD_FLD",
    "Polyps",
]
df_diseases = df[target_diseases]

# 3. Use only the species from the provided figure, in exact left-to-right order.
species_order = [
    ("Anaerostipes_hadrus", "Anaerostipes_hadrus"),
    ("Dorea_longicatena", "Dorea_longicatena"),
    ("Coprococcus_comes", "Coprococcus_comes"),
    ("Bacteroides_ovatus", "Bacteroides_ovatus"),
    ("Bacteroides_vulgatus", "Bacteroides_vulgatus"),
    ("Parasutterella_excrementihominis", "Parasutterella_excrementihominis"),
    ("Clostridium_innocuum", "Clostridium_innocuum"),
    ("Streptococcus_anginosus", "Streptococcus_anginosus_group"),
    ("Streptococcus_oralis", "Streptococcus_oralis"),
    ("Veillonella_parvula", "Veillonella_parvula"),
    ("Eisenbergiella_tayi", "Eisenbergiella_tayi"),
    ("Bifidobacterium_dentium", "Bifidobacterium_dentium"),
    ("Streptococcus_mutans", "Streptococcus_mutans"),
    ("Streptococcus_vestibularis", "Streptococcus_vestibularis"),
    ("Streptococcus_gordonii", "Streptococcus_gordonii"),
    ("Klebsiella_variicola", "Klebsiella_variicola"),
    ("Rothia_mucilaginosa", "Rothia_mucilaginosa"),
    ("Ruminococcaceae_bacterium_D5", "Ruminococcaceae_bacterium_D5"),
    ("Coprobacillus_cateniformis", "Coprobacillus_cateniformis"),
    ("Veillonella_sp_T11011_6", "Veillonella_sp_T11011_6"),
    ("Veillonella_atypica", "Veillonella_atypica"),
]

display_labels = [display_name for display_name, _ in species_order]
source_labels = [source_name for _, source_name in species_order]

missing_species = [sp for sp in source_labels if sp not in df_diseases.index]
if missing_species:
    raise ValueError(
        "These required species were not found in Final_Association_Scores_Age_Adjusted.csv: "
        + ", ".join(missing_species)
    )

df_filtered = df_diseases.loc[source_labels].copy()
df_filtered.index = display_labels

# Remove the exact age-negative taxa shown in Plot A.
# The reference figure marks these six taxa as the orange age-negative group.
removed_taxa = [
    "Anaerostipes_hadrus",
    "Dorea_longicatena",
    "Coprococcus_comes",
    "Bacteroides_ovatus",
    "Bacteroides_vulgatus",
    "Parasutterella_excrementihominis",
]
removed_display_labels = [display_name for display_name, source_name in species_order if source_name in removed_taxa]
kept_labels = [label for label in display_labels if label not in removed_display_labels]
df_filtered = df_filtered.loc[kept_labels]

print(f"Removed {len(removed_display_labels)} age-negative taxa from Plot A.")
print("Removed taxa:")
for taxon in removed_display_labels:
    print(f"  - {taxon}")
print(f"Plotting {len(kept_labels)} remaining taxa in exact figure order.")

# 4. Transpose so diseases are rows and species are columns
df_heatmap = df_filtered.T

# Export carpet tables that mirror the plotted matrix.
carpet_export_dir = "Heatmap_Carpet_Exports"
carpet_wide = df_heatmap.copy()
carpet_wide.insert(0, "Display_Row", carpet_wide.index)

carpet_long = (
    df_heatmap.reset_index()
    .rename(columns={"index": "Display_Row"})
    .melt(id_vars="Display_Row", var_name="Taxon", value_name="Association_Score")
)

carpet_wide_path = f"{carpet_export_dir}/Figure3_AgeAdjusted_Carpet_Wide.csv"
carpet_long_path = f"{carpet_export_dir}/Figure3_AgeAdjusted_Carpet_Long.csv"

if not pd.io.common.file_exists(carpet_export_dir):
    import os
    os.makedirs(carpet_export_dir, exist_ok=True)

carpet_wide.to_csv(carpet_wide_path, index=False)
carpet_long.to_csv(carpet_long_path, index=False)

print(f"Exported carpet tables to '{carpet_export_dir}'.")

print("==================================================")
print("PART 2: DRAWING THE STRICT HEATMAP")
print("==================================================")

fig, (cax, ax) = plt.subplots(
    ncols=2,
    figsize=(16, 9),
    gridspec_kw={"width_ratios": [1.3, 28], "wspace": 0.12},
)

# --- NEW: Define custom colormap to match the exact image colors ---
# Use deeper endpoints so the legend and cells match the reference figure more closely.
custom_cmap = LinearSegmentedColormap.from_list(
    "custom_pink_green",
    [
        (0.0, "#8E005C"),
        (0.35, "#C01C76"),
        (0.50, "#F7F7F7"),
        (0.65, "#5FA230"),
        (1.0, "#2F6F16"),
    ],
)

ax = sns.heatmap(
    df_heatmap,
    cmap=custom_cmap,         # <-- APPLIED CUSTOM COLORMAP HERE
    vmin=-1.0,
    vmax=1.0,
    linewidths=1.5,           # Increased slightly for visual weight
    linecolor="black",        # <-- CHANGED TO BLACK TO MATCH IMAGE GRID
    cbar_ax=cax,
    cbar_kws={"label": "Association Score (AS)"},
)

cax.yaxis.set_ticks_position("left")
cax.yaxis.set_label_position("left")
cax.tick_params(labelsize=10)

# --- NEW: Make the taxa (X-axis labels) bolded ---
ax.set_xticklabels(
    ax.get_xticklabels(),
    rotation=60,
    ha="right",
    rotation_mode="anchor",
    fontsize=11,          
    fontstyle="italic",   
    fontweight="bold",    
    color="black",
)

ax.set_yticklabels(
    ax.get_yticklabels(),
    rotation=0,
    fontsize=11,
    fontweight="bold",
    color="black",
)
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
ax.tick_params(axis="y", labelright=True, labelleft=False, pad=4)

ax.set_title("Age-Adjusted Disease Association Scores", fontsize=15, fontweight="bold", pad=20)
plt.xlabel("")
plt.ylabel("")

fig.subplots_adjust(left=0.10, right=0.88, bottom=0.34, top=0.86, wspace=0.16)
output_filename = "Figure3_AgeAdjusted_Positive_Only.png"
plt.savefig(output_filename, dpi=300, bbox_inches="tight", pad_inches=0.25)

print("==================================================")
print(f"SUCCESS! Strict-ordered heatmap saved as '{output_filename}'")