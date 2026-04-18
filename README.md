# Medication Data Figure Generation README

This document gives the proper from-scratch pipeline order for generating both figures and their exported tables.

Figures covered:
1. Age-adjusted disease association heatmap
2. Cohort-specific species-medication binary heatmap

All commands below are workspace-tested and run from the repository root.

## 1) Initialize Environment

```bash
cd /home/debarshi/medication_data
conda activate base
```

## 2) Check Required Inputs (From Scratch)

Before running scripts, confirm these core inputs exist:

```bash
ls -1 \
  42_medications.xlsx \
  drug_matrix.csv \
  AllCombinedSpProfile.csv \
  Final_Association_Scores_Age_Adjusted.csv \
  Focal_Species_Cohort_Medication_Table_TSG.xlsx
```

Notes:
- If `AllCombinedSpProfile.csv` is empty/missing, `Run_Regressions.py` can load it from `Required_Data_FEBS_Pipeline.RData` or `Required_Data_FEBS_Pipepline.RData` when available.
- Figure B uses `Focal_Species_Cohort_Medication_Table_TSG.xlsx` directly.

## 3) Build Upstream Regression Outputs

### Script
- `Run_Regressions.py`

### Command
```bash
python Run_Regressions.py
```

### Outputs
- `All_Cohorts_Regression_Results.csv`
- `Raw_Multivariate_Regression_Results.csv`

## 4) Build Focal Cohort-Species-Medication CSV (Linked Upstream Step)

This is linked to the workflow and useful for traceability, although Figure B consumes the curated XLSX file.

### Script
- `generate_focal_species_table.py`

### Command
```bash
python generate_focal_species_table.py
```

### Output
- `Focal_Species_Cohort_Medication_Table.csv`

Optional utility: append binary indicator column to the focal CSV.

```bash
set -e
in='Focal_Species_Cohort_Medication_Table.csv'
out='Focal_Species_Cohort_Medication_Table.tmp.csv'
{
  IFS= read -r header
  printf '%s\n' "$header,Binary_Value"
  while IFS= read -r line; do
    printf '%s\n' "$line,1"
  done
} < "$in" > "$out"
mv "$out" "$in"
```

## 5) Generate Figure A (Age-Adjusted Heatmap)

### Script
- `Generate_Fig3Ex_Heatmap.py`

### Required Input
- `Final_Association_Scores_Age_Adjusted.csv`

### Command
```bash
python Generate_Fig3Ex_Heatmap.py
```

### Outputs
- `Figure3_AgeAdjusted_Positive_Only.png`
- `Heatmap_Carpet_Exports/Figure3_AgeAdjusted_Carpet_Wide.csv`
- `Heatmap_Carpet_Exports/Figure3_AgeAdjusted_Carpet_Long.csv`

## 6) Generate Figure B (Cohort-Specific Binary Heatmap)

### Script
- `Generate_Cohort_Specific_Binary_Heatmaps.R`

### Required Input
- `Focal_Species_Cohort_Medication_Table_TSG.xlsx`

### Command
```bash
Rscript Generate_Cohort_Specific_Binary_Heatmaps.R
```

### Outputs
- `Binary_Heatmap_Cohort_Outputs/Single_Cohort_Stacked_Binary_Heatmap.png`
- `Binary_Heatmap_Cohort_Outputs/Single_Cohort_Stacked_Binary_Heatmap.pdf`
- `Binary_Heatmap_Cohort_Outputs/Single_Heatmap_Binary_Carpet_Long.csv`
- `Binary_Heatmap_Cohort_Outputs/Single_Heatmap_Binary_Carpet_Wide.csv`

## 7) One-Block End-to-End Run (From Scratch)

```bash
cd /home/debarshi/medication_data
conda activate base

# Input checks
ls -1 \
  42_medications.xlsx \
  drug_matrix.csv \
  AllCombinedSpProfile.csv \
  Final_Association_Scores_Age_Adjusted.csv \
  Focal_Species_Cohort_Medication_Table_TSG.xlsx

# Upstream linked generation
python Run_Regressions.py
python generate_focal_species_table.py

# Figure generation
python Generate_Fig3Ex_Heatmap.py
Rscript Generate_Cohort_Specific_Binary_Heatmaps.R
```

## 8) Verify Final Outputs

```bash
ls -1 \
  All_Cohorts_Regression_Results.csv \
  Raw_Multivariate_Regression_Results.csv \
  Focal_Species_Cohort_Medication_Table.csv \
  Figure3_AgeAdjusted_Positive_Only.png \
  Heatmap_Carpet_Exports/Figure3_AgeAdjusted_Carpet_Wide.csv \
  Heatmap_Carpet_Exports/Figure3_AgeAdjusted_Carpet_Long.csv \
  Binary_Heatmap_Cohort_Outputs/Single_Cohort_Stacked_Binary_Heatmap.png \
  Binary_Heatmap_Cohort_Outputs/Single_Cohort_Stacked_Binary_Heatmap.pdf \
  Binary_Heatmap_Cohort_Outputs/Single_Heatmap_Binary_Carpet_Long.csv \
  Binary_Heatmap_Cohort_Outputs/Single_Heatmap_Binary_Carpet_Wide.csv
```

If all paths print, the pipeline completed successfully.

## 9) Validation Status

- `python Generate_Fig3Ex_Heatmap.py` completed successfully and produced all listed Figure A outputs.
- `Rscript Generate_Cohort_Specific_Binary_Heatmaps.R` completed successfully and produced all listed Figure B outputs.
- No script internals were modified during validation.
