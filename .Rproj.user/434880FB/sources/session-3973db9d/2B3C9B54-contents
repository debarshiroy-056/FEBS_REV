# -----------------------------------------------------------------
# MASTER SCRIPT: Run Full Analysis (New Age Constraints)
# -----------------------------------------------------------------
# This script runs the complete 6-part pipeline.
#
# IT APPLIES YOUR NEW STUDY-LEVEL FILTER:
# 1. Finds studies where MINIMUM age (from 41-80 cohort) is 41-45
# 2. AND MAXIMUM age (from 41-80 cohort) is 75-80
#
# It then takes this list of eligible studies and
# re-runs the entire analysis pipeline.
# -----------------------------------------------------------------

# --- 1. Load Required Libraries ---
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}
library(dplyr)
library(tidyr)

load("Required_Data_FEBS_Pipepline.RData")

# --- 2. Check for All Required Objects ---
print("--- Checking for all required data objects...")
if (!exists("AllCombinedMetadata")) {
  stop("Error: 'AllCombinedMetadata' is not loaded.")
}
if (!exists("AllCombinedSpProfile")) {
  stop("Error: 'AllCombinedSpProfile' is not loaded.")
}
if (!exists("SpeciesScores_20241021")) {
  stop("Error: 'SpeciesScores_20241021' not found. Please load it.")
}
if (!exists("diseaseSampleList")) {
  stop("Error: 'diseaseSampleList' is not loaded.")
}
if (!exists("study_disease_age_summary")) {
  stop("Error: 'study_disease_age_summary' is not loaded. Please run 'create_study_disease_summary.R' first.")
}
print("All required data objects are loaded.")


# =================================================================
# PART A: Find Eligible Studies (Your New Constraints)
# =================================================================
print("--------------------------------------------------")
print("PART A: Finding eligible studies based on new age constraints...")

# Define your constraints
MIN_AGE_LOWER <- 41
MIN_AGE_UPPER <- 45
MAX_AGE_LOWER <- 75
MAX_AGE_UPPER <- 80

# Group by study and find the TRUE min and max age for each one
# from the *already filtered* 41-80 summary
study_age_ranges <- study_disease_age_summary %>%
  group_by(study_name) %>%
  summarise(
    Overall_Min_Age = min(`Minimum age`, na.rm = TRUE),
    Overall_Max_Age = max(`Max age`, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(!is.infinite(Overall_Min_Age) & !is.infinite(Overall_Max_Age))

# Apply your new filter to find the eligible studies
eligible_studies_df <- study_age_ranges %>%
  filter(
    Overall_Min_Age >= MIN_AGE_LOWER & Overall_Min_Age <= MIN_AGE_UPPER &
      Overall_Max_Age >= MAX_AGE_LOWER & Overall_Max_Age <= MAX_AGE_UPPER
  )

# Get the final list of study names to keep
eligible_studies_list <- eligible_studies_df$study_name

print(paste("Found", length(eligible_studies_list), "studies that meet the constraints (MinAge 41-45 AND MaxAge 75-80):"))
if (length(eligible_studies_list) > 0) {
  print(eligible_studies_list)
} else {
  stop("No studies matched your new age criteria. Stopping script.")
}


# =================================================================
# PART B: Prepare Filtered Data (Modified Step 2)
# =================================================================
print("--------------------------------------------------")
print("PART B: Preparing cohorts using *only* these eligible studies...")

# --- 3. Define Taxonomic Renaming Map ---
old_names <- c("Ruminococcus_obeum", "Eubacterium_biforme", "Clostridium_bartlettii", "Clostridium_hathewayi", "Clostridium_ramosum")
new_names <- c("Blautia_obeum", "Holdemanella_biformis", "Intestinibacter_bartlettii", "Hungatella_hathewayi", "Erysipelatoclostridium_ramosum")
renaming_map <- setNames(new_names, old_names)
profile_taxa_names <- colnames(AllCombinedSpProfile)
taxa_to_rename <- intersect(old_names, profile_taxa_names)

if (length(taxa_to_rename) > 0) {
  for (old_name in taxa_to_rename) {
    new_name <- renaming_map[old_name]
    colnames(AllCombinedSpProfile)[which(colnames(AllCombinedSpProfile) == old_name)] <- new_name
  }
}

# --- 4. Filter Species Profile (for 201 taxa) ---
all_201_taxa_names <- rownames(SpeciesScores_20241021)
profile_taxa_names_updated <- colnames(AllCombinedSpProfile)
common_taxa_to_use <- intersect(all_201_taxa_names, profile_taxa_names_updated)
filtered_species_profile <- AllCombinedSpProfile[, common_taxa_to_use, drop = FALSE]

# --- 5. Select, Merge, and Filter Metadata ---
cols_to_keep <- c("age", "study_name", "diseaseCat")
selected_metadata <- AllCombinedMetadata[, cols_to_keep]
master_data <- merge(selected_metadata, filtered_species_profile, by = "row.names")
rownames(master_data) <- master_data$Row.names
master_data$Row.names <- NULL

# --- 6. Filter by NEW Study List AND Age (41-80) ---
master_data$age_num <- as.numeric(master_data$age)

# Filter 1: By the eligible studies list
master_data_filtered_studies <- subset(master_data, study_name %in% eligible_studies_list)
# Filter 2: By the 41-80 age range
master_data_filtered_age <- subset(master_data_filtered_studies, age_num > 40 & age_num <= 80)

print(paste("After filtering for eligible studies AND age 41-80,", nrow(master_data_filtered_age), "samples remain."))

# --- 7. Split into Control and Disease Cohorts ---
target_disease_names <- names(diseaseSampleList)

control_cohort_new <- subset(master_data_filtered_age, diseaseCat == "control")
control_cohort_new$diseaseCat <- NULL

disease_cohort_new <- subset(master_data_filtered_age,
                             diseaseCat %in% target_disease_names | grepl(";", diseaseCat)
)
colnames(disease_cohort_new)[colnames(disease_cohort_new) == 'diseaseCat'] <- 'disease_name'

# --- 8. Un-nest Multi-Disease Rows ---
disease_cohort_new$SampleID <- rownames(disease_cohort_new)
disease_cohort_clean <- disease_cohort_new %>%
  separate_rows(disease_name, sep = ";")
disease_cohort_clean <- subset(disease_cohort_clean,
                               disease_name %in% target_disease_names
)
print(paste("Created 'control_cohort_new' with", nrow(control_cohort_new), "samples."))
print(paste("Created 'disease_cohort_clean' with", nrow(disease_cohort_clean), "sample-disease pairs."))




# -----------------------------------------------------------------
# Step 2.6: Check for Common Studies in NEW Filtered Data
# -----------------------------------------------------------------
# This script checks for the intersection of studies between your
# new, clean, age-filtered control and disease data frames.
#
# INPUTS:
# - 'control_cohort_new'
# - 'disease_cohort_clean'
#
# OUTPUT:
# - A printout of the common studies.
# -----------------------------------------------------------------

print("--------------------------------------------------")
print("Checking for common studies in new 41-80 data...")

# --- 1. Check for required data ---
if (!exists("control_cohort_new")) {
  stop("Error: 'control_cohort_new' data frame is missing. Please run Step 2 (v8) first.")
}
if (!exists("disease_cohort_clean")) {
  stop("Error: 'disease_cohort_clean' data frame is missing. Please run Step 2 (v8) first.")
}

# --- 2. Get Unique Study Lists ---

# Get unique studies from the control data
filtered_control_studies <- unique(control_cohort_new$study_name)
print(paste("Found", length(filtered_control_studies), "unique studies in 'control_cohort_new'."))

# Get unique studies from the disease data
filtered_disease_studies <- unique(disease_cohort_clean$study_name)
print(paste("Found", length(filtered_disease_studies), "unique studies in 'disease_cohort_clean'."))


# --- 3. Find the Intersection ---
common_studies_new <- intersect(filtered_control_studies, filtered_disease_studies)


# --- 4. Print the Results ---
print("--------------------------------------------------")
print(paste("Found", length(common_studies_new), "common studies present in BOTH new cohorts:"))
print(common_studies_new)
print("--------------------------------------------------")
print("This is the final list of studies that have both control and disease samples in the 41-80 age range.")





# =================================================================
# PART C: Structure Data (Old Steps 3 & 4)
# =================================================================
print("--------------------------------------------------")
print("PART C: Structuring data into lists...")

# --- Create Control List (by study) ---
unique_control_studies <- unique(control_cohort_new$study_name)
control_cohort_list_new <- list()
for (study in unique_control_studies) {
  control_cohort_list_new[[study]] <- control_cohort_new[control_cohort_new$study_name == study, ]
}
print(paste("Created 'control_cohort_list_new' with", length(control_cohort_list_new), "studies."))

# --- Create Disease List (nested) ---
unique_disease_names <- unique(disease_cohort_clean$disease_name)
disease_cohorts_nested_list_new <- list()
for (disease in unique_disease_names) {
  disease_df <- disease_cohort_clean[disease_cohort_clean$disease_name == disease, ]
  studies_for_this_disease <- unique(disease_df$study_name)
  inner_study_list <- list()
  for (study in studies_for_this_disease) {
    study_df <- disease_df[disease_df$study_name == study, ]
    rownames(study_df) <- study_df$SampleID
    inner_study_list[[study]] <- study_df
  }
  disease_cohorts_nested_list_new[[disease]] <- inner_study_list
}
print(paste("Created 'disease_cohorts_nested_list_new' with", length(disease_cohorts_nested_list_new), "diseases."))


# =================================================================
# PART D: Filter Diseases for >= 3 Studies (Old Step 4.5)
# =================================================================
print("--------------------------------------------------")
print("PART D: Filtering for diseases with >= 3 common studies...")

MIN_STUDIES_FINAL <- 3 # Using >= 3 studies as per our last successful pipeline
all_control_studies <- names(control_cohort_list_new)
all_disease_names <- names(disease_cohorts_nested_list_new)
disease_cohorts_final_list <- list()
diseases_kept <- c()
diseases_skipped <- c()

for (disease in all_disease_names) {
  studies_for_this_disease <- names(disease_cohorts_nested_list_new[[disease]])
  common_studies <- intersect(studies_for_this_disease, all_control_studies)
  common_study_count <- length(common_studies)
  
  if (common_study_count >= MIN_STUDIES_FINAL) {
    diseases_kept <- c(diseases_kept, disease)
    studies_to_keep <- disease_cohorts_nested_list_new[[disease]][common_studies]
    disease_cohorts_final_list[[disease]] <- studies_to_keep
  } else {
    diseases_skipped <- c(diseases_skipped, disease)
  }
}
print(paste("KEPT", length(diseases_kept), "diseases:", paste(diseases_kept, collapse=", ")))
print(paste("SKIPPED", length(diseases_skipped), "diseases."))


# =================================================================
# PART E: Create Direction Matrices (Old Step 5/6)
# =================================================================
print("--------------------------------------------------")
print("PART E: Creating All Disease Direction Matrices...")

# This will be the list of direction matrices
all_disease_direction_matrices_NEW <- list()

for (disease_name in names(disease_cohorts_final_list)) {
  print(paste("... Processing Disease:", disease_name))
  study_names_for_disease <- names(disease_cohorts_final_list[[disease_name]])
  common_studies <- intersect(study_names_for_disease, all_control_studies)
  disease_results_list <- list()
  
  for (study_name in common_studies) {
    disease_df <- disease_cohorts_final_list[[disease_name]][[study_name]]
    control_df <- control_cohort_list_new[[study_name]]
    
    if (nrow(disease_df) == 0 || nrow(control_df) == 0) { next }
    
    study_score_vec <- setNames(rep(0, length(all_201_taxa_names)), all_201_taxa_names)
    taxa_cols_in_data <- intersect(all_201_taxa_names, colnames(disease_df))
    
    for (taxon in taxa_cols_in_data) {
      x_disease <- disease_df[[taxon]]
      y_control <- control_df[[taxon]]
      p_val <- NA
      median_disease <- NA
      median_control <- NA
      
      tryCatch({
        if (length(unique(x_disease)) > 1 || length(unique(y_control)) > 1) {
          p_val <- wilcox.test(x_disease, y_control)$p.value
        }
        median_disease <- median(x_disease, na.rm = TRUE)
        median_control <- median(y_control, na.rm = TRUE)
      }, error = function(e) {})
      
      if (!is.na(p_val) && !is.na(median_disease) && !is.na(median_control)) {
        if (p_val <= 0.1) {
          study_score_vec[taxon] <- sign(median_disease - median_control)
        }
      }
    }
    disease_results_list[[study_name]] <- study_score_vec
  }
  
  if (length(disease_results_list) > 0) {
    disease_matrix_df <- as.data.frame(do.call(cbind, disease_results_list))
    disease_matrix_df <- disease_matrix_df[all_201_taxa_names, ]
    disease_matrix_df[is.na(disease_matrix_df)] <- 0
    all_disease_direction_matrices_NEW[[disease_name]] <- disease_matrix_df
  }
}
print("SUCCESS: Part E complete. 'all_disease_direction_matrices_NEW' is ready.")


# =================================================================
# PART F: Calculate Association Scores (Old Step 7)
# =================================================================
print("--------------------------------------------------")
print("PART F: Calculating Final Association Scores...")

all_disease_association_scores_NEW <- list()

for (disease_name in names(all_disease_direction_matrices_NEW)) {
  
  direction_matrix <- all_disease_direction_matrices_NEW[[disease_name]]
  if (ncol(direction_matrix) == 0) { next }
  
  SP_vec <- apply(direction_matrix, 1, function(row) sum(row == 1, na.rm = TRUE))
  SN_vec <- apply(direction_matrix, 1, function(row) sum(row == -1, na.rm = TRUE))
  T_val <- ncol(direction_matrix)
  
  pattern_df <- data.frame(
    SP = SP_vec[all_201_taxa_names],
    SN = SN_vec[all_201_taxa_names],
    T = T_val,
    row.names = all_201_taxa_names
  )
  pattern_df[is.na(pattern_df)] <- 0
  
  pmax_val <- pmax(pattern_df$SP, pattern_df$SN) + 0.00001
  pmin_val <- pmin(pattern_df$SP, pattern_df$SN) + 0.00001
  
  consensus_score <- (1 - (pmin_val / pmax_val))
  consensus_score[pmax_val == 0.00001] <- 1 
  
  association_score <- ((pattern_df$SP - pattern_df$SN) / pattern_df$T) * consensus_score
  
  association_score_df <- data.frame(
    SP = pattern_df$SP,
    SN = pattern_df$SN,
    T = T_val,
    AS = association_score,
    row.names = all_201_taxa_names
  )
  
  all_disease_association_scores_NEW[[disease_name]] <- association_score_df
}
print("SUCCESS: Part F complete. 'all_disease_association_scores_NEW' is ready.")


# =================================================================
# PART G: Create Final Summary Matrix (Old Step 8)
# =================================================================
print("--------------------------------------------------")
print("PART G: Combining all Association Scores into a final matrix...")

# Extract the 'AS' columns
as_scores_list <- lapply(all_disease_association_scores_NEW, function(df) {
  as_vector <- df$AS
  names(as_vector) <- rownames(df)
  return(as_vector)
})

final_association_matrix_NEW_CONSTRAINTS <- as.data.frame(do.call(cbind, as_scores_list))
final_association_matrix_NEW_CONSTRAINTS <- final_association_matrix_NEW_CONSTRAINTS[all_201_taxa_names, ]
final_association_matrix_NEW_CONSTRAINTS[is.na(final_association_matrix_NEW_CONSTRAINTS)] <- 0

# --- FINAL CONFIRMATION ---
print("--------------------------------------------------")
print("SUCCESS: PIPELINE COMPLETE.")
print("Created the final 'final_association_matrix_NEW_CONSTRAINTS'.")
print(paste("Dimensions:",
            nrow(final_association_matrix_NEW_CONSTRAINTS), "taxa (rows) x",
            ncol(final_association_matrix_NEW_CONSTRAINTS), "diseases (cols)"))
print("--------------------------------------------------")
print("Preview of the final matrix (first 6 taxa, first 5 diseases):")
if (ncol(final_association_matrix_NEW_CONSTRAINTS) > 0) {
  print(head(final_association_matrix_NEW_CONSTRAINTS[, 1:min(5, ncol(final_association_matrix_NEW_CONSTRAINTS))]))
} else {
  print("No diseases had >= 2 studies, so the final matrix is empty.")
}




# -----------------------------------------------------------------
# Manual Check for "T2D" (New Data, Age 41-80 Constraints)
# -----------------------------------------------------------------
# This script performs the complete AS calculation for 'T2D'
# using your supervisor's median-based formula on the NEW,
# age-constrained data.
#
# INPUTS:
# - 'disease_cohorts_final_list' (from Part D of your master script)
# - 'control_cohort_list_new' (from Part C of your master script)
# - 'SpeciesScores_20241021'
#
# OUTPUTS:
# 1. 'T2D_direction_matrix_NEW'
# 2. 'T2D_association_scores_df_NEW'
# -----------------------------------------------------------------

# --- A. Setup ---

print("--------------------------------------------------")
print("Running Manual Check for 'T2D' (New Age-Constrained Data)...")

# Define the target disease
target_disease <- "T2D"

# Check for required data
if (!exists("disease_cohorts_final_list")) {
  stop("Error: 'disease_cohorts_final_list' is missing. Please run the master script (Part D) first.")
}
if (!exists("control_cohort_list_new")) {
  stop("Error: 'control_cohort_list_new' is missing. Please run the master script (Part C) first.")
}
if (!exists("SpeciesScores_20241021")) {
  stop("Error: 'SpeciesScores_20241021' not found. Please load it first.")
}
if (!target_disease %in% names(disease_cohorts_final_list)) {
  stop(paste("Error: Target disease", target_disease, "not found in the filtered disease list 'disease_cohorts_final_list'."))
}

# Get the list of all 201 taxa
all_201_taxa_names <- rownames(SpeciesScores_20241021)
# Get the list of taxa that are *actually* in the new profile
common_taxa_to_use <- intersect(all_201_taxa_names, colnames(control_cohort_list_new[[1]]))


# --- B. Part 1: Create the Direction Matrix for T2D ---
print(paste("--- Part 1: Creating Direction Matrix for:", target_disease, "---"))

# Get the list of studies for this disease
study_names_for_disease <- names(disease_cohorts_final_list[[target_disease]])

# Initialize a list to store the results vector for each study
T2D_results_list <- list()

print(paste("Found", length(study_names_for_disease), "common studies for", target_disease, ". Processing..."))

# Loop through each study
for (study_name in study_names_for_disease) {
  
  # --- Get Data ---
  disease_df <- disease_cohorts_final_list[[target_disease]][[study_name]]
  control_df <- control_cohort_list_new[[study_name]]
  
  # Check for empty data frames
  if (nrow(disease_df) == 0 || nrow(control_df) == 0) {
    print(paste("     - Skipping", study_name, "due to 0 samples."))
    next
  }
  
  # --- Initialize Score Vector ---
  study_score_vec <- rep(0, length(all_201_taxa_names))
  names(study_score_vec) <- all_201_taxa_names
  
  # --- Loop through each *available* taxon ---
  # (We only loop through the 196, the other 5 will stay 0)
  for (taxon in common_taxa_to_use) {
    
    x_disease <- disease_df[[taxon]]
    y_control <- control_df[[taxon]]
    
    p_val <- NA
    median_disease <- NA
    median_control <- NA
    
    # Run the test safely
    tryCatch({
      if (length(unique(x_disease)) > 1 || length(unique(y_control)) > 1) {
        p_val <- wilcox.test(x_disease, y_control)$p.value
      }
      median_disease <- median(x_disease, na.rm = TRUE)
      median_control <- median(y_control, na.rm = TRUE)
    }, error = function(e) {
      # On error, p_val remains NA
    })
    
    # --- Apply Supervisor's Formula ---
    if (!is.na(p_val) && !is.na(median_disease) && !is.na(median_control)) {
      if (p_val <= 0.1) {
        median_diff <- median_disease - median_control
        study_score_vec[taxon] <- sign(median_diff)
      }
    }
  } # End taxon loop
  
  # --- Store Results for this study ---
  T2D_results_list[[study_name]] <- study_score_vec
  
} # End study loop

# --- Combine into the final matrix ---
T2D_direction_matrix_NEW <- as.data.frame(do.call(cbind, T2D_results_list))

# Ensure it has all 201 taxa, with NAs for the 5 missing
T2D_direction_matrix_NEW <- T2D_direction_matrix_NEW[all_201_taxa_names, ]
# Replace all NAs (from missing taxa) with 0
T2D_direction_matrix_NEW[is.na(T2D_direction_matrix_NEW)] <- 0

print("... Part 1 Complete. 'T2D_direction_matrix_NEW' is ready.")
print(paste("Dimensions:", nrow(T2D_direction_matrix_NEW), "taxa x", ncol(T2D_direction_matrix_NEW), "studies"))


# --- C. Part 2: Calculate Association Score for T2D ---
print(paste("--- Part 2: Calculating Association Score for:", target_disease, "---"))

# 1. This is your 'df_T2D_pattern' command
# SP (Increased) = count of +1
# SN (Decreased) = count of -1
pattern_df <- data.frame(
  "Increased" = apply(T2D_direction_matrix_NEW, 1, function(x) length(x[x == 1])),
  "Decreased" = apply(T2D_direction_matrix_NEW, 1, function(x) length(x[x == -1]))
)

# 2. This is your 'apply(...)' command
# T = total number of studies
T_val <- ncol(T2D_direction_matrix_NEW)

as_vector <- apply(pattern_df, 1, function(x) {
  SP <- x[1] # Increased
  SN <- x[2] # Decreased
  
  # Calculate consensus score
  pmax_val <- max(SP, SN) + 0.00001
  pmin_val <- min(SP, SN) + 0.00001
  
  consensus_score <- (1 - (pmin_val / pmax_val))
  if (pmax_val == 0.00001) { consensus_score <- 1 } # Handle SP=0, SN=0
  
  # Calculate final AS
  association_score <- ((SP - SN) / T_val) * consensus_score
  return(association_score)
})

# 3. Create the final results data frame
T2D_association_scores_df_NEW <- data.frame(
  SP = pattern_df$Increased,
  SN = pattern_df$Decreased,
  T = T_val,
  AS = as_vector,
  row.names = all_201_taxa_names # Ensure row names are set
)

# --- D. Confirmation ---
print("--------------------------------------------------")
print("SUCCESS: Manual check for T2D complete.")
print("Created 'T2D_direction_matrix_NEW' and 'T2D_association_scores_df_NEW'.")
print("--------------------------------------------------")
print("Preview of final Association Scores for T2D:")
print(head(T2D_association_scores_df_NEW))



### Disease-wise Association Score Plots ###########

library(ggplot2)

# --- Combine into a data frame ---
df_combined_scores <- data.frame(
  age_association_controls = df_combined_scores$age_association_controls,
  T2D = final_association_matrix_NEW_CONSTRAINTS[,"T2D"]
)

# --- Plot ---
ggplot(df_combined_scores, aes(x = age_association_controls, y = T2D)) +
  geom_point(
    size = 4,
    colour = ifelse(
      (df_combined_scores$age_association_controls > 0 & df_combined_scores$T2D > 0), "slateblue1",
      ifelse(
        (df_combined_scores$age_association_controls < 0 & df_combined_scores$T2D < 0), "maroon2",
        ifelse(
          sign(df_combined_scores$age_association_controls) != sign(df_combined_scores$T2D),
          "lightgoldenrod3",
          "grey"
        )
      )
    )
  ) +
  geom_hline(yintercept = 0, linewidth = 2) +
  geom_vline(xintercept = 0, linewidth = 2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15)
  ) +
  xlab("Age Association Score (Controls)") +
  ylab("T2D Association Score")


# Load necessary libraries
library(ggplot2)
library(dplyr)

# 1. Disease list
disease_list <- c("CRC", "CVD", "NAFLD_FLD", "Parkinsons", "Polyps", "T2D")

# 2. Output folder (no legends)
dir.create("individual_disease_plots_NO_LEGEND", showWarnings = FALSE)

# 3. Loop through diseases
for (disease_name in disease_list) {
  
  # Dynamic color groups
  df_temp_with_groups <- df_combined_scores %>%
    mutate(color_group = case_when(
      age_association_controls > 0 & .data[[disease_name]] > 0 ~ "Positive Concordance",
      age_association_controls < 0 & .data[[disease_name]] < 0 ~ "Negative Concordance",
      sign(age_association_controls) != sign(.data[[disease_name]]) ~ "Discordant",
      TRUE ~ "On Axis"
    ))
  
  # Plot without legend
  p <- ggplot(df_temp_with_groups,
              aes(x = age_association_controls,
                  y = .data[[disease_name]],
                  colour = color_group)) +
    geom_point(size = 4) +
    scale_color_manual(values = c(
      "Positive Concordance" = "slateblue1",
      "Negative Concordance" = "maroon2",
      "Discordant"           = "lightgoldenrod3",
      "On Axis"              = "grey"
    )) +
    geom_hline(yintercept = 0, linewidth = 2) +
    geom_vline(xintercept = 0, linewidth = 2) +
    labs(
      x = "Age Association Score (Controls)",
      y = paste("Association Score (", disease_name, ")")
    ) +
    theme_bw(base_size = 14, base_family = "Arial") +
    theme(
      legend.position = "none"   # <- legend removed completely
    )
  
  # Save to new folder
  file_path <- paste0("individual_disease_plots_NO_LEGEND/plot_age_vs_", disease_name, ".png")
  ggsave(file_path, plot = p, width = 8, height = 6)
}

###############################################################


# ==============================================================================
# Script Name: Final Disease Study Summary Generator
# Description: 
#   This script aggregates sample counts and study information for the final 
#   meta-analysis dataset. It cross-references disease cohorts with control 
#   cohorts to ensure accurate sample counts for the 41-80 age group.
#
# Inputs:
#   1. disease_cohorts_nested_list_new: Nested list containing disease dataframes per study.
#   2. control_cohort_list_new: List containing control dataframes per study.
#
# Output:
#   - A dataframe 'final_disease_study_summary' sorted by study count.
#   - An Excel file 'final_disease_study_summary.xlsx'.
# ==============================================================================

library(dplyr)
library(writexl)

# ------------------------------------------------------------------------------
# 1. Verification Step
# ------------------------------------------------------------------------------
# Ensure the necessary data structures from previous pipeline steps exist 
# in the global environment before proceeding.
if (!exists("disease_cohorts_nested_list_new")) {
  stop("Error: 'disease_cohorts_nested_list_new' not found. Please run the data preparation step first.")
}
if (!exists("control_cohort_list_new")) {
  stop("Error: 'control_cohort_list_new' not found. Please run the control cohort generation step first.")
}

# Initialize an empty list to store the summary row for each disease
summary_list <- list()

# ------------------------------------------------------------------------------
# 2. Iterative Summary Generation
# ------------------------------------------------------------------------------
# Loop through every disease present in the nested list
for (disease_name in names(disease_cohorts_nested_list_new)) {
  
  # Retrieve the list of studies associated with the current disease
  study_list <- disease_cohorts_nested_list_new[[disease_name]]
  study_names <- names(study_list)
  
  # -- A. Calculate Study Metrics --
  num_studies <- length(study_names)
  
  # Create a semicolon-separated string of study names for the report (e.g., "StudyA;StudyB")
  study_names_str <- paste(sort(study_names), collapse = ";")
  
  # -- B. Count Disease Samples --
  # Sum the number of rows (samples) in every study dataframe for this disease
  num_disease_samples <- sum(sapply(study_list, nrow))
  
  # -- C. Count Control Samples --
  # We must only count controls from studies that are actually used for this disease.
  # This avoids counting controls from studies that might have been filtered out.
  num_control_samples <- 0
  for (study in study_names) {
    # Check if this specific study exists in our clean control list
    if (study %in% names(control_cohort_list_new)) {
      num_control_samples <- num_control_samples + 
        nrow(control_cohort_list_new[[study]])
    }
  }
  
  # -- D. Store Result --
  # Create a single-row dataframe for this disease and add it to the list
  summary_list[[disease_name]] <- data.frame(
    Disease = disease_name,
    `Number of Studies` = num_studies,
    `Study Names` = study_names_str,
    `Number of Disease Samples` = num_disease_samples,
    `Number of Control Samples` = num_control_samples,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# 3. Data Aggregation and Formatting
# ------------------------------------------------------------------------------
# Combine all individual disease rows into a single master dataframe
final_disease_study_summary <- bind_rows(summary_list)

# (Optional) First sort alphabetically by Disease name for tidiness
final_disease_study_summary <- final_disease_study_summary %>% 
  arrange(Disease)

# Print initial view
print("Initial Summary Table:")
print(final_disease_study_summary)

# ------------------------------------------------------------------------------
# 4. Sorting by Importance (Number of Studies)
# ------------------------------------------------------------------------------
# Dynamically find the column name containing "Studies" to handle potential naming variations
study_col <- grep("Studies", colnames(final_disease_study_summary), value = TRUE)

# Sort the dataframe in Descending order of study count (most studied diseases on top)
final_disease_study_summary <- final_disease_study_summary %>% 
  arrange(desc(.data[[study_col]]))

# Display the sorted table in the RStudio viewer
View(final_disease_study_summary)

# ------------------------------------------------------------------------------
# 5. Export to Excel
# ------------------------------------------------------------------------------
# Check if writexl is installed, install if missing (safety check)
if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl")
}

# Save the final sorted dataframe to an Excel file in the working directory
write_xlsx(final_disease_study_summary, path = "final_disease_study_summary.xlsx")

# Confirmation message
cat("SUCCESS: Summary table created and saved as 'final_disease_study_summary.xlsx'\n")


