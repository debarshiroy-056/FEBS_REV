library(readxl)


load("ElderMet_MetaCardis_NuAge_Drug_Data.RData")
load("Required_Data_FEBS_Pipepline.RData")


# Read first Excel file
medications_df <- read_excel("42_medications.xlsx")

# Read second Excel file
metadata_df <- read_excel("heitz_bushartA_2016_metadata.xlsx")

str(medications_df)
str(metadata_df)
str(ElderMet_Drug_data)
str(MetaCardis_drug_data)
str(NuAge_Drug_data)
str(AllCombinedMetadata)
str(AllCombinedSpProfile)
dir()


############################################################
# STEP 1 — Load required libraries
############################################################

library(dplyr)
library(stringr)
library(readxl)


############################################################
# STEP 2 — Identify species belonging to the target genera
############################################################
# Target genera from professor:
# Streptococcus
# Veillonella
# Rothia

species_names <- colnames(AllCombinedSpProfile)

target_species <- species_names[
  str_detect(species_names,
             "^Streptococcus_|^Veillonella_|^Rothia_")
]

length(target_species)
head(target_species)


############################################################
# STEP 3 — Extract abundance of selected species
############################################################

selected_species_df <- AllCombinedSpProfile[, target_species]

dim(selected_species_df)


############################################################
# STEP 4 — Calculate genus-level abundance
############################################################

# Streptococcus abundance
streptococcus_abundance <- rowSums(
  selected_species_df[, str_detect(target_species, "^Streptococcus_")],
  na.rm = TRUE
)

# Veillonella abundance
veillonella_abundance <- rowSums(
  selected_species_df[, str_detect(target_species, "^Veillonella_")],
  na.rm = TRUE
)

# Rothia abundance
rothia_abundance <- rowSums(
  selected_species_df[, str_detect(target_species, "^Rothia_")],
  na.rm = TRUE
)


############################################################
# STEP 5 — Create microbial feature matrix
############################################################

microbial_matrix <- data.frame(
  
  Sample = AllCombinedMetadata$sample_id,
  
  Streptococcus = streptococcus_abundance,
  Veillonella   = veillonella_abundance,
  Rothia        = rothia_abundance
  
)

head(microbial_matrix)


############################################################
# STEP 6 — Combine all medication datasets
############################################################

# Combine the three cohorts

all_drug_data <- bind_rows(
  ElderMet_Drug_data,
  MetaCardis_drug_data,
  NuAge_Drug_data
)

dim(all_drug_data)


############################################################
# STEP 7 — Convert medication counts to binary
############################################################
# One-hot encoding (presence / absence)

drug_matrix <- all_drug_data %>%
  mutate(
    across(
      -Sample,
      ~ ifelse(. > 0, 1, 0)
    )
  )


############################################################
# STEP 8 — Merge microbiome and medication data
############################################################

final_matrix <- microbial_matrix %>%
  inner_join(
    drug_matrix,
    by = c("Sample" = "Sample")
  )

dim(final_matrix)
head(final_matrix)


############################################################
# STEP 9 — Optional: keep only the 42 medications
############################################################

medications_df <- read_excel("42_medications.xlsx")

drug_list <- colnames(medications_df)[seq(2,85,2)]

selected_drugs <- drug_list[drug_list %in% colnames(final_matrix)]

final_matrix_42 <- final_matrix %>%
  select(
    Sample,
    Streptococcus,
    Veillonella,
    Rothia,
    all_of(selected_drugs)
  )

final_matrix_42[is.na(final_matrix_42)] <- 0
sum(is.na(final_matrix_42))
grep("^Rothia_", colnames(AllCombinedSpProfile), value = TRUE)
summary(final_matrix_42[,2:4])


############################################################
# STEP 10 — Save final dataset
############################################################

write.csv(
  final_matrix_42,
  "Microbe_Medication_Matrix.csv",
  row.names = FALSE
)




# ==============================================================================
# SCRIPT: Medication_Microbiome_Footprint_PCoA.R
# DESCRIPTION: Runs multivariate regression for target oral species against 
# 42 medications simultaneously (adjusting for polypharmacy). Translates results
# into a +2 to -2 scoring matrix, computes Jaccard distance, and plots PCoA.
# ==============================================================================

library(dplyr)
library(stringr)
library(vegan) # For Jaccard distance and PCoA
library(ggplot2)
library(ggrepel) # For nice text labels on the PCoA plot

print("==================================================")
print("STARTING MULTIVARIATE MEDICATION SCORING PIPELINE")
print("==================================================")

# --- 1. Get Individual Species Data ---
# The whiteboard shows individual ticks for species under the genera.
# We need individual species abundances, not the genus sum.

target_species <- species_names[str_detect(species_names, "^Streptococcus_|^Veillonella_|^Rothia_")]

# Extract these specific columns from the main profile
species_abundance_df <- AllCombinedSpProfile[, target_species, drop = FALSE]
species_abundance_df$Sample <- rownames(AllCombinedSpProfile)

# Merge with your previously generated binary drug matrix
# (Assuming 'drug_matrix' and 'selected_drugs' exist from your Step 7 & 9)
analysis_df <- species_abundance_df %>%
  inner_join(drug_matrix[, c("Sample", selected_drugs)], by = "Sample")

# --- 2. Initialize the Scoring Matrix ---
# Rows = Species, Columns = Medications
scoring_matrix <- matrix(0, nrow = length(target_species), ncol = length(selected_drugs))
rownames(scoring_matrix) <- target_species
colnames(scoring_matrix) <- selected_drugs


# --- 3. Multivariate Regression Loop ---
print("Running Multivariate Linear Models for each species...")

for (sp in target_species) {
  
  # Build the formula: Species ~ Drug1 + Drug2 + ... + Drug42
  formula_str <- paste0("`", sp, "` ~ `", paste(selected_drugs, collapse = "` + `"), "`")
  
  # Run the model
  model <- lm(as.formula(formula_str), data = analysis_df)
  mod_summary <- summary(model)$coefficients
  coef_names <- rownames(mod_summary) # Get the actual names R generated
  
  # --- 4. Apply the Professor's Scoring Rules ---
  for (drug in selected_drugs) {
    
    # Safely find the matching coefficient row, regardless of backticks
    # We use fixed() to prevent errors if drug names have parentheses
    matching_coef <- coef_names[str_detect(coef_names, fixed(drug))]
    
    if (length(matching_coef) > 0) {
      
      actual_name <- matching_coef[1] # Grab the matched name
      
      p_val <- mod_summary[actual_name, "Pr(>|t|)"]
      estimate <- mod_summary[actual_name, "Estimate"]
      
      # The +2 to -2 Logic
      score <- 0
      if (!is.na(p_val)) {
        if (p_val <= 0.05 && estimate > 0) {
          score <- 2
        } else if (p_val > 0.05 && p_val <= 0.1 && estimate > 0) {
          score <- 1
        } else if (p_val > 0.05 && p_val <= 0.1 && estimate < 0) {
          score <- -1
        } else if (p_val <= 0.05 && estimate < 0) {
          score <- -2
        }
      }
      scoring_matrix[sp, drug] <- score
    }
  }
}

print("Scoring Matrix Generated Successfully!")

# --- 5. Prepare Matrix for Distance Calculation ---
# Transpose so Medications are ROWS for clustering
medication_matrix <- t(scoring_matrix)

# Remove medications that have ZERO effect on ALL target species
active_meds <- rowSums(abs(medication_matrix)) > 0
medication_matrix_clean <- medication_matrix[active_meds, , drop = FALSE]

print(paste("Dropped", sum(!active_meds), "medications that showed no associations."))
print(paste("Proceeding to PCoA with", nrow(medication_matrix_clean), "active medications."))

# Failsafe: Only run PCoA if we have at least 3 active medications
if (nrow(medication_matrix_clean) >= 3) {
  
  # --- 6. Calculate Jaccard Distance & Run PCoA ---
  dist_matrix <- vegdist(abs(medication_matrix_clean), method = "jaccard")
  pcoa_res <- cmdscale(dist_matrix, k = 2, eig = TRUE)
  var_explained <- round(pcoa_res$eig / sum(pcoa_res$eig) * 100, 1)
  
  pcoa_df <- data.frame(
    Medication = rownames(pcoa_res$points),
    PCoA1 = pcoa_res$points[, 1],
    PCoA2 = pcoa_res$points[, 2]
  )
  
  # --- 7. Plot the PCoA ---
  print("Generating PCoA Scatter Plot...")
  
  p_pcoa <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, label = Medication)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
    geom_point(color = "dodgerblue4", size = 3, alpha = 0.8) +
    geom_text_repel(size = 4, fontface = "bold", max.overlaps = 30) +
    theme_bw(base_size = 14) +
    labs(
      title = "PCoA of Medications Based on Microbiome Footprint",
      subtitle = "Adjusted for polypharmacy (Distance = Quantitative Jaccard)",
      x = paste0("PCoA 1 (", var_explained[1], "%)"),
      y = paste0("PCoA 2 (", var_explained[2], "%)")
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  print(p_pcoa)
  ggsave("Medication_Footprint_PCoA.png", plot = p_pcoa, width = 10, height = 8, dpi = 300)
  
  print("SUCCESS: Pipeline Complete! Plot Saved.")
  
} else {
  print("ERROR: Not enough active medications survived the scoring to run a PCoA.")
}




# ==============================================================================
# SCRIPT: Extract Raw Regression Results
# DESCRIPTION: Re-runs the multivariate models to capture raw statistics 
# and exports them to a CSV file for review.
# ==============================================================================

print("Extracting raw multivariate regression statistics...")

raw_stats_list <- list()

for (sp in target_species) {
  
  # Build and run the model
  formula_str <- paste0("`", sp, "` ~ `", paste(selected_drugs, collapse = "` + `"), "`")
  model <- lm(as.formula(formula_str), data = analysis_df)
  mod_summary <- summary(model)$coefficients
  coef_names <- rownames(mod_summary)
  
  # Extract stats for each drug
  for (drug in selected_drugs) {
    
    matching_coef <- coef_names[str_detect(coef_names, fixed(drug))]
    
    if (length(matching_coef) > 0) {
      actual_name <- matching_coef[1] 
      
      # Grab the raw numbers
      p_val    <- mod_summary[actual_name, "Pr(>|t|)"]
      estimate <- mod_summary[actual_name, "Estimate"]
      std_err  <- mod_summary[actual_name, "Std. Error"]
      t_val    <- mod_summary[actual_name, "t value"]
      
      # Re-calculate the score just so the Prof can see how it translated
      score <- 0
      if (!is.na(p_val)) {
        if (p_val <= 0.05 && estimate > 0) score <- 2
        else if (p_val > 0.05 && p_val <= 0.1 && estimate > 0) score <- 1
        else if (p_val > 0.05 && p_val <= 0.1 && estimate < 0) score <- -1
        else if (p_val <= 0.05 && estimate < 0) score <- -2
      }
      
      # Save to dataframe row
      row_data <- data.frame(
        Target_Species = sp,
        Medication = drug,
        Beta_Estimate = estimate,
        Std_Error = std_err,
        t_value = t_val,
        P_value = p_val,
        Assigned_Score = score,
        stringsAsFactors = FALSE
      )
      
      raw_stats_list[[paste(sp, drug, sep="_")]] <- row_data
    }
  }
}

# Combine everything into one master table
master_regression_stats <- bind_rows(raw_stats_list)

# Sort it so the most significant hits (lowest p-values) are at the top
master_regression_stats <- master_regression_stats %>%
  arrange(P_value)

# Export to CSV
write.csv(master_regression_stats, "Raw_Multivariate_Regression_Results.csv", row.names = FALSE)

print("SUCCESS: 'Raw_Multivariate_Regression_Results.csv' has been saved to your working directory!")



# ==============================================================================
# SCRIPT: Generate_Fig4_Medication_Enrichment.R (UPDATED)
# DESCRIPTION: Runs multivariate regression on ALL species to identify
# significant positive associations (FDR <= 0.05). Computes clade-specific
# percentages and runs Fisher's Exact Test against "All Others" background.
# Plots the stacked Figure 4 heatmaps.
# ==============================================================================

library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork) # For stacking the heatmaps cleanly

print("==================================================")
print("STARTING FIGURE 4 ENRICHMENT PIPELINE")
print("==================================================")

# --- 1. Prepare Data for ALL Species ---
print("Preparing master dataframe for all species...")

all_abundance_df <- AllCombinedSpProfile

# THE FIX: Force all column names to be unique to resolve duplicates 
# caused by previous taxonomic renaming steps.
colnames(all_abundance_df) <- make.unique(colnames(all_abundance_df))

# Grab the newly unique column names
all_species <- colnames(all_abundance_df)
all_abundance_df$Sample <- rownames(all_abundance_df)

# Merge ALL species with the medication matrix
analysis_df_all <- all_abundance_df %>%
  inner_join(drug_matrix[, c("Sample", selected_drugs)], by = "Sample")


# --- 2. Run Multivariate Regressions on All Species ---
print(paste("Running multivariate models for all", length(all_species), "species. This may take ~30 seconds..."))

results_list <- list()

for (sp in all_species) {
  # Formula: Species ~ Drug1 + Drug2 + ... + Drug42
  formula_str <- paste0("`", sp, "` ~ `", paste(selected_drugs, collapse = "` + `"), "`")
  
  model <- lm(as.formula(formula_str), data = analysis_df_all)
  mod_summary <- summary(model)$coefficients
  coef_names <- rownames(mod_summary)
  
  for (drug in selected_drugs) {
    # Safely find the coefficient name, ignoring backtick quirks
    matching_coef <- coef_names[str_detect(coef_names, fixed(drug))]
    
    if (length(matching_coef) > 0) {
      actual_name <- matching_coef[1]
      
      results_list[[paste(sp, drug, sep="_")]] <- data.frame(
        Species = sp,
        Drug = drug,
        Estimate = mod_summary[actual_name, "Estimate"],
        Pvalue = mod_summary[actual_name, "Pr(>|t|)"],
        stringsAsFactors = FALSE
      )
    }
  }
}

all_results <- bind_rows(results_list)
print("All models complete!")


# --- 3. Apply FDR Correction Globally ---
print("Applying Benjamini-Hochberg FDR correction...")
all_results$FDR <- p.adjust(all_results$Pvalue, method = "BH")

# Define "Significant Positive Association" (Estimate > 0 AND FDR <= 0.05)
all_results$SigPos <- ifelse(all_results$Estimate > 0 & all_results$FDR <= 0.05, 1, 0)


# --- 4. Categorize Species into Genera ---
all_results <- all_results %>%
  mutate(Group = case_when(
    str_detect(Species, "^Streptococcus_") ~ "Streptococcus",
    str_detect(Species, "^Veillonella_") ~ "Veillonella",
    str_detect(Species, "^Rothia_") ~ "Rothia",
    TRUE ~ "All Others"
  ))

# Order the factors to match the manuscript top-to-bottom
group_order <- c("Streptococcus", "Veillonella", "Rothia", "All Others")
all_results$Group <- factor(all_results$Group, levels = group_order)


# --- 5. Calculate Extent of Association (TOP HEATMAP DATA) ---
print("Calculating Clade percentages...")

percent_df <- all_results %>%
  group_by(Drug, Group) %>%
  summarise(
    Total_Species = n(),
    Sig_Count = sum(SigPos),
    Percentage = Sig_Count / Total_Species,
    .groups = 'drop'
  )

# Sort drugs based on total impact (so highly impactful drugs group together)
drug_order <- percent_df %>%
  group_by(Drug) %>%
  summarise(Total_Sig = sum(Sig_Count)) %>%
  arrange(desc(Total_Sig)) %>%
  pull(Drug)

percent_df$Drug <- factor(percent_df$Drug, levels = drug_order)


# --- 6. Fisher's Exact Test for Enrichment (BOTTOM HEATMAP DATA) ---
print("Running Fisher's Exact Tests for Enrichment vs. Background...")

fisher_results <- list()
focal_groups <- c("Streptococcus", "Veillonella", "Rothia")

for (d in selected_drugs) {
  drug_data <- all_results %>% filter(Drug == d)
  
  # Background ("All Others") frequencies
  others_sig <- sum(drug_data$Group == "All Others" & drug_data$SigPos == 1)
  others_nonsig <- sum(drug_data$Group == "All Others" & drug_data$SigPos == 0)
  
  for (g in focal_groups) {
    g_sig <- sum(drug_data$Group == g & drug_data$SigPos == 1)
    g_nonsig <- sum(drug_data$Group == g & drug_data$SigPos == 0)
    
    # 2x2 Contingency Table for Fisher's Test
    contingency_table <- matrix(c(g_sig, others_sig, g_nonsig, others_nonsig), nrow = 2)
    
    # We want to know if the focal group has a GREATER proportion of hits than background
    p_val <- fisher.test(contingency_table, alternative = "greater")$p.value
    
    fisher_results[[paste(d, g, sep="_")]] <- data.frame(
      Drug = d,
      Group = g,
      Fisher_P = p_val,
      stringsAsFactors = FALSE
    )
  }
}

fisher_df <- bind_rows(fisher_results)
fisher_df$Drug <- factor(fisher_df$Drug, levels = drug_order)
fisher_df$Group <- factor(fisher_df$Group, levels = focal_groups)

# Apply FDR to the Fisher P-values
fisher_df$Fisher_FDR <- p.adjust(fisher_df$Fisher_P, method = "BH")

# Create categories matching the manuscript's legend
fisher_df <- fisher_df %>%
  mutate(Significance = case_when(
    Fisher_FDR <= 0.05 ~ "FDR <= 0.05",
    Fisher_FDR > 0.05 & Fisher_FDR <= 0.1 ~ "0.1 >= FDR > 0.05",
    TRUE ~ "Not Significant"
  )) %>%
  mutate(Label = case_when(
    Fisher_FDR <= 0.05 ~ "*",
    Fisher_FDR > 0.05 & Fisher_FDR <= 0.1 ~ "@",
    TRUE ~ ""
  ))

fisher_df$Significance <- factor(fisher_df$Significance, 
                                 levels = c("FDR <= 0.05", "0.1 >= FDR > 0.05", "Not Significant"))


# --- 7. Generate Heatmaps ---
print("Generating Figure 4 Plots...")

# TOP PLOT: Extent of Association
p_top <- ggplot(percent_df, aes(x = Drug, y = Group, fill = Percentage)) +
  geom_tile(color = "gray80", linewidth = 0.2) +
  scale_fill_gradient(low = "white", high = "#cc4c02", limits = c(0, 1), name = "Percentage\n(FDR <= 0.05)") +
  scale_y_discrete(limits = rev(group_order)) + # Reverse so Streptococcus is on top
  theme_minimal(base_size = 12) +
  labs(x = NULL, y = "Microbial Groups", 
       title = "Percentage of detected microbes within each clade showing significant positive association") +
  theme(
    axis.text.x = element_blank(), # Hide x text because it will be shared with the bottom plot
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", color = "black"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )

# BOTTOM PLOT: Fisher Enrichment Significance
p_bottom <- ggplot(fisher_df, aes(x = Drug, y = Group, fill = Significance)) +
  geom_tile(color = "gray80", linewidth = 0.2) +
  geom_text(aes(label = Label), color = "white", size = 6, vjust = 0.75) +
  scale_fill_manual(values = c("FDR <= 0.05" = "#08519c", 
                               "0.1 >= FDR > 0.05" = "#4292c6", 
                               "Not Significant" = "gray95"),
                    name = "Enrichment\nSignificance") +
  scale_y_discrete(limits = rev(focal_groups)) +
  theme_minimal(base_size = 12) +
  labs(x = "Medications", y = NULL, 
       title = "Significance of positive association enrichment (vs. All Other Lineages)") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )

# Combine the plots using patchwork
final_fig4 <- p_top / p_bottom + plot_layout(heights = c(1.5, 1))

print(final_fig4)
ggsave("Figure4_Medication_Enrichment.png", plot = final_fig4, width = 16, height = 10, dpi = 300)

print("==================================================")
print("SUCCESS: Figure 4 Generated and Saved!")