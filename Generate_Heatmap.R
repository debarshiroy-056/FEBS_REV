# ==============================================================================
# SCRIPT 2 (R): FDR, Fisher's Exact Tests, and Crystal Clear Heatmap Generation
# ==============================================================================
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(tidyr)

if (!requireNamespace("metap", quietly = TRUE)) {
  stop("Package 'metap' is required. Please install it with install.packages('metap').")
}

print("==================================================")
print("STARTING FIGURE 4 HEATMAP GENERATION")
print("==================================================")

# --- User-tunable plotting controls ---
TOP_N_DRUGS <- 42
FOCAL_GROUPS <- c("Streptococcus", "Veillonella", "Rothia")
USE_REFERENCE_DRUG_ORDER <- TRUE
REFERENCE_FIGURE_ASPECT_RATIO <- 1600 / 937
MIN_FIGURE_WIDTH_IN <- 18
EXTRA_WIDTH_FOR_LEGENDS_IN <- 6.5
PER_DRUG_WIDTH_IN <- 0.28
REFERENCE_DRUG_ORDER <- c(
  "Liraglutide", "Irbesartan", "Gliclazide", "Fluvastatin", "Sitagliptin", "Torasemide",
  "Pravastatin", "Azilsartan", "Acarbose", "Alogliptin", "Indapamide", "Repaglinide",
  "Glibenclamide", "Atenolol", "Fenofibrate", "Verapamil", "Allopurinol", "Dapagliflozin",
  "Amiodarone", "Metoprolol", "Atorvastatin", "Bisoprolol", "Apixaban", "Amlodipine",
  "Furosemide", "Enalapril", "Vildagliptin", "Carvedilol", "Nifedipine", "Warfarin",
  "Rivaroxaban", "Rosuvastatin", "Olmesartan", "Esomeprazole", "Metformin", "Prasugrel",
  "Lansoprazole", "Omeprazole", "Diltiazem", "Aspirin", "Rabeprazole", "Clopidogrel"
)

build_dendrogram_segments <- function(node, leaf_positions) {
  if (is.leaf(node)) {
    label <- attr(node, "label")
    return(list(
      x = leaf_positions[[label]],
      y = 0,
      segments = data.frame(x = numeric(0), y = numeric(0), xend = numeric(0), yend = numeric(0))
    ))
  }

  child_results <- lapply(node, build_dendrogram_segments, leaf_positions = leaf_positions)
  child_x <- vapply(child_results, `[[`, numeric(1), "x")
  child_y <- vapply(child_results, `[[`, numeric(1), "y")
  parent_y <- attr(node, "height")

  vertical_segments <- bind_rows(lapply(seq_along(child_results), function(i) {
    data.frame(x = child_x[i], y = child_y[i], xend = child_x[i], yend = parent_y)
  }))

  horizontal_segment <- data.frame(
    x = min(child_x), y = parent_y, xend = max(child_x), yend = parent_y
  )

  list(
    x = mean(range(child_x)),
    y = parent_y,
    segments = bind_rows(
      bind_rows(lapply(child_results, `[[`, "segments")),
      vertical_segments,
      horizontal_segment
    )
  )
}

# Fisher's method implementation using metap::sumlog.
# This combines the cohort-specific Fisher enrichment p-values for the same drug
# into one meta-analytic p-value while preserving the exact method used.
combine_pvalues_with_sumlog <- function(p_values) {
  valid_p <- p_values[is.finite(p_values) & !is.na(p_values)]

  if (length(valid_p) == 0) {
    return(NA_real_)
  }

  # Avoid log(0) when extremely small p-values are present.
  valid_p <- pmax(valid_p, .Machine$double.xmin)
  metap::sumlog(valid_p)$p
}

# --- 1. Load the merged master data created by the Python script ---
print("Loading Python regression results...")
all_results <- read.csv("All_Cohorts_Regression_Results.csv", stringsAsFactors = FALSE)

if (file.exists("Raw_Multivariate_Regression_Results.csv")) {
  raw_results <- read.csv("Raw_Multivariate_Regression_Results.csv", stringsAsFactors = FALSE)
  missing_from_all <- base::setdiff(unique(raw_results$Drug), unique(all_results$Drug))
  missing_from_raw <- base::setdiff(unique(all_results$Drug), unique(raw_results$Drug))
  print(paste0(
    "Data integrity check (Raw vs All): raw drugs = ", length(unique(raw_results$Drug)),
    ", all drugs = ", length(unique(all_results$Drug)),
    ", missing in All = ", length(missing_from_all),
    ", missing in Raw = ", length(missing_from_raw)
  ))
}

# --- 2. FDR Correction & Categorization ---
print("Applying FDR correction...")
all_results <- all_results %>%
  group_by(Cohort) %>%
  mutate(FDR = p.adjust(Pvalue, method = "BH")) %>%
  ungroup() %>%
  mutate(
    SigPos = ifelse(Estimate > 0 & FDR <= 0.05, 1, 0),
    Group = case_when(
      str_detect(Species, "(^|_)Streptococcus(_|$)") ~ "Streptococcus",
      str_detect(Species, "(^|_)Veillonella(_|$)") ~ "Veillonella",
      str_detect(Species, "(^|_)Rothia(_|$)") ~ "Rothia",
      TRUE ~ "All Others"
    )
  )

all_drugs <- unique(all_results$Drug)

# --- 3. Top Heatmap Data ---
print("Calculating clade percentages...")
percent_df <- all_results %>%
  group_by(Cohort, Drug, Group) %>%
  summarise(
    Total_Species = n(),
    Sig_Count = sum(SigPos),
    Percentage = Sig_Count / Total_Species,
    .groups = 'drop'
  ) %>%
  filter(Total_Species > 1 | Group == "All Others") %>%
  mutate(Display_Row = paste0(Group, " (", Cohort, ")"))

top_row_levels <- c(
  "Streptococcus (MetaCardis)", "Veillonella (MetaCardis)", "Rothia (MetaCardis)",
  "Streptococcus (Japanese-4D)", "All Others (MetaCardis)", "All Others (Japanese-4D)"
)
percent_df$Display_Row <- factor(percent_df$Display_Row, levels = rev(top_row_levels))

# --- 4. Bottom Heatmap Data (Fisher's Exact Test) ---
print("Running Fisher's Exact Tests...")
fisher_results <- list()

for (c in unique(all_results$Cohort)) {
  cohort_data <- all_results %>% filter(Cohort == c)
  for (d in all_drugs) {
    drug_data <- cohort_data %>% filter(Drug == d)

    # Compare the pooled focal lineages against all other lineages.
    focal_sig <- sum(drug_data$Group %in% FOCAL_GROUPS & drug_data$SigPos == 1)
    focal_nonsig <- sum(drug_data$Group %in% FOCAL_GROUPS & drug_data$SigPos == 0)
    others_sig <- sum(drug_data$Group == "All Others" & drug_data$SigPos == 1)
    others_nonsig <- sum(drug_data$Group == "All Others" & drug_data$SigPos == 0)

    if ((focal_sig + focal_nonsig) > 0 && (others_sig + others_nonsig) > 0) {
      contingency_table <- matrix(c(focal_sig, others_sig, focal_nonsig, others_nonsig), nrow = 2)
      p_val <- fisher.test(contingency_table, alternative = "greater")$p.value

      fisher_results[[length(fisher_results) + 1]] <- data.frame(
        Cohort = c, Drug = d, Fisher_P = p_val, stringsAsFactors = FALSE
      )
    }
  }
}
fisher_df <- bind_rows(fisher_results)

# Calculate the combined cross-cohort p-value for each drug using Fisher's method.
# We use metap::sumlog() rather than a hand-written formula so the meta-analysis
# step is explicit in the script and reproducible for review.
combined_fisher <- fisher_df %>%
  group_by(Drug) %>%
  filter(n() == 2) %>%
  summarise(
    Fisher_P = combine_pvalues_with_sumlog(Fisher_P),
    Cohort = "Combined both Cohorts (Fisher's Method)",
    .groups = "drop"
  )

fisher_df_final <- bind_rows(fisher_df, combined_fisher)
fisher_df_final$Fisher_FDR <- p.adjust(fisher_df_final$Fisher_P, method = "BH")

fisher_df_final <- fisher_df_final %>%
  mutate(
    Significance = case_when(
      Fisher_FDR <= 0.05 ~ "FDR <= 0.05",
      Fisher_FDR > 0.05 & Fisher_FDR <= 0.1 ~ "0.05 < FDR <= 0.1",
      TRUE ~ "Not Significant"
    ),
    Label = case_when(
      Fisher_FDR <= 0.05 ~ "*",
      Fisher_FDR > 0.05 & Fisher_FDR <= 0.1 ~ "@",
      TRUE ~ ""
    ),
    Display_Row = Cohort
  )

fisher_df_final$Significance <- factor(fisher_df_final$Significance, levels = c("FDR <= 0.05", "0.05 < FDR <= 0.1", "Not Significant"))

bottom_row_levels <- c(
  "Japanese-4D", "MetaCardis", "Combined both Cohorts (Fisher's Method)"
)
fisher_df_final$Display_Row <- factor(fisher_df_final$Display_Row, levels = rev(bottom_row_levels))

# Export a dedicated all-drugs Fisher-method table.
# This file keeps the original cohort-specific Fisher enrichment p-values and
# adds the Fisher's-method combined p-value from metap::sumlog() to make the
# cross-cohort aggregation fully transparent.
all_drugs_fisher_method_table <- fisher_df %>%
  select(Cohort, Drug, Fisher_P) %>%
  mutate(Cohort = paste0(Cohort, "_Fisher_P")) %>%
  tidyr::pivot_wider(names_from = Cohort, values_from = Fisher_P) %>%
  left_join(
    combined_fisher %>%
      rename(Fisher_Method_P = Fisher_P) %>%
      select(Drug, Fisher_Method_P),
    by = "Drug"
  ) %>%
  mutate(
    Fisher_Method_FDR = p.adjust(Fisher_Method_P, method = "BH"),
    Significance = case_when(
      Fisher_Method_FDR <= 0.05 ~ "FDR <= 0.05",
      Fisher_Method_FDR > 0.05 & Fisher_Method_FDR <= 0.1 ~ "0.05 < FDR <= 0.1",
      TRUE ~ "Not Significant"
    ),
    Combination_Method = "metap::sumlog"
  ) %>%
  arrange(Drug)

write.csv(
  all_drugs_fisher_method_table,
  "All_Drugs_Fisher_Method_Sumlog_Table.csv",
  row.names = FALSE
)

# --- 4b. Select the most informative drugs for a readable figure ---
signal_scores <- percent_df %>%
  filter(Group %in% FOCAL_GROUPS) %>%
  group_by(Drug) %>%
  summarise(Top_Percentage = max(Percentage, na.rm = TRUE), .groups = "drop")

fisher_scores <- fisher_df_final %>%
  filter(Display_Row == "Combined both Cohorts (Fisher's Method)") %>%
  transmute(
    Drug,
    Fisher_Score = case_when(
      Fisher_FDR <= 0.05 ~ 2,
      Fisher_FDR <= 0.1 ~ 1,
      TRUE ~ 0
    )
  )

drug_rank <- full_join(signal_scores, fisher_scores, by = "Drug") %>%
  mutate(
    Top_Percentage = ifelse(is.na(Top_Percentage), 0, Top_Percentage),
    Fisher_Score = ifelse(is.na(Fisher_Score), 0, Fisher_Score),
    Rank_Score = Top_Percentage + Fisher_Score
  ) %>%
  arrange(desc(Rank_Score), desc(Top_Percentage), Drug)

selected_drugs <- drug_rank %>%
  filter(Rank_Score > 0) %>%
  slice_head(n = TOP_N_DRUGS) %>%
  pull(Drug)

if (USE_REFERENCE_DRUG_ORDER) {
  missing_reference_drugs <- base::setdiff(REFERENCE_DRUG_ORDER, all_drugs)
  if (length(missing_reference_drugs) > 0) {
    stop(paste0(
      "Reference drug order includes drugs not found in All_Cohorts_Regression_Results.csv: ",
      paste(missing_reference_drugs, collapse = ", ")
    ))
  }
  selected_drugs <- REFERENCE_DRUG_ORDER
}

if (length(selected_drugs) == 0) {
  stop("No informative drugs were found after ranking. Please check input data.")
}

percent_df <- percent_df %>% filter(Drug %in% selected_drugs)
fisher_df_final <- fisher_df_final %>% filter(Drug %in% selected_drugs)

# Create complete grid to ensure all rows have all drugs (no missing boxes)
percent_df_complete <- as.data.frame(expand.grid(
  Display_Row = unique(percent_df$Display_Row),
  Drug = selected_drugs,
  stringsAsFactors = FALSE
)) %>%
  left_join(percent_df, by = c("Display_Row", "Drug")) %>%
  mutate(Percentage = ifelse(is.na(Percentage), 0, Percentage))
percent_df <- percent_df_complete

fisher_df_complete <- as.data.frame(expand.grid(
  Display_Row = unique(fisher_df_final$Display_Row),
  Drug = selected_drugs,
  stringsAsFactors = FALSE
)) %>%
  left_join(fisher_df_final, by = c("Display_Row", "Drug")) %>%
  mutate(
    Significance = as.character(Significance),
    Significance = ifelse(is.na(Significance), "Not Significant", Significance),
    Label = ifelse(is.na(Label), "", Label)
  ) %>%
  mutate(
    Significance = factor(Significance, levels = c("FDR <= 0.05", "0.05 < FDR <= 0.1", "Not Significant"), ordered = FALSE)
  )
fisher_df_final <- fisher_df_complete

# --- 4c. Export the heatmap carpet data for CSV/Excel delivery ---
print("Exporting heatmap carpet data...")

heatmap_export_dir <- "Heatmap_Carpet_Exports"
if (!dir.exists(heatmap_export_dir)) {
  dir.create(heatmap_export_dir, recursive = TRUE)
}

top_carpet_export <- percent_df %>%
  arrange(Display_Row, Drug) %>%
  select(Display_Row, Drug, Percentage, Cohort, Group, Total_Species, Sig_Count)

bottom_carpet_export <- fisher_df_final %>%
  arrange(Display_Row, Drug) %>%
  select(Display_Row, Drug, Cohort, Fisher_P, Fisher_FDR, Significance, Label)

top_carpet_wide <- percent_df %>%
  mutate(
    Display_Row = as.character(Display_Row),
    Drug = factor(as.character(Drug), levels = selected_drugs)
  ) %>%
  arrange(Display_Row, Drug) %>%
  select(Display_Row, Drug, Percentage) %>%
  tidyr::pivot_wider(names_from = Drug, values_from = Percentage, values_fill = 0) %>%
  select(Display_Row, dplyr::all_of(selected_drugs)) %>%
  arrange(match(Display_Row, rev(top_row_levels)))

bottom_carpet_wide <- fisher_df_final %>%
  mutate(
    Display_Row = as.character(Display_Row),
    Drug = factor(as.character(Drug), levels = selected_drugs)
  ) %>%
  arrange(Display_Row, Drug) %>%
  select(Display_Row, Drug, Significance) %>%
  tidyr::pivot_wider(names_from = Drug, values_from = Significance, values_fill = "Not Significant") %>%
  select(Display_Row, dplyr::all_of(selected_drugs)) %>%
  arrange(match(Display_Row, rev(bottom_row_levels)))

write.csv(top_carpet_export, file.path(heatmap_export_dir, "Heatmap_Carpet_Top.csv"), row.names = FALSE)
write.csv(bottom_carpet_export, file.path(heatmap_export_dir, "Heatmap_Carpet_Bottom.csv"), row.names = FALSE)
write.csv(top_carpet_wide, file.path(heatmap_export_dir, "Heatmap_Carpet_Top_Wide.csv"), row.names = FALSE)
write.csv(bottom_carpet_wide, file.path(heatmap_export_dir, "Heatmap_Carpet_Bottom_Wide.csv"), row.names = FALSE)

if (requireNamespace("openxlsx", quietly = TRUE)) {
  workbook_path <- file.path(heatmap_export_dir, "Heatmap_Carpet.xlsx")
  workbook <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(workbook, "Top_Carpet")
  openxlsx::writeData(workbook, "Top_Carpet", top_carpet_export)
  openxlsx::addWorksheet(workbook, "Bottom_Carpet")
  openxlsx::writeData(workbook, "Bottom_Carpet", bottom_carpet_export)
  openxlsx::addWorksheet(workbook, "Top_Carpet_Wide")
  openxlsx::writeData(workbook, "Top_Carpet_Wide", top_carpet_wide)
  openxlsx::addWorksheet(workbook, "Bottom_Carpet_Wide")
  openxlsx::writeData(workbook, "Bottom_Carpet_Wide", bottom_carpet_wide)
  openxlsx::addWorksheet(workbook, "Selected_Drugs")
  openxlsx::writeData(workbook, "Selected_Drugs", drug_rank %>% filter(Drug %in% selected_drugs))
  openxlsx::saveWorkbook(workbook, workbook_path, overwrite = TRUE)
}

# Cluster columns from the top-panel percentages and render a dendrogram.
clustering_matrix <- percent_df %>%
  filter(Display_Row %in% top_row_levels) %>%
  select(Drug, Display_Row, Percentage) %>%
  tidyr::pivot_wider(names_from = Display_Row, values_from = Percentage, values_fill = 0) %>%
  as.data.frame()

rownames(clustering_matrix) <- clustering_matrix$Drug
clustering_matrix$Drug <- NULL

hc <- hclust(dist(clustering_matrix))

if (USE_REFERENCE_DRUG_ORDER) {
  # Keep the publication drug order while still showing clustering structure.
  drug_order <- selected_drugs
} else {
  drug_order <- hc$labels[hc$order]
}

dendro <- as.dendrogram(hc)
leaf_positions <- as.list(stats::setNames(seq_along(drug_order), drug_order))
dendro_segments <- build_dendrogram_segments(dendro, leaf_positions)$segments

percent_df$Drug <- factor(percent_df$Drug, levels = drug_order)
fisher_df_final$Drug <- factor(fisher_df_final$Drug, levels = drug_order)

print("==================================================")
print("PART 4: GENERATING CRYSTAL CLEAR HEATMAPS")
print("==================================================")

# --- 5. Generate Crystal Clear Heatmaps ---
print("Generating Final Figure 4 Plots...")

significance_levels <- c("FDR <= 0.05", "0.05 < FDR <= 0.1", "Not Significant")
legend_dummy <- data.frame(
  Drug = rep(drug_order[1], length(significance_levels)),
  Display_Row = rep(as.character(rev(bottom_row_levels)[1]), length(significance_levels)),
  Significance = factor(significance_levels, levels = significance_levels)
)

p_dendro <- ggplot(dendro_segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.35, color = "black") +
  scale_x_continuous(limits = c(0.5, length(drug_order) + 0.5), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(
    plot.margin = margin(t = 2, r = 8, b = -18, l = 5)
  )

p_top <- ggplot(percent_df, aes(x = Drug, y = Display_Row, fill = Percentage)) +
  geom_tile(color = "gray50", linewidth = 0.4) +
  scale_fill_gradient(
    low = "white", high = "#cc4c02", limits = c(0, 1), name = "Percentage\n(FDR <= 0.05)",
    guide = guide_colorbar(
      order = 1,
      direction = "vertical",
      barwidth = grid::unit(0.6, "cm"),
      barheight = grid::unit(4, "cm"),
      title.position = "left",
      ticks = FALSE,
      label.theme = element_text(size = 10)
    )
  ) +
  
  # Removes the blank white padding around the edges of the tiles
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0), position = "right") +
  
  theme_minimal(base_size = 15) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", color = "black", size = 12),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", color = "black", size = 13),
    axis.text.y.right = element_text(face = "bold", color = "black", size = 13),
    legend.position = "left",
    legend.direction = "vertical",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.margin = margin(r = 8),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.1),
    plot.margin = margin(t = 0, r = 8, b = -4, l = 5)
  )

p_bottom <- ggplot(fisher_df_final, aes(x = Drug, y = Display_Row, fill = Significance)) +
  geom_tile(color = "gray50", linewidth = 0.4) +
  geom_text(aes(label = Label), color = "white", size = 7, vjust = 0.5, hjust = 0.5, fontface = "bold", show.legend = FALSE) + 
  geom_point(
    data = legend_dummy,
    aes(x = Drug, y = Display_Row, fill = Significance),
    inherit.aes = FALSE,
    shape = 22,
    size = 5,
    alpha = 0,
    show.legend = TRUE
  ) +
  scale_fill_manual(
    values = c("FDR <= 0.05" = "#08519c", "0.05 < FDR <= 0.1" = "#4292c6", "Not Significant" = "#d9e5f0"),
    limits = significance_levels,
    breaks = significance_levels,
    drop = FALSE,
    na.translate = FALSE,
    name = "Significance Status",
    guide = guide_legend(
      order = 1,
      ncol = 1,
      byrow = TRUE,
      override.aes = list(
        fill = c("#08519c", "#4292c6", "#d9e5f0"),
        alpha = 1,
        shape = 22,
        size = 8,
        color = "#3f6fb0"
      ),
      label.theme = element_text(size = 14, face = "bold"),
      title.theme = element_text(size = 14, face = "bold")
    )
  ) +
  
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0), position = "right") +
  
  theme_minimal(base_size = 15) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", color = "black", size = 13),
    axis.text.y.right = element_text(face = "bold", color = "black", size = 13),
    legend.position = "left",
    legend.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.key.width = grid::unit(1.5, "cm"),
    legend.key.height = grid::unit(1.5, "cm"),
    legend.background = element_rect(fill = "white", color = "#3f6fb0", linewidth = 0.8),
    legend.margin = margin(8, 8, 8, 8),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.1),
    plot.margin = margin(t = -30, r = 8, b = 8, l = 5)
  )

# Combine the plots.
final_fig4 <- p_dendro / p_top / p_bottom + plot_layout(guides = "keep", heights = c(0.45, 1.35, 0.9))

# Keep the output in the same wide landscape proportion as the reference figure.
# Width scales with the number of medications; height follows the reference aspect ratio.
plot_width <- max(MIN_FIGURE_WIDTH_IN, EXTRA_WIDTH_FOR_LEGENDS_IN + length(selected_drugs) * PER_DRUG_WIDTH_IN)
plot_height <- plot_width / REFERENCE_FIGURE_ASPECT_RATIO
ggsave("Figure4_Crystal_Clear.png", plot = final_fig4, width = plot_width, height = plot_height, dpi = 300, bg = "white")

print("==================================================")
print(paste0("SUCCESS: Crystal Clear Pipeline Complete! Check 'Figure4_Crystal_Clear.png' and '", heatmap_export_dir, "'."))
