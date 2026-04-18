# ==============================================================================
# SINGLE COHORT-STACKED BINARY HEATMAP FROM FOCAL SPECIES-MEDICATION TABLE
# ==============================================================================
# This script reads the curated Excel table:
#   Focal_Species_Cohort_Medication_Table_TSG.xlsx
#
# Expected logical columns in the worksheet:
#   1. Cohort
#   2. Species
#   3. Corresponding Medication
#
# Outputs:
#   - one single stacked binary heatmap with Japanese-4D on top and MetaCardis below
#   - one PNG and one PDF figure
#   - one long-format carpet CSV
#   - one wide-format carpet CSV
# ==============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

print("==================================================")
print("STARTING SINGLE COHORT-STACKED BINARY HEATMAP")
print("==================================================")

input_file <- "Focal_Species_Cohort_Medication_Table_TSG.xlsx"
output_dir <- "Binary_Heatmap_Cohort_Outputs"

if (!file.exists(input_file)) {
  stop(paste0("Input file not found: ", input_file))
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# The workbook contains a title row above the true header row, so we skip 1 row.
heatmap_table <- read_excel(input_file, sheet = 1, skip = 1) %>%
  setNames(c("Cohort", "Species", "Medication")) %>%
  mutate(
    Cohort = trimws(as.character(Cohort)),
    Species = trimws(as.character(Species)),
    Medication = trimws(as.character(Medication))
  ) %>%
  filter(
    !is.na(Cohort), Cohort != "",
    !is.na(Species), Species != "",
    !is.na(Medication), Medication != ""
  ) %>%
  distinct(Cohort, Species, Medication)

print("Input table summary:")
print(table(heatmap_table$Cohort))

cohort_order <- c("Japanese-4D", "MetaCardis")
present_color <- "#6E0F1F"
absent_color <- "#F5EFE8"
grid_color <- "black"
title_color <- "#1D2733"

medication_order <- heatmap_table %>%
  count(Medication, name = "Species_Count") %>%
  arrange(desc(Species_Count), Medication) %>%
  pull(Medication)

species_order_by_cohort <- lapply(cohort_order, function(cohort_name) {
  heatmap_table %>%
    filter(Cohort == cohort_name) %>%
    count(Species, name = "Medication_Count") %>%
    arrange(desc(Medication_Count), Species) %>%
    pull(Species)
})
names(species_order_by_cohort) <- cohort_order

binary_heatmap_df <- bind_rows(lapply(cohort_order, function(cohort_name) {
  cohort_species_order <- species_order_by_cohort[[cohort_name]]

  heatmap_table %>%
    filter(Cohort == cohort_name) %>%
    mutate(Present = 1L) %>%
    complete(
      Species = cohort_species_order,
      Medication = medication_order,
      fill = list(Present = 0L)
    ) %>%
    mutate(Cohort = cohort_name)
})) %>%
  mutate(
    Row_Label = paste0(Cohort, " | ", Species),
    Row_Label = factor(
      Row_Label,
      levels = rev(
        paste0(
          rep(cohort_order, times = vapply(species_order_by_cohort, length, integer(1))),
          " | ",
          unlist(species_order_by_cohort)
        )
      )
    ),
    Medication = factor(Medication, levels = medication_order),
    Cohort = factor(Cohort, levels = cohort_order)
  )

# Carpet exports: long format and wide format.
binary_carpet_long <- binary_heatmap_df %>%
  mutate(
    Row_Label = as.character(Row_Label),
    Medication = as.character(Medication)
  ) %>%
  select(Cohort, Species, Row_Label, Medication, Present) %>%
  arrange(match(Cohort, cohort_order), Species, Medication)

binary_carpet_wide <- binary_heatmap_df %>%
  mutate(
    Row_Label = as.character(Row_Label),
    Medication = as.character(Medication)
  ) %>%
  select(Cohort, Species, Row_Label, Medication, Present) %>%
  pivot_wider(names_from = Medication, values_from = Present) %>%
  mutate(Cohort = factor(Cohort, levels = cohort_order),
         Row_Label = factor(Row_Label, levels = rev(levels(binary_heatmap_df$Row_Label)))) %>%
  arrange(Cohort, desc(Row_Label)) %>%
  mutate(
    Cohort = as.character(Cohort),
    Row_Label = as.character(Row_Label)
  )

write.csv(
  binary_carpet_long,
  file.path(output_dir, "Single_Heatmap_Binary_Carpet_Long.csv"),
  row.names = FALSE
)

write.csv(
  binary_carpet_wide,
  file.path(output_dir, "Single_Heatmap_Binary_Carpet_Wide.csv"),
  row.names = FALSE
)

n_species <- length(unique(binary_heatmap_df$Row_Label))
n_medications <- length(medication_order)

single_heatmap_plot <- ggplot(binary_heatmap_df, aes(x = Medication, y = Row_Label, fill = factor(Present))) +
  geom_tile(color = grid_color, linewidth = 0.45) +
  scale_fill_manual(
    values = c("0" = absent_color, "1" = present_color),
    breaks = c("1", "0"),
    labels = c("Present", "Absent"),
    name = "Binary Status"
  ) +
  labs(
    title = "Cohort-Specific Species-Medication Binary Heatmap",
    subtitle = "Japanese-4D on top, MetaCardis below",
    x = "Corresponding Medication",
    y = "Cohort and Species"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(size = 17, face = "bold", color = title_color, hjust = 0),
    plot.subtitle = element_text(size = 12, face = "bold", color = "#4A5563", hjust = 0),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 10, face = "bold", color = "black"),
    axis.text.y = element_text(size = ifelse(n_species > 25, 8.2, 10), face = "italic", color = "black"),
    axis.title.x = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "top",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    plot.margin = margin(t = 10, r = 12, b = 10, l = 10)
  )

plot_width <- max(14, 5 + 0.42 * n_medications)
plot_height <- max(9, 3 + 0.26 * n_species)

ggsave(
  file.path(output_dir, "Single_Cohort_Stacked_Binary_Heatmap.png"),
  plot = single_heatmap_plot,
  width = plot_width,
  height = plot_height,
  dpi = 400,
  bg = "white"
)

ggsave(
  file.path(output_dir, "Single_Cohort_Stacked_Binary_Heatmap.pdf"),
  plot = single_heatmap_plot,
  width = plot_width,
  height = plot_height,
  bg = "white"
)

print("==================================================")
print(paste0("SUCCESS: Single stacked heatmap and carpets saved in '", output_dir, "'"))
print("Generated files include PNG/PDF figure plus long and wide carpet CSV exports.")
