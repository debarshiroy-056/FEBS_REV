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

# ==============================================================
# SEX-SPECIFIC AGE ASSOCIATION PIPELINE
# ==============================================================

# --------------------------------------------------------------
# 1. Load libraries
# --------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)

# --------------------------------------------------------------
# 2. Clean metadata
# --------------------------------------------------------------
AllCombinedMetadata$age_num <- as.numeric(AllCombinedMetadata$age)

AllCombinedMetadata$gender <- tolower(AllCombinedMetadata$gender)

AllCombinedMetadata$gender[
  AllCombinedMetadata$gender %in% c("m","male")
] <- "male"

AllCombinedMetadata$gender[
  AllCombinedMetadata$gender %in% c("f","female")
] <- "female"

# --------------------------------------------------------------
# 3. Select metadata columns
# --------------------------------------------------------------
metadata_selected <- AllCombinedMetadata[,c(
  "sample_id",
  "age_num",
  "gender",
  "study_name",
  "diseaseCat"
)]

# --------------------------------------------------------------
# 4. Select the 201 taxa used in pipeline
# --------------------------------------------------------------
all_201_taxa <- rownames(SpeciesScores_20241021)

common_taxa <- intersect(
  all_201_taxa,
  colnames(AllCombinedSpProfile)
)

species_profile_filtered <- AllCombinedSpProfile[,common_taxa]

# --------------------------------------------------------------
# 5. Merge metadata + species
# --------------------------------------------------------------
master_data <- merge(
  metadata_selected,
  species_profile_filtered,
  by.x="sample_id",
  by.y="row.names"
)

# --------------------------------------------------------------
# 6. Age filter (41–80)
# --------------------------------------------------------------
master_data <- subset(
  master_data,
  age_num > 40 & age_num <= 80
)

# --------------------------------------------------------------
# 7. Control samples only
# --------------------------------------------------------------
controls <- subset(master_data, diseaseCat=="control")

# --------------------------------------------------------------
# 8. Split Male vs Female
# --------------------------------------------------------------
controls_male <- subset(controls, gender=="male")
controls_female <- subset(controls, gender=="female")

print(paste("Male controls:",nrow(controls_male)))
print(paste("Female controls:",nrow(controls_female)))

# ==============================================================
# FUNCTION: CREATE AGE DIRECTION MATRIX
# ==============================================================

calculate_age_direction <- function(data){
  
  study_names <- unique(data$study_name)
  
  results_list <- list()
  
  for(study in study_names){
    
    study_df <- data[data$study_name==study,]
    
    if(nrow(study_df) < 10) next
    
    age_vec <- study_df$age_num
    
    score_vec <- setNames(rep(0,length(common_taxa)),common_taxa)
    
    for(taxon in common_taxa){
      
      abundance <- study_df[[taxon]]
      
      rho <- NA
      pval <- NA
      
      tryCatch({
        
        test <- cor.test(
          abundance,
          age_vec,
          method="spearman"
        )
        
        rho <- test$estimate
        pval <- test$p.value
        
      }, error=function(e){ })
      
      if(!is.na(pval)){
        
        if(pval <= 0.1){
          
          score_vec[taxon] <- sign(rho)
          
        }
        
      }
      
    }
    
    results_list[[study]] <- score_vec
    
  }
  
  direction_matrix <- as.data.frame(do.call(cbind,results_list))
  
  direction_matrix[is.na(direction_matrix)] <- 0
  
  return(direction_matrix)
  
}

# --------------------------------------------------------------
# 9. Create direction matrices
# --------------------------------------------------------------
male_direction_matrix <- calculate_age_direction(controls_male)
female_direction_matrix <- calculate_age_direction(controls_female)

# ==============================================================
# FUNCTION: ASSOCIATION SCORE
# ==============================================================

calculate_AS <- function(direction_matrix){
  
  SP <- apply(direction_matrix,1,function(x) sum(x==1))
  SN <- apply(direction_matrix,1,function(x) sum(x==-1))
  
  T_val <- ncol(direction_matrix)
  
  pmax_val <- pmax(SP,SN) + 0.00001
  pmin_val <- pmin(SP,SN) + 0.00001
  
  consensus <- 1 - (pmin_val/pmax_val)
  consensus[pmax_val==0.00001] <- 1
  
  AS <- ((SP-SN)/T_val) * consensus
  
  result <- data.frame(
    SP=SP,
    SN=SN,
    T=T_val,
    AS=AS
  )
  
  return(result)
  
}

# --------------------------------------------------------------
# 10. Calculate Male and Female Age Association Scores
# --------------------------------------------------------------
male_AS_df <- calculate_AS(male_direction_matrix)
female_AS_df <- calculate_AS(female_direction_matrix)

male_AS <- male_AS_df$AS
female_AS <- female_AS_df$AS

names(male_AS) <- rownames(male_AS_df)
names(female_AS) <- rownames(female_AS_df)

# --------------------------------------------------------------
# 11. Combine scores
# --------------------------------------------------------------
df_sex_age <- data.frame(
  taxa = names(male_AS),
  male_AS = male_AS,
  female_AS = female_AS
)

# ==============================================================
# 12. Correlation
# ==============================================================

cor_result <- cor.test(
  df_sex_age$male_AS,
  df_sex_age$female_AS,
  method="spearman"
)

print(cor_result)

# ==============================================================
# 13. Quadrant classification
# ==============================================================

df_sex_age$group <- "neutral"

df_sex_age$group[
  df_sex_age$male_AS > 0 & df_sex_age$female_AS > 0
] <- "both_increase"

df_sex_age$group[
  df_sex_age$male_AS < 0 & df_sex_age$female_AS < 0
] <- "both_decrease"

df_sex_age$group[
  df_sex_age$male_AS > 0 & df_sex_age$female_AS < 0
] <- "male_specific"

df_sex_age$group[
  df_sex_age$male_AS < 0 & df_sex_age$female_AS > 0
] <- "female_specific"

# --------------------------------------------------------------
# 14. Count microbes per quadrant
# --------------------------------------------------------------

quadrant_counts <- table(df_sex_age$group)

print("Microbe counts per quadrant:")
print(quadrant_counts)

# ==============================================================
# 15. Publication-quality plot
# ==============================================================

ggplot(df_sex_age,
       aes(x=male_AS,y=female_AS,color=group)) +
  
  geom_point(size=4,alpha=0.85) +
  
  geom_vline(xintercept=0,linewidth=1.5) +
  geom_hline(yintercept=0,linewidth=1.5) +
  
  scale_color_manual(values=c(
    both_increase="blue",
    both_decrease="red",
    male_specific="orange",
    female_specific="green",
    neutral="grey70"
  )) +
  
  theme_bw() +
  
  labs(
    x="Male Age Association Score",
    y="Female Age Association Score",
    title="Gender-specific Age Association of Gut Microbiome"
  ) +
  
  theme(
    axis.text=element_text(size=14),
    axis.title=element_text(size=16),
    plot.title=element_text(size=18,hjust=0.5)
  )


# ==============================================================
# LIFESTYLE AGE ASSOCIATION PIPELINE
# ==============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# --------------------------------------------------------------
# 1. Clean age
# --------------------------------------------------------------
AllCombinedMetadata$age_num <- as.numeric(AllCombinedMetadata$age)

# --------------------------------------------------------------
# 2. Create lifestyle groups
# --------------------------------------------------------------
AllCombinedMetadata$lifestyle_group <- NA

AllCombinedMetadata$lifestyle_group[
  AllCombinedMetadata$cohort_life_style == "IndustrializedUrban"
] <- "Industrialized"

AllCombinedMetadata$lifestyle_group[
  AllCombinedMetadata$cohort_life_style %in%
    c("RuralTribal","UrbanRuralMixed")
] <- "NonIndustrialized"

table(AllCombinedMetadata$lifestyle_group)

# --------------------------------------------------------------
# 3. Select metadata
# --------------------------------------------------------------
metadata_selected <- AllCombinedMetadata[,c(
  "sample_id",
  "age_num",
  "study_name",
  "diseaseCat",
  "lifestyle_group"
)]

# --------------------------------------------------------------
# 4. Select 201 taxa
# --------------------------------------------------------------
all_201_taxa <- rownames(SpeciesScores_20241021)

common_taxa <- intersect(
  all_201_taxa,
  colnames(AllCombinedSpProfile)
)

species_profile_filtered <- AllCombinedSpProfile[,common_taxa]

# --------------------------------------------------------------
# 5. Merge metadata + species
# --------------------------------------------------------------
master_data <- merge(
  metadata_selected,
  species_profile_filtered,
  by.x="sample_id",
  by.y="row.names"
)

# --------------------------------------------------------------
# 6. Age filter
# --------------------------------------------------------------
master_data <- subset(
  master_data,
  age_num > 40 & age_num <= 80
)

# --------------------------------------------------------------
# 7. Controls only
# --------------------------------------------------------------
controls <- subset(master_data, diseaseCat=="control")

# --------------------------------------------------------------
# 8. Split lifestyle groups
# --------------------------------------------------------------
controls_industrialized <- subset(
  controls,
  lifestyle_group=="Industrialized"
)

controls_nonindustrialized <- subset(
  controls,
  lifestyle_group=="NonIndustrialized"
)

print(paste("Industrialized:",nrow(controls_industrialized)))
print(paste("NonIndustrialized:",nrow(controls_nonindustrialized)))

# ==============================================================
# FUNCTION: AGE DIRECTION MATRIX
# ==============================================================

calculate_age_direction <- function(data){
  
  study_names <- unique(data$study_name)
  
  results_list <- list()
  
  for(study in study_names){
    
    study_df <- data[data$study_name==study,]
    
    if(nrow(study_df) < 10) next
    
    age_vec <- study_df$age_num
    
    score_vec <- setNames(rep(0,length(common_taxa)),common_taxa)
    
    for(taxon in common_taxa){
      
      abundance <- study_df[[taxon]]
      
      rho <- NA
      pval <- NA
      
      tryCatch({
        
        test <- cor.test(
          abundance,
          age_vec,
          method="spearman"
        )
        
        rho <- test$estimate
        pval <- test$p.value
        
      }, error=function(e){ })
      
      if(!is.na(pval)){
        
        if(pval <= 0.1){
          
          score_vec[taxon] <- sign(rho)
          
        }
        
      }
      
    }
    
    results_list[[study]] <- score_vec
    
  }
  
  direction_matrix <- as.data.frame(do.call(cbind,results_list))
  
  direction_matrix[is.na(direction_matrix)] <- 0
  
  return(direction_matrix)
  
}

# --------------------------------------------------------------
# 9. Create direction matrices
# --------------------------------------------------------------
ind_direction_matrix <- calculate_age_direction(controls_industrialized)
nonind_direction_matrix <- calculate_age_direction(controls_nonindustrialized)

# ==============================================================
# FUNCTION: ASSOCIATION SCORE
# ==============================================================

calculate_AS <- function(direction_matrix){
  
  SP <- apply(direction_matrix,1,function(x) sum(x==1))
  SN <- apply(direction_matrix,1,function(x) sum(x==-1))
  
  T_val <- ncol(direction_matrix)
  
  pmax_val <- pmax(SP,SN) + 0.00001
  pmin_val <- pmin(SP,SN) + 0.00001
  
  consensus <- 1 - (pmin_val/pmax_val)
  consensus[pmax_val==0.00001] <- 1
  
  AS <- ((SP-SN)/T_val) * consensus
  
  result <- data.frame(
    SP=SP,
    SN=SN,
    T=T_val,
    AS=AS
  )
  
  return(result)
  
}

# --------------------------------------------------------------
# 10. Calculate association scores
# --------------------------------------------------------------
ind_AS_df <- calculate_AS(ind_direction_matrix)
nonind_AS_df <- calculate_AS(nonind_direction_matrix)

ind_AS <- ind_AS_df$AS
nonind_AS <- nonind_AS_df$AS

names(ind_AS) <- rownames(ind_AS_df)
names(nonind_AS) <- rownames(nonind_AS_df)

# --------------------------------------------------------------
# 11. Combine results
# --------------------------------------------------------------
df_lifestyle <- data.frame(
  taxa=names(ind_AS),
  industrialized_AS=ind_AS,
  nonindustrialized_AS=nonind_AS
)

# --------------------------------------------------------------
# 12. Correlation
# --------------------------------------------------------------
cor_result <- cor.test(
  df_lifestyle$industrialized_AS,
  df_lifestyle$nonindustrialized_AS,
  method="spearman"
)

print(cor_result)

# --------------------------------------------------------------
# 13. Quadrant classification
# --------------------------------------------------------------
df_lifestyle$group <- "neutral"

df_lifestyle$group[
  df_lifestyle$industrialized_AS > 0 &
    df_lifestyle$nonindustrialized_AS > 0
] <- "both_increase"

df_lifestyle$group[
  df_lifestyle$industrialized_AS < 0 &
    df_lifestyle$nonindustrialized_AS < 0
] <- "both_decrease"

df_lifestyle$group[
  df_lifestyle$industrialized_AS > 0 &
    df_lifestyle$nonindustrialized_AS < 0
] <- "industrialized_specific"

df_lifestyle$group[
  df_lifestyle$industrialized_AS < 0 &
    df_lifestyle$nonindustrialized_AS > 0
] <- "nonindustrialized_specific"

# --------------------------------------------------------------
# 14. Count microbes
# --------------------------------------------------------------
print(table(df_lifestyle$group))

# --------------------------------------------------------------
# 15. Plot
# --------------------------------------------------------------
ggplot(df_lifestyle,
       aes(x=industrialized_AS,
           y=nonindustrialized_AS,
           color=group)) +
  
  geom_point(size=4) +
  
  geom_vline(xintercept=0,linewidth=1.5) +
  geom_hline(yintercept=0,linewidth=1.5) +
  
  scale_color_manual(values=c(
    both_increase="blue",
    both_decrease="red",
    industrialized_specific="orange",
    nonindustrialized_specific="green",
    neutral="grey70"
  )) +
  
  theme_bw() +
  
  labs(
    x="Industrialized Age Association Score",
    y="Non-Industrialized Age Association Score",
    title="Lifestyle-specific Age Association of Gut Microbiome"
  )


# ==============================================================
# SEQUENCING TYPE AGE ASSOCIATION PIPELINE
# ==============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# --------------------------------------------------------------
# 1. Clean age
# --------------------------------------------------------------
AllCombinedMetadata$age_num <- as.numeric(AllCombinedMetadata$age)

# --------------------------------------------------------------
# 2. Select metadata
# --------------------------------------------------------------
metadata_selected <- AllCombinedMetadata[,c(
  "sample_id",
  "age_num",
  "study_name",
  "diseaseCat",
  "seq_type"
)]

# --------------------------------------------------------------
# 3. Select 201 taxa used in pipeline
# --------------------------------------------------------------
all_201_taxa <- rownames(SpeciesScores_20241021)

common_taxa <- intersect(
  all_201_taxa,
  colnames(AllCombinedSpProfile)
)

species_profile_filtered <- AllCombinedSpProfile[,common_taxa]

# --------------------------------------------------------------
# 4. Merge metadata + species
# --------------------------------------------------------------
master_data <- merge(
  metadata_selected,
  species_profile_filtered,
  by.x="sample_id",
  by.y="row.names"
)

# --------------------------------------------------------------
# 5. Age filter (41–80)
# --------------------------------------------------------------
master_data <- subset(
  master_data,
  age_num > 40 & age_num <= 80
)

# --------------------------------------------------------------
# 6. Control samples only
# --------------------------------------------------------------
controls <- subset(master_data, diseaseCat=="control")

# --------------------------------------------------------------
# 7. Split sequencing types
# --------------------------------------------------------------
controls_16s <- subset(controls, seq_type=="16s")
controls_wgs <- subset(controls, seq_type=="WGS")

print(paste("16S samples:",nrow(controls_16s)))
print(paste("WGS samples:",nrow(controls_wgs)))

# ==============================================================
# FUNCTION: CREATE AGE DIRECTION MATRIX
# ==============================================================

calculate_age_direction <- function(data){
  
  study_names <- unique(data$study_name)
  
  results_list <- list()
  
  for(study in study_names){
    
    study_df <- data[data$study_name==study,]
    
    if(nrow(study_df) < 10) next
    
    age_vec <- study_df$age_num
    
    score_vec <- setNames(rep(0,length(common_taxa)),common_taxa)
    
    for(taxon in common_taxa){
      
      abundance <- study_df[[taxon]]
      
      rho <- NA
      pval <- NA
      
      tryCatch({
        
        test <- cor.test(
          abundance,
          age_vec,
          method="spearman"
        )
        
        rho <- test$estimate
        pval <- test$p.value
        
      }, error=function(e){ })
      
      if(!is.na(pval)){
        
        if(pval <= 0.1){
          
          score_vec[taxon] <- sign(rho)
          
        }
        
      }
      
    }
    
    results_list[[study]] <- score_vec
    
  }
  
  direction_matrix <- as.data.frame(do.call(cbind,results_list))
  
  direction_matrix[is.na(direction_matrix)] <- 0
  
  return(direction_matrix)
  
}

# --------------------------------------------------------------
# 8. Create direction matrices
# --------------------------------------------------------------
dir_16s <- calculate_age_direction(controls_16s)
dir_wgs <- calculate_age_direction(controls_wgs)

# ==============================================================
# FUNCTION: ASSOCIATION SCORE
# ==============================================================

calculate_AS <- function(direction_matrix){
  
  SP <- apply(direction_matrix,1,function(x) sum(x==1))
  SN <- apply(direction_matrix,1,function(x) sum(x==-1))
  
  T_val <- ncol(direction_matrix)
  
  pmax_val <- pmax(SP,SN) + 0.00001
  pmin_val <- pmin(SP,SN) + 0.00001
  
  consensus <- 1 - (pmin_val/pmax_val)
  consensus[pmax_val==0.00001] <- 1
  
  AS <- ((SP-SN)/T_val) * consensus
  
  result <- data.frame(
    SP=SP,
    SN=SN,
    T=T_val,
    AS=AS
  )
  
  return(result)
  
}

# --------------------------------------------------------------
# 9. Calculate association scores
# --------------------------------------------------------------
AS_16s_df <- calculate_AS(dir_16s)
AS_wgs_df <- calculate_AS(dir_wgs)

AS_16s <- AS_16s_df$AS
AS_wgs <- AS_wgs_df$AS

names(AS_16s) <- rownames(AS_16s_df)
names(AS_wgs) <- rownames(AS_wgs_df)

# --------------------------------------------------------------
# 10. Combine results
# --------------------------------------------------------------
df_seq <- data.frame(
  taxa=names(AS_16s),
  AS_16s=AS_16s,
  AS_wgs=AS_wgs
)

# --------------------------------------------------------------
# 11. Correlation
# --------------------------------------------------------------
cor_result <- cor.test(
  df_seq$AS_16s,
  df_seq$AS_wgs,
  method="spearman"
)

print(cor_result)

# --------------------------------------------------------------
# 12. Quadrant classification
# --------------------------------------------------------------
df_seq$group <- "neutral"

df_seq$group[
  df_seq$AS_16s > 0 & df_seq$AS_wgs > 0
] <- "both_increase"

df_seq$group[
  df_seq$AS_16s < 0 & df_seq$AS_wgs < 0
] <- "both_decrease"

df_seq$group[
  df_seq$AS_16s > 0 & df_seq$AS_wgs < 0
] <- "16s_specific"

df_seq$group[
  df_seq$AS_16s < 0 & df_seq$AS_wgs > 0
] <- "wgs_specific"

# --------------------------------------------------------------
# 13. Count microbes
# --------------------------------------------------------------
print(table(df_seq$group))

# --------------------------------------------------------------
# 14. Plot
# --------------------------------------------------------------
ggplot(df_seq,
       aes(x=AS_16s,
           y=AS_wgs,
           color=group)) +
  
  geom_point(size=4) +
  
  geom_vline(xintercept=0,linewidth=1.5) +
  geom_hline(yintercept=0,linewidth=1.5) +
  
  scale_color_manual(values=c(
    both_increase="blue",
    both_decrease="red",
    "16s_specific"="orange",
    "wgs_specific"="green",
    neutral="grey70"
  )) +
  
  theme_bw() +
  
  labs(
    x="16S Age Association Score",
    y="WGS Age Association Score",
    title="Sequencing Technology Effect on Age Association"
  )


## Updated Plotting Script

### For gender specific
ggplot(df_sex_age,
       aes(x=male_AS, y=female_AS)) +
  
  # Background (non-highlighted points)
  geom_point(
    data=subset(df_sex_age,
                !(group %in% c("both_increase","both_decrease"))),
    color="grey80",
    alpha=0.4,
    size=3
  ) +
  
  # Highlighted points (Q1 & Q3)
  geom_point(
    data=subset(df_sex_age,
                group %in% c("both_increase","both_decrease")),
    aes(color=group),
    size=4,
    alpha=0.9
  ) +
  
  # Best fit line
  geom_smooth(method="lm", se=FALSE, color="black", linewidth=1.2) +
  
  geom_vline(xintercept=0, linewidth=1.5) +
  geom_hline(yintercept=0, linewidth=1.5) +
  
  scale_color_manual(values=c(
    both_increase="blue",
    both_decrease="red"
  )) +
  
  theme_bw() +
  
  labs(
    x="Male Age Association Score",
    y="Female Age Association Score",
    title="Gender-specific Age Association of Gut Microbiome"
  )

## For cohort life style

ggplot(df_lifestyle,
       aes(x=industrialized_AS,
           y=nonindustrialized_AS)) +
  
  geom_point(
    data=subset(df_lifestyle,
                !(group %in% c("both_increase","both_decrease"))),
    color="grey80",
    alpha=0.4,
    size=3
  ) +
  
  geom_point(
    data=subset(df_lifestyle,
                group %in% c("both_increase","both_decrease")),
    aes(color=group),
    size=4,
    alpha=0.9
  ) +
  
  geom_smooth(method="lm", se=FALSE, color="black", linewidth=1.2) +
  
  geom_vline(xintercept=0, linewidth=1.5) +
  geom_hline(yintercept=0, linewidth=1.5) +
  
  scale_color_manual(values=c(
    both_increase="blue",
    both_decrease="red"
  )) +
  
  theme_bw() +
  
  labs(
    x="Industrialized Age Association Score",
    y="Non-Industrialized Age Association Score",
    title="Lifestyle-specific Age Association of Gut Microbiome"
  )

## For Sequence Type

ggplot(df_seq,
       aes(x=AS_16s,
           y=AS_wgs)) +
  
  geom_point(
    data=subset(df_seq,
                !(group %in% c("both_increase","both_decrease"))),
    color="grey80",
    alpha=0.4,
    size=3
  ) +
  
  geom_point(
    data=subset(df_seq,
                group %in% c("both_increase","both_decrease")),
    aes(color=group),
    size=4,
    alpha=0.9
  ) +
  
  geom_smooth(method="lm", se=FALSE, color="black", linewidth=1.2) +
  
  geom_vline(xintercept=0, linewidth=1.5) +
  geom_hline(yintercept=0, linewidth=1.5) +
  
  scale_color_manual(values=c(
    both_increase="blue",
    both_decrease="red"
  )) +
  
  theme_bw() +
  
  labs(
    x="16S Age Association Score",
    y="WGS Age Association Score",
    title="Sequencing Technology Effect on Age Association"
  )




# ==============================================================
# EXACT SUMMARY OF STUDIES & SAMPLES FOR ALL 6 COHORTS
# ==============================================================

# Helper function to extract the counts
get_cohort_stats <- function(df, cohort_name) {
  n_samples <- nrow(df)
  n_studies <- length(unique(df$study_name))
  return(data.frame(
    Cohort = cohort_name, 
    Number_of_Studies = n_studies, 
    Number_of_Samples = n_samples
  ))
}

# Compile the statistics into a single table
summary_table <- do.call(rbind, list(
  get_cohort_stats(controls_male, "Male"),
  get_cohort_stats(controls_female, "Female"),
  get_cohort_stats(controls_industrialized, "Industrialized"),
  get_cohort_stats(controls_nonindustrialized, "Non-Industrialized"),
  get_cohort_stats(controls_16s, "16S"),
  get_cohort_stats(controls_wgs, "WGS")
))

# Print the final result
print("--------------------------------------------------")
print("EXACT SUMMARY: COHORT STATISTICS (Age 41-80 Controls)")
print("--------------------------------------------------")
print(summary_table)
