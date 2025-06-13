#!/usr/bin/env Rscript
# ==============================================================================
# Data Loading and Preprocessing Script for Lithospermum erythrorhizon
# Network Pharmacology Analysis
# ==============================================================================

set.seed(42)
options(warn = 1)

# Load required packages
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(jsonlite)
  library(validate)
})

# Set working directory to project root (two levels up from scripts/R/)
project_root <- file.path(dirname(dirname(dirname(rstudioapi::getActiveDocumentContext()$path))))
setwd(project_root)

cat("=== Lithospermum erythrorhizon Network Pharmacology Analysis ===\n")
cat("Start time:", as.character(Sys.time()), "\n\n")

# ==============================================================================
# 1. Setup Data Paths
# ==============================================================================

cat("1. Setting up data file paths...\n")

# CMAUP database file paths (relative to project root)
data_files <- list(
  plants = "../../zwsjk/CMAUPv2.0_download_Plants.txt",
  ingredients = "../../zwsjk/CMAUPv2.0_download_Ingredients_onlyActive.txt", 
  plant_ingredients = "../../zwsjk/CMAUPv2.0_download_Plant_Ingredient_Associations_onlyActiveIngredients.txt",
  targets = "../../zwsjk/CMAUPv2.0_download_Targets.txt",
  ingredient_targets = "../../zwsjk/CMAUPv2.0_download_Ingredient_Target_Associations_ActivityValues_References.txt",
  admet = "../../zwsjk/CMAUPv2.0_download_Human_Oral_Bioavailability_information_of_Ingredients_All.txt"
)

# Check if all files exist
missing_files <- c()
for(name in names(data_files)) {
  if(!file.exists(data_files[[name]])) {
    missing_files <- c(missing_files, data_files[[name]])
  } else {
    cat("✓", data_files[[name]], "\n")
  }
}

if(length(missing_files) > 0) {
  stop("Missing data files: ", paste(missing_files, collapse = ", "))
}

# ==============================================================================
# 2. Load and Filter Lithospermum Data
# ==============================================================================

cat("\n2. Loading and filtering Lithospermum-related data...\n")

# Load plant data and search for Lithospermum
plants <- read_tsv(data_files$plants, show_col_types = FALSE)
lithospermum <- plants %>%
  filter(str_detect(tolower(Species_Name), "lithospermum") & 
         str_detect(tolower(Species_Name), "erythrorhizon"))

if(nrow(lithospermum) == 0) {
  stop("Lithospermum erythrorhizon data not found")
}

cat("✓ Found Lithospermum data, Plant ID:", lithospermum$Plant_ID, "\n")
cat("  Scientific name:", lithospermum$Species_Name, "\n")
cat("  Chinese name:", lithospermum$Plant_Name, "\n")

# Load plant-ingredient association data (no column names)
plant_ingredients_raw <- read_tsv(data_files$plant_ingredients, 
                                  col_names = c("Plant_ID", "Ingredient_ID"), 
                                  show_col_types = FALSE)

# Filter ingredients for Lithospermum
lithospermum_ingredient_ids <- plant_ingredients_raw %>%
  filter(Plant_ID == lithospermum$Plant_ID) %>%
  pull(Ingredient_ID)

cat("✓ Number of Lithospermum-related ingredients:", length(lithospermum_ingredient_ids), "\n")

# Load ingredient detailed information
ingredients_raw <- read_tsv(data_files$ingredients, show_col_types = FALSE)
lithospermum_ingredients <- ingredients_raw %>%
  filter(np_id %in% lithospermum_ingredient_ids)

cat("✓ Number of ingredient detail records:", nrow(lithospermum_ingredients), "\n")

# ==============================================================================
# 3. ADMET Screening
# ==============================================================================

cat("\n3. Performing ADMET drug-likeness screening...\n")

# Load ADMET data
admet_raw <- read_tsv(data_files$admet, show_col_types = FALSE)

# Merge ingredient and ADMET data - note ADMET file uses Ingredient_ID, not np_id
ingredients_with_admet <- lithospermum_ingredients %>%
  left_join(admet_raw, by = c("np_id" = "Ingredient_ID"), suffix = c("_ingredient", "_admet"))

# Debug: check merged column names
cat("Merged column names:", paste(names(ingredients_with_admet), collapse = ", "), "\n")

# Apply strict ADMET filtering criteria - use correct column names
admet_filtered <- ingredients_with_admet %>%
  filter(
    !is.na(MW_admet), !is.na(XLOGP3), !is.na(TPSA_admet), !is.na(RTB),
    MW_admet >= 150 & MW_admet <= 800,        # Molecular weight range
    XLOGP3 >= -2 & XLOGP3 <= 5,              # Lipophilicity
    TPSA_admet <= 140,                       # Topological polar surface area
    RTB <= 10                                # Number of rotatable bonds
  ) %>%
  distinct(np_id, .keep_all = TRUE) %>%
  # Rename MW column to avoid conflicts
  mutate(MW = MW_admet, LogP = XLOGP3, TPSA = TPSA_admet, nRot = RTB) %>%
  select(-MW_ingredient, -MW_admet, -TPSA_ingredient, -TPSA_admet)

cat("✓ ADMET screening results:\n")
cat("  Original compound count:", nrow(lithospermum_ingredients), "\n")
cat("  Filtered compound count:", nrow(admet_filtered), "\n")
cat("  Success rate:", round(nrow(admet_filtered)/nrow(lithospermum_ingredients)*100, 1), "%\n")

# ==============================================================================
# 4. Load Target Data
# ==============================================================================

cat("\n4. Loading compound-target interaction data...\n")

# Load ingredient-target associations
ingredient_targets_raw <- read_tsv(data_files$ingredient_targets, show_col_types = FALSE)

# Filter targets related to ADMET-filtered compounds
lithospermum_targets <- ingredient_targets_raw %>%
  filter(Ingredient_ID %in% admet_filtered$np_id)

cat("✓ Number of compound-target interactions:", nrow(lithospermum_targets), "\n")

# Load target detailed information
targets_raw <- read_tsv(data_files$targets, show_col_types = FALSE)

# Merge target details - note compound-target file uses Target_ID, target file may use different ID column
targets_with_details <- lithospermum_targets %>%
  left_join(targets_raw, by = c("Target_ID" = "Target_ID")) %>%
  filter(!is.na(Gene_Symbol))

cat("✓ Valid target records:", nrow(targets_with_details), "\n")
cat("✓ Unique targets:", length(unique(targets_with_details$Gene_Symbol)), "\n")

# ==============================================================================
# 5. Data Quality Validation
# ==============================================================================

cat("\n5. Performing data quality validation...\n")

# Define validation rules
validation_rules <- validator(
  # Ingredient data validation
  ingredient_id_not_na = !is.na(np_id),
  ingredient_name_not_empty = nchar(as.character(pref_name)) > 0,
  mw_positive = MW > 0,
  logp_range = LogP >= -10 & LogP <= 10,
  tpsa_positive = TPSA >= 0,
  
  # Target data validation
  gene_symbol_not_empty = nchar(as.character(Gene_Symbol)) > 0,
  protein_name_not_empty = nchar(as.character(Protein_Name)) > 0
)

# Validate ingredient data
ingredient_validation <- confront(admet_filtered, validation_rules[1:5])
target_validation <- confront(targets_with_details, validation_rules[6:7])

cat("✓ Ingredient data validation pass rate:", round(mean(values(ingredient_validation), na.rm = TRUE)*100, 1), "%\n")
cat("✓ Target data validation pass rate:", round(mean(values(target_validation), na.rm = TRUE)*100, 1), "%\n")

# ==============================================================================
# 6. Create Output Directory and Save Data
# ==============================================================================

cat("\n6. Saving preprocessed data...\n")

# Create output directory
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# Save filtered ingredient data
write_tsv(admet_filtered, "data/processed/lithospermum_ingredients_filtered.tsv")

# Save target data
write_tsv(targets_with_details, "data/processed/lithospermum_targets.tsv")

# Save plant information
write_tsv(lithospermum, "data/processed/lithospermum_plant_info.tsv")

# Create data quality report
quality_report <- list(
  timestamp = Sys.time(),
  plant_info = list(
    plant_id = lithospermum$Plant_ID,
    scientific_name = lithospermum$Species_Name,
    chinese_name = lithospermum$Plant_Name
  ),
  ingredient_stats = list(
    total_ingredients = nrow(lithospermum_ingredients),
    admet_filtered = nrow(admet_filtered),
    filter_rate = round(nrow(admet_filtered)/nrow(lithospermum_ingredients)*100, 1)
  ),
  target_stats = list(
    total_interactions = nrow(lithospermum_targets),
    valid_interactions = nrow(targets_with_details),
    unique_targets = length(unique(targets_with_details$Gene_Symbol))
  ),
  admet_criteria = list(
    MW_range = "150-800 Da",
    LogP_range = "-2 to 5",
    TPSA_max = "≤140 Ų",
    nRot_max = "≤10"
  ),
  validation_results = list(
    ingredient_pass_rate = round(mean(values(ingredient_validation), na.rm = TRUE)*100, 1),
    target_pass_rate = round(mean(values(target_validation), na.rm = TRUE)*100, 1)
  )
)

write_json(quality_report, "data/processed/data_quality_report.json", pretty = TRUE)

cat("\n=== Data Loading and Preprocessing Complete ===\n")
cat("✓ Filtered active ingredients:", nrow(admet_filtered), "\n")
cat("✓ Related targets:", length(unique(targets_with_details$Gene_Symbol)), "\n")
cat("✓ Ingredient-target interactions:", nrow(targets_with_details), "\n")
cat("✓ Data quality report saved\n")
cat("Completion time:", as.character(Sys.time()), "\n") 