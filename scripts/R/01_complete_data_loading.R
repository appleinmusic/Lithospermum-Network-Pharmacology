#!/usr/bin/env Rscript
# Data loading and preprocessing script for Lithospermum erythrorhizon network pharmacology analysis

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

setwd("/Users/lgmoon/Desktop/zdhky")

cat("Loading and preprocessing data for network pharmacology analysis...\n")
cat("Start time:", as.character(Sys.time()), "\n\n")

# 1. Data file paths
cat("1. Setting up data file paths...\n")

# CMAUP database file paths
data_files <- list(
  plants = "zwsjk/CMAUPv2.0_download_Plants.txt",
  ingredients = "zwsjk/CMAUPv2.0_download_Ingredients_onlyActive.txt", 
  plant_ingredients = "zwsjk/CMAUPv2.0_download_Plant_Ingredient_Associations_onlyActiveIngredients.txt",
  targets = "zwsjk/CMAUPv2.0_download_Targets.txt",
  ingredient_targets = "zwsjk/CMAUPv2.0_download_Ingredient_Target_Associations_ActivityValues_References.txt",
  admet = "zwsjk/CMAUPv2.0_download_Human_Oral_Bioavailability_information_of_Ingredients_All.txt"
)

# 检查所有文件是否存在
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


cat("\n2. Loading and filtering Lithospermum data...\n")

plants <- read_tsv(data_files$plants, show_col_types = FALSE)
lithospermum <- plants %>%
  filter(str_detect(tolower(Species_Name), "lithospermum") & 
         str_detect(tolower(Species_Name), "erythrorhizon"))

if(nrow(lithospermum) == 0) {
  stop("Lithospermum not found(Lithospermum erythrorhizon)数据")
}

cat("✓ Found Lithospermum data，Plant ID:", lithospermum$Plant_ID, "\n")
cat("  学名:", lithospermum$Species_Name, "\n")
cat("  中文名:", lithospermum$Plant_Name, "\n")

# 加载植物-成分关联数据（无列名）
plant_ingredients_raw <- read_tsv(data_files$plant_ingredients, 
                                  col_names = c("Plant_ID", "Ingredient_ID"), 
                                  show_col_types = FALSE)

lithospermum_ingredient_ids <- plant_ingredients_raw %>%
  filter(Plant_ID == lithospermum$Plant_ID) %>%
  pull(Ingredient_ID)

cat("✓ 紫草相关成分数:", length(lithospermum_ingredient_ids), "\n")

# 加载成分详细信息
ingredients_raw <- read_tsv(data_files$ingredients, show_col_types = FALSE)
lithospermum_ingredients <- ingredients_raw %>%
  filter(np_id %in% lithospermum_ingredient_ids)

cat("✓ 成分详细信息记录数:", nrow(lithospermum_ingredients), "\n")

# 3. ADMET筛选

cat("\n3. 进行ADMET药物相似性筛选...\n")

# 加载ADMET数据
admet_raw <- read_tsv(data_files$admet, show_col_types = FALSE)

# 合并成分和ADMET数据 - 注意ADMET文件使用的是 Ingredient_ID，不是 np_id
ingredients_with_admet <- lithospermum_ingredients %>%
  left_join(admet_raw, by = c("np_id" = "Ingredient_ID"), suffix = c("_ingredient", "_admet"))

# 调试：检查合并后的列名
cat("合并后的列名:", paste(names(ingredients_with_admet), collapse = ", "), "\n")

# 应用严格的ADMET筛选标准 - 使用正确的列名
admet_filtered <- ingredients_with_admet %>%
  filter(
    !is.na(MW_admet), !is.na(XLOGP3), !is.na(TPSA_admet), !is.na(RTB),
    MW_admet >= 150 & MW_admet <= 800,        # 分子量范围
    XLOGP3 >= -2 & XLOGP3 <= 5,              # 脂水分配系数
    TPSA_admet <= 140,                       # 极性表面积
    RTB <= 10                                # 可旋转键数
  ) %>%
  distinct(np_id, .keep_all = TRUE) %>%
  # 重命名MW列避免冲突
  mutate(MW = MW_admet, LogP = XLOGP3, TPSA = TPSA_admet, nRot = RTB) %>%
  select(-MW_ingredient, -MW_admet, -TPSA_ingredient, -TPSA_admet)

cat("✓ ADMET筛选结果:\n")
cat("  原始成分数:", nrow(lithospermum_ingredients), "\n")
cat("  筛选后成分数:", nrow(admet_filtered), "\n")
cat("  筛选成功率:", round(nrow(admet_filtered)/nrow(lithospermum_ingredients)*100, 1), "%\n")

# 4. 加载靶点数据

cat("\n4. 加载成分-靶点相互作用数据...\n")

# 加载成分-靶点关联
ingredient_targets_raw <- read_tsv(data_files$ingredient_targets, show_col_types = FALSE)

# 筛选与ADMET筛选后成分相关的靶点
lithospermum_targets <- ingredient_targets_raw %>%
  filter(Ingredient_ID %in% admet_filtered$np_id)

cat("✓ 成分-靶点相互作用数:", nrow(lithospermum_targets), "\n")

# 加载靶点详细信息
targets_raw <- read_tsv(data_files$targets, show_col_types = FALSE)

# 合并靶点详细信息 - 注意成分-靶点文件使用 Target_ID，靶点文件可能使用不同的ID列
targets_with_details <- lithospermum_targets %>%
  left_join(targets_raw, by = c("Target_ID" = "Target_ID")) %>%
  filter(!is.na(Gene_Symbol))

cat("✓ 有效靶点记录数:", nrow(targets_with_details), "\n")
cat("✓ 独特靶点数:", length(unique(targets_with_details$Gene_Symbol)), "\n")

# 5. 数据质量验证

cat("\n5. 进行数据质量验证...\n")

# 定义验证规则
validation_rules <- validator(
  # 成分数据验证
  ingredient_id_not_na = !is.na(np_id),
  ingredient_name_not_empty = nchar(as.character(pref_name)) > 0,
  mw_positive = MW > 0,
  logp_range = LogP >= -10 & LogP <= 10,
  tpsa_positive = TPSA >= 0,
  
  # 靶点数据验证
  gene_symbol_not_empty = nchar(as.character(Gene_Symbol)) > 0,
  protein_name_not_empty = nchar(as.character(Protein_Name)) > 0
)

# 验证成分数据
ingredient_validation <- confront(admet_filtered, validation_rules[1:5])
target_validation <- confront(targets_with_details, validation_rules[6:7])

cat("✓ 成分数据验证通过率:", round(mean(values(ingredient_validation), na.rm = TRUE)*100, 1), "%\n")
cat("✓ 靶点数据验证通过率:", round(mean(values(target_validation), na.rm = TRUE)*100, 1), "%\n")

# 6. 创建输出目录并保存数据

cat("\n6. 保存预处理数据...\n")

# 创建输出目录
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# 保存筛选后的成分数据
write_tsv(admet_filtered, "data/processed/lithospermum_ingredients_filtered.tsv")

# 保存靶点数据
write_tsv(targets_with_details, "data/processed/lithospermum_targets.tsv")

# 保存植物信息
write_tsv(lithospermum, "data/processed/lithospermum_plant_info.tsv")

# 创建数据质量报告
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

cat("✓ 筛选后活性成分数:", nrow(admet_filtered), "\n")
cat("✓ 相关靶点数:", length(unique(targets_with_details$Gene_Symbol)), "\n")
cat("✓ 成分-靶点相互作用数:", nrow(targets_with_details), "\n")
cat("✓ 数据质量报告已保存\n")
cat("End time:", as.character(Sys.time()), "\n")