#!/usr/bin/env Rscript
# ==============================================================================
# 数据加载、清理与预处理脚本 (V3 - Curated)
# 针对紫草(Lithospermum erythrorhizon)网络药理学研究
#
# 职责:
# 1. 加载与紫草关联的原始数据。
# 2. **执行关键的数据审查与清理步骤，剔除已知的非植物源污染物。**
# 3. 提取清理后化合物列表的靶点信息。
# 4. 为所有下游分析 (ADMET, 网络构建等) 提供一个统一、干净的输入源。
# ==============================================================

set.seed(42)
# options(warn = -1) # 在调试时建议注释掉此行以查看所有警告

# 加载必需包
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(jsonlite)
})

cat("=== 紫草网络药理学数据加载、清理与预处理 (V3 - Curated) ===\n")
cat("开始时间:", as.character(Sys.time()), "\n\n")

# ==============================================================================
# 1. 设置数据路径并加载原始数据
# ==============================================================================

cat("### 步骤 1: 加载CMAUP v2.0原始数据文件 ###\n")

# 定义文件路径
data_files <- list(
  plants = "zwsjk/CMAUPv2.0_download_Plants.txt",
  ingredients = "zwsjk/CMAUPv2.0_download_Ingredients_onlyActive.txt",
  plant_ingredients = "zwsjk/CMAUPv2.0_download_Plant_Ingredient_Associations_onlyActiveIngredients.txt",
  targets = "zwsjk/CMAUPv2.0_download_Targets.txt",
  ingredient_targets = "zwsjk/CMAUPv2.0_download_Ingredient_Target_Associations_ActivityValues_References.txt"
)

# 检查文件是否存在
missing_files <- Filter(Negate(file.exists), data_files)
if (length(missing_files) > 0) {
  stop("错误：缺失数据文件: ", paste(names(missing_files), collapse = ", "), "\n请确保在项目根目录运行。")
}
cat("✓ 所有必需的数据文件均存在。\n")

# 加载数据
plants_raw <- read_tsv(data_files$plants, show_col_types = FALSE)
plant_ingredients_raw <- read_tsv(data_files$plant_ingredients, col_names = c("Plant_ID", "Ingredient_ID"), show_col_types = FALSE)
ingredients_raw <- read_tsv(data_files$ingredients, show_col_types = FALSE)
ingredient_targets_raw <- read_tsv(data_files$ingredient_targets, show_col_types = FALSE)
targets_raw <- read_tsv(data_files$targets, show_col_types = FALSE)

cat("✓ 原始数据文件加载完成。\n")

# ==============================================================================
# 2. 识别并清理数据源中的已知污染物
# 这是本分析流程最关键的质量控制步骤
# ==============================================================================

cat("\n### 步骤 2: 审查并清理紫草的初始化合物列表 ###\n")

# 查找紫草ID
lithospermum_id <- "NPO11776"
lithospermum_info <- plants_raw %>% filter(Plant_ID == lithospermum_id)
if(nrow(lithospermum_info) == 0) stop("错误: 未在Plants.txt中找到紫草ID '", lithospermum_id, "'")
cat("INFO: 目标植物: ", lithospermum_info$Species_Name, " (ID: ", lithospermum_id, ")\n")

# 获取所有与紫草关联的原始成分ID
initial_ingredient_ids <- plant_ingredients_raw %>%
  filter(Plant_ID == lithospermum_id) %>%
  pull(Ingredient_ID)
cat("INFO: CMAUP数据库将", length(initial_ingredient_ids), "个化合物与紫草关联。\n")

# 定义一个基于审查发现的污染物黑名单
contaminant_ids <- c("NPC159589", "NPC234053") 
contaminant_names <- c("Chloramphenicol", "Griseofulvin")

# **执行清理**
curated_ingredient_ids <- setdiff(initial_ingredient_ids, contaminant_ids)
cat("INFO: 基于科学审查，已移除", length(contaminant_ids), "个已知的非植物源污染物:", paste(contaminant_names, collapse=", "), "\n")

# 获取最终用于所有分析的、干净的化合物信息列表
curated_ingredients <- ingredients_raw %>%
  filter(np_id %in% curated_ingredient_ids)
cat("✓ 清理完成。最终将对", nrow(curated_ingredients), "个真实的紫草化合物进行所有后续分析。\n")

# ==============================================================================
# 3. 提取清理后化合物列表的靶点信息
# ==============================================================================

cat("\n### 步骤 3: 提取", nrow(curated_ingredients), "个真实化合物的靶点信息 ###\n")

# 从清理后的成分ID列表中筛选靶点
curated_targets <- ingredient_targets_raw %>%
  filter(Ingredient_ID %in% curated_ingredients$np_id)

# 合并靶点详细信息
targets_with_details <- curated_targets %>%
  left_join(targets_raw, by = "Target_ID") %>%
  # 过滤掉没有基因名或基因名为空的记录
  filter(!is.na(Gene_Symbol) & Gene_Symbol != "") %>%
  # 去重，保留每个化合物-靶点的唯一记录
  distinct(Ingredient_ID, Gene_Symbol, .keep_all = TRUE)

cat("✓ 发现", nrow(targets_with_details), "个有效的化合物-靶点相互作用。\n")
cat("✓ 对应", length(unique(targets_with_details$Gene_Symbol)), "个独特的蛋白质靶点。\n")

# ==============================================================================
# 4. 创建输出目录并保存清理后的数据，作为下游分析的“单一数据源”
# ==============================================================================

cat("\n### 步骤 4: 保存清理后的数据以供下游脚本使用 ###\n")

# 创建输出目录
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# **重要**: 保存的文件现在是经过审查和清理的，文件名反映了这一点。
# 这个文件是下游所有分析的起点。
write_tsv(curated_ingredients, "data/processed/lithospermum_ingredients_curated.tsv")
cat("✓ 清理后的化合物列表已保存至: data/processed/lithospermum_ingredients_curated.tsv\n")

# 保存对应的靶点数据
write_tsv(targets_with_details, "data/processed/lithospermum_targets_curated.tsv")
cat("✓ 清理后的靶点列表已保存至: data/processed/lithospermum_targets_curated.tsv\n")

# (可选) 保存一份简化的数据质量报告
quality_report <- list(
  timestamp = Sys.time(),
  script_version = "V3 - Curated",
  plant_info = list(
    plant_id = lithospermum_info$Plant_ID,
    scientific_name = lithospermum_info$Species_Name
  ),
  curation_step = list(
    initial_compounds_from_db = length(initial_ingredient_ids),
    removed_contaminants_ids = contaminant_ids,
    final_authentic_compounds = nrow(curated_ingredients)
  ),
  downstream_stats = list(
    total_unique_targets = length(unique(targets_with_details$Gene_Symbol)),
    total_interactions = nrow(targets_with_details)
  ),
  next_step = "使用 'lithospermum_ingredients_curated.tsv' 作为输入，运行 02_ADMET_filtering.R"
)
write_json(quality_report, "data/processed/data_curation_report.json", pretty = TRUE, auto_unbox = TRUE)
cat("✓ 数据清理报告已保存。\n")

cat("\n=== 数据加载、清理与预处理完成 ===\n")
cat("下一步：请运行 02_ADMET_filtering.R 脚本，它将自动使用新生成的干净数据。\n")
cat("完成时间:", as.character(Sys.time()), "\n")