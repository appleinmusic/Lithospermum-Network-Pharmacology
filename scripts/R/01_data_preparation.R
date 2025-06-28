#!/usr/bin/env Rscript
# ==============================================================================
# 数据加载和预处理脚本 (V2 - 简化版)
# 针对紫草(Lithospermum erythrorhizon)网络药理学研究
# 职责：仅加载原始数据，不进行任何筛选，为后续ADMET分析提供干净的输入。
# ==============================================================================

set.seed(42)
options(warn = -1) # 暂时关闭警告

# 加载必需包
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(jsonlite)
})

# 使用相对路径，移除 setwd()
# setwd("/Users/lgmoon/Desktop/zdhky") 

cat("=== 紫草网络药理学数据加载和预处理 (V2) ===
")
cat("开始时间:", as.character(Sys.time()), "

")

# ==============================================================================
# 1. 设置数据路径
# ==============================================================================

cat("1. 检查并设置数据文件路径...
")

# CMAUP数据库文件路径 (使用相对路径)
data_files <- list(
  plants = "zwsjk/CMAUPv2.0_download_Plants.txt",
  ingredients = "zwsjk/CMAUPv2.0_download_Ingredients_onlyActive.txt", 
  plant_ingredients = "zwsjk/CMAUPv2.0_download_Plant_Ingredient_Associations_onlyActiveIngredients.txt",
  targets = "zwsjk/CMAUPv2.0_download_Targets.txt",
  ingredient_targets = "zwsjk/CMAUPv2.0_download_Ingredient_Target_Associations_ActivityValues_References.txt"
)

# 检查所有文件是否存在
missing_files <- c()
for(name in names(data_files)) {
  if(!file.exists(data_files[[name]])) {
    missing_files <- c(missing_files, data_files[[name]])
  }
}

if(length(missing_files) > 0) {
  stop("错误：缺失数据文件: ", paste(missing_files, collapse = ", "), "
请确保在项目根目录运行。")
}
cat("✓ 所有必需的数据文件均存在。
")

# ==============================================================================
# 2. 加载和筛选紫草相关数据
# ==============================================================================

cat("
2. 加载紫草的初始化合物和靶点数据...
")

# 加载植物数据，查找紫草
plants <- read_tsv(data_files$plants, show_col_types = FALSE)
lithospermum <- plants %>%
  filter(str_detect(tolower(Species_Name), "lithospermum") & 
         str_detect(tolower(Species_Name), "erythrorhizon"))

if(nrow(lithospermum) == 0) {
  stop("未找到紫草(Lithospermum erythrorhizon)数据")
}

cat("✓ 找到紫草数据，Plant ID:", lithospermum$Plant_ID, "
")

# 加载植物-成分关联数据
plant_ingredients_raw <- read_tsv(data_files$plant_ingredients, 
                                  col_names = c("Plant_ID", "Ingredient_ID"), 
                                  show_col_types = FALSE)

# 筛选紫草的成分ID
lithospermum_ingredient_ids <- plant_ingredients_raw %>%
  filter(Plant_ID == lithospermum$Plant_ID) %>%
  pull(Ingredient_ID)

cat("✓ 发现紫草初始化合物ID数量:", length(lithospermum_ingredient_ids), "
")

# 加载成分详细信息
ingredients_raw <- read_tsv(data_files$ingredients, show_col_types = FALSE)
lithospermum_initial_ingredients <- ingredients_raw %>%
  filter(np_id %in% lithospermum_ingredient_ids)

# **核心修改**：不再进行任何ADMET筛选，保留所有初始化合物
cat("✓ 保留所有", nrow(lithospermum_initial_ingredients), "个初始化合物用于后续分析。
")

# ==============================================================================
# 3. 加载与所有初始化合物相关的靶点数据
# ==============================================================================

cat("
3. 加载与所有初始化合物相关的靶点数据...
")

# 加载成分-靶点关联
ingredient_targets_raw <- read_tsv(data_files$ingredient_targets, show_col_types = FALSE)

# 筛选与所有紫草成分相关的靶点
all_lithospermum_targets <- ingredient_targets_raw %>%
  filter(Ingredient_ID %in% lithospermum_initial_ingredients$np_id)

cat("✓ 初始化合物-靶点相互作用总数:", nrow(all_lithospermum_targets), "
")

# 加载靶点详细信息
targets_raw <- read_tsv(data_files$targets, show_col_types = FALSE)

# 合并靶点详细信息
targets_with_details <- all_lithospermum_targets %>%
  left_join(targets_raw, by = "Target_ID") %>%
  filter(!is.na(Gene_Symbol) & Gene_Symbol != "") %>%
  distinct(Ingredient_ID, Gene_Symbol, .keep_all = TRUE) # 去重，保留每个化合物-靶点的唯一记录

cat("✓ 有效且唯一的靶点记录数:", nrow(targets_with_details), "
")
cat("✓ 独特靶点基因符号数:", length(unique(targets_with_details$Gene_Symbol)), "
")

# ==============================================================================
# 4. 创建输出目录并保存原始数据以供下一步使用
# ==============================================================================

cat("
4. 保存预处理的原始数据...
")

# 创建输出目录
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# 保存所有初始化合物数据
# **重要**：重命名输出文件以反映其内容（未经过滤）
write_tsv(lithospermum_initial_ingredients, "data/processed/lithospermum_ingredients_initial.tsv")
cat("✓ 所有初始化合物数据已保存至: data/processed/lithospermum_ingredients_initial.tsv
")

# 保存所有靶点数据
write_tsv(targets_with_details, "data/processed/lithospermum_targets_initial.tsv")
cat("✓ 所有靶点数据已保存至: data/processed/lithospermum_targets_initial.tsv
")

# 保存植物信息
write_tsv(lithospermum, "data/processed/lithospermum_plant_info.tsv")

# 创建简化的数据报告
quality_report <- list(
  timestamp = Sys.time(),
  plant_info = list(
    plant_id = lithospermum$Plant_ID,
    scientific_name = lithospermum$Species_Name
  ),
  initial_data_stats = list(
    total_initial_ingredients = nrow(lithospermum_initial_ingredients),
    total_unique_targets = length(unique(targets_with_details$Gene_Symbol)),
    total_interactions = nrow(targets_with_details)
  ),
  next_step = "Perform de novo ADMET prediction and filtering using 02_ADMET_filtering.R"
)

write_json(quality_report, "data/processed/data_quality_report.json", pretty = TRUE, auto_unbox = TRUE)
cat("✓ 数据质量报告已更新并保存。
")

cat("
=== 数据加载和预处理完成 ===
")
cat("此脚本已成功提取所有初始化合物和靶点，未进行任何ADMET筛选。
")
cat("下一步：请运行 02_ADMET_filtering.R 脚本进行标准的Lipinski规则筛选。
")
cat("完成时间:", as.character(Sys.time()), "
")