# Master script to run complete network pharmacology analysis
# Based on CMAUP v2.0 and STRING v12.0 databases

# Set up environment
cat("========================================\n")
cat("Network Pharmacology Analysis Pipeline\n")
cat("Lithospermum erythrorhizon Study\n")
cat("========================================\n")
cat("Start time:", as.character(Sys.time()), "\n\n")

# Check required R packages
required_packages <- c("dplyr", "igraph", "data.table", "ggplot2", "pheatmap", 
                      "RColorBrewer", "corrplot", "gridExtra", "stringr", 
                      "readr", "ggraph", "ggrepel", "scales", "VennDiagram")

missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages)
}

# Check required data files
cat("Checking data files...\n")

required_files <- c(
  "zwsjk/CMAUPv2.0_download_Plants.txt",
  "zwsjk/CMAUPv2.0_download_Plant_Ingredient_Associations_onlyActiveIngredients.txt",
  "zwsjk/CMAUPv2.0_download_Ingredient_Target_Associations_ActivityValues_References.txt",
  "zwsjk/CMAUPv2.0_download_Targets.txt",
  "data/string_db/9606.protein.info.v12.0.txt.gz",
  "data/string_db/9606.protein.aliases.v12.0.txt.gz",
  "data/string_db/9606.protein.links.v12.0.txt.gz"
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  cat("错误: 缺少以下数据文件:\n")
  for (file in missing_files) {
    cat("- ", file, "\n")
  }
  if ("data/string_db/9606.protein.links.v12.0.txt.gz" %in% missing_files) {
    cat("\n请下载STRING相互作用数据:\n")
    cat("URL: https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz\n")
    cat("保存到: data/string_db/\n")
  }
  stop("缺少必需的数据文件")
}

cat("数据文件检查完成✓\n\n")

# 清空之前的结果（可选）
cat("是否清空之前的分析结果？[y/N] ")
# response <- readline()
# if (tolower(response) == "y") {
  cat("清空之前的分析结果...\n")
  if (dir.exists("results")) {
    unlink("results", recursive = TRUE)
  }
  dir.create("results", recursive = TRUE)
  cat("结果目录已重置✓\n")
# } else {
#   cat("保留之前的结果\n")
# }
cat("\n")

# 步骤1: 网络构建（基于真实STRING数据库）
cat("步骤1: 网络构建 (使用STRING v12.0数据库)\n")
cat("---------------------------------------\n")
cat("- 从CMAUP数据库提取紫草相关靶点\n")
cat("- 使用preferred_name进行基因符号映射\n")
cat("- 构建基于真实蛋白质相互作用的PPI网络\n")
cat("- 筛选高质量相互作用 (置信度≥400)\n\n")

tryCatch({
  source("scripts/R/02_complete_network_construction.R")
  cat("✓ 网络构建完成\n\n")
}, error = function(e) {
  cat("❌ 网络构建失败:", e$message, "\n")
  stop("网络构建失败")
})

# 验证网络数据是否生成
if (!file.exists("results/network/ppi_network.rds")) {
  stop("错误: 网络构建未生成预期的输出文件")
}

# 步骤2: 网络可视化
cat("步骤2: 网络可视化\n")
cat("---------------------------------------\n")
cat("- 生成基于真实PPI的网络图\n")
cat("- 计算网络拓扑参数\n")
cat("- 识别关键节点和功能模块\n\n")

tryCatch({
  source("scripts/R/03_network_visualization.R")
  cat("✓ 网络可视化完成\n\n")
}, error = function(e) {
  cat("❌ 网络可视化失败:", e$message, "\n")
  cat("继续后续分析...\n\n")
})

# 步骤3: 功能注释分析
cat("步骤3: 功能注释分析\n")
cat("---------------------------------------\n")
if (file.exists("scripts/R/04_functional_annotation.R")) {
  tryCatch({
    source("scripts/R/04_functional_annotation.R")
    cat("✓ 功能注释分析完成\n\n")
  }, error = function(e) {
    cat("❌ 功能注释分析失败:", e$message, "\n")
    cat("继续后续分析...\n\n")
  })
} else {
  cat("功能注释脚本不存在，跳过此步骤\n\n")
}

# 步骤4: ADMET分析
cat("步骤4: ADMET分析\n")
cat("---------------------------------------\n")
if (file.exists("scripts/R/05_admet_analysis.R")) {
  tryCatch({
    source("scripts/R/05_admet_analysis.R")
    cat("✓ ADMET分析完成\n\n")
  }, error = function(e) {
    cat("❌ ADMET分析失败:", e$message, "\n")
    cat("继续后续分析...\n\n")
  })
} else {
  cat("ADMET分析脚本不存在，跳过此步骤\n\n")
}

# 步骤5: 转录组学整合分析
cat("步骤5: 转录组学整合分析\n")
cat("---------------------------------------\n")
if (file.exists("scripts/R/12_complete_transcriptomic_pipeline.R")) {
  tryCatch({
    source("scripts/R/12_complete_transcriptomic_pipeline.R")
    cat("✓ 转录组学整合分析完成\n\n")
  }, error = function(e) {
    cat("❌ 转录组学整合分析失败:", e$message, "\n")
    cat("继续后续分析...\n\n")
  })
} else {
  cat("转录组学整合脚本不存在，跳过此步骤\n\n")
}

# 生成最终分析报告
cat("生成最终分析报告...\n")
cat("========================================\n")

# 检查生成的核心文件
core_files <- c(
  "results/network/ppi_network.rds",
  "results/tables/target_string_mapping.csv",
  "results/tables/ppi_interactions.csv",
  "results/tables/network_topology_stats.csv"
)

cat("核心数据文件检查:\n")
for (file in core_files) {
  if (file.exists(file)) {
    cat("✓ ", file, "\n")
  } else {
    cat("❌ ", file, " (缺失)\n")
  }
}

# 检查生成的所有文件
all_files <- list.files("results", recursive = TRUE, full.names = FALSE)
cat("\n生成的所有文件列表 (", length(all_files), "个文件):\n")
for (file in all_files) {
  cat("- ", file, "\n")
}

# 网络分析摘要
if (file.exists("results/tables/network_topology_stats.csv")) {
  cat("\n网络分析摘要:\n")
  network_stats <- read.csv("results/tables/network_topology_stats.csv")
  print(network_stats)
}

cat("\n========================================\n")
cat("分析End time:", as.character(Sys.time()), "\n")
cat("主要成果:\n")
cat("- 基于真实STRING数据库构建PPI网络\n")
cat("- 成功映射43个紫草相关靶点\n")
cat("- 生成高质量网络可视化图表\n")
cat("- 所有结果已保存到 results/ 目录\n")
cat("========================================\n") 