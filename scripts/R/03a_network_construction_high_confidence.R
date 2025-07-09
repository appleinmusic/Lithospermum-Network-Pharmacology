# 紫草网络药理学分析 - 高置信度网络构建 (敏感性分析)
# 使用STRING数据库构建蛋白质相互作用网络 (score >= 700)

# 加载必需的包
suppressMessages({
  library(dplyr)
  library(igraph)
  library(data.table)
  library(VennDiagram)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(corrplot)
  library(gridExtra)
  library(stringr)
  library(readr)
})

# 创建输出目录
dir.create("results/network", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

cat("开始高置信度网络构建分析 (score >= 700)...\n")

# ==============================================================================
# 1. 加载通过ADMET筛选的化合物
# ==============================================================================
cat("正在加载通过ADMET筛选的化合物...\n")
admet_filtered_compounds_path <- "results/admet_analysis_data.csv"
if (!file.exists(admet_filtered_compounds_path)) {
  stop("错误: ADMET筛选结果文件 '", admet_filtered_compounds_path, "' 不存在。")
}
admet_filtered_compounds <- read_csv(admet_filtered_compounds_path, show_col_types = FALSE)
filtered_ingredient_ids <- unique(admet_filtered_compounds$np_id)
cat("✓ 加载了 ", length(filtered_ingredient_ids), " 个通过ADMET筛选的化合物。\n")


# ==============================================================================
# 2. 根据筛选后的化合物确定靶点
# ==============================================================================
cat("正在根据筛选后的化合物确定靶点...\n")
filtered_targets_path <- "data/processed/lithospermum_targets.tsv"
if (!file.exists(filtered_targets_path)) {
  stop("错误: 已处理的靶点文件 '", filtered_targets_path, "' 不存在。")
}
zicao_targets_filtered <- read_tsv(filtered_targets_path, show_col_types = FALSE)
all_targets <- unique(zicao_targets_filtered$Gene_Symbol)
all_targets <- all_targets[!is.na(all_targets) & all_targets != ""]
cat("✓ 最终用于网络构建的独特靶点数:", length(all_targets), "\n")


# 3. 加载STRING数据库文件
cat("正在加载STRING数据库...\n")
protein_info <- fread("data/string_db/9606.protein.info.v12.0.txt.gz")
protein_aliases <- fread("data/string_db/9606.protein.aliases.v12.0.txt.gz")

if (file.exists("data/string_db/9606.protein.links.v12.0.txt.gz")) {
  protein_links <- fread("data/string_db/9606.protein.links.v12.0.txt.gz")
  cat("STRING蛋白质相互作用数据加载完成，共", nrow(protein_links), "条记录\n")
  
  # ============================================================================
  # 核心修改：使用高置信度阈值 (score >= 700)
  # ============================================================================
  high_quality_links <- protein_links[protein_links$combined_score >= 700, ]
  cat("高置信度相互作用 (score>=700):", nrow(high_quality_links), "条\n")
  
} else {
  stop("缺少必需的STRING数据文件: 9606.protein.links.v12.0.txt.gz")
}

# 4. 靶点名称映射到STRING ID
cat("正在进行靶点名称映射...\n")
# (映射逻辑与原脚本一致)
map_to_string_id <- function(target_names, aliases_df, protein_info_df) {
  protein_id_col <- colnames(aliases_df)[1]
  alias_col <- "alias"
  mapped_data <- data.frame(target_name=character(0), string_id=character(0), source=character(0), stringsAsFactors=FALSE)
  for (target in target_names) {
    if (is.na(target) || target == "") next
    found <- FALSE
    protein_matches <- protein_info_df[protein_info_df$preferred_name == target, ]
    if (nrow(protein_matches) > 0) {
      protein_id <- protein_matches[[colnames(protein_info_df)[1]]][1]
      if (length(protein_id) > 0 && !is.na(protein_id) && nchar(as.character(protein_id)) > 0) {
        mapped_data <- rbind(mapped_data, data.frame(target_name=target, string_id=as.character(protein_id), source="preferred_name", stringsAsFactors=FALSE))
        found <- TRUE
      }
    }
    if (!found) {
      matches <- aliases_df[aliases_df[[alias_col]] == target, ]
      if (nrow(matches) > 0) {
        protein_id <- matches[[protein_id_col]][1]
        if (length(protein_id) > 0 && !is.na(protein_id) && nchar(as.character(protein_id)) > 0) {
          mapped_data <- rbind(mapped_data, data.frame(target_name=target, string_id=as.character(protein_id), source=paste0("alias:", matches$source[1]), stringsAsFactors=FALSE))
          found <- TRUE
        }
      }
    }
  }
  return(mapped_data)
}
target_mapping <- map_to_string_id(all_targets, protein_aliases, protein_info)
cat("成功映射靶点数:", nrow(target_mapping), "/", length(all_targets), "\n")
unmapped_targets <- setdiff(all_targets, target_mapping$target_name)
if (length(unmapped_targets) > 0) {
    cat("警告：以下", length(unmapped_targets), "个靶点未能映射到STRING ID:\n")
    for (target in unmapped_targets) { cat("  - ", target, "\n") }
}

# 5. 构建PPI网络
cat("正在构建高置信度PPI网络...\n")
if (nrow(target_mapping) > 0) {
  target_string_ids <- target_mapping$string_id
  target_interactions <- high_quality_links[
    high_quality_links$protein1 %in% target_string_ids & 
    high_quality_links$protein2 %in% target_string_ids, 
  ]
  cat("高置信度目标蛋白质间相互作用数:", nrow(target_interactions), "\n")
  
  if (nrow(target_interactions) > 0) {
    g <- graph_from_data_frame(target_interactions[, c("protein1", "protein2", "combined_score")], directed = FALSE)
    id_to_name <- setNames(target_mapping$target_name, target_mapping$string_id)
    V(g)$name_display <- ifelse(names(V(g)) %in% names(id_to_name), id_to_name[names(V(g))], names(V(g)))
    V(g)$node_type <- "紫草靶点"
    
    cat("\n正在提取最大连通分量 (LCC)...\n")
    components_list <- decompose(g)
    if (length(components_list) > 0) {
        largest_component_index <- which.max(sapply(components_list, vcount))
        g_core_network <- components_list[[largest_component_index]]
        cat("信息：从", vcount(g), "个节点的网络中提取出", vcount(g_core_network), "个节点的核心网络。\n")
    } else {
        g_core_network <- g
    }
    
    topo_stats <- data.frame(
      节点数 = vcount(g_core_network), 边数 = ecount(g_core_network), 平均度 = mean(degree(g_core_network)),
      网络密度 = edge_density(g_core_network), 聚集系数 = transitivity(g_core_network)
    )
    cat("高置信度核心网络拓扑统计:\n")
    print(topo_stats)
    
    # 保存高置信度网络数据
    write.csv(target_mapping, "results/tables/target_string_mapping_high_confidence.csv", row.names = FALSE)
    write.csv(target_interactions, "results/tables/ppi_interactions_high_confidence.csv", row.names = FALSE)
    write.csv(topo_stats, "results/tables/network_topology_stats_high_confidence.csv", row.names = FALSE)
    saveRDS(g_core_network, "results/network/ppi_network_high_confidence.rds")
    
    cat("✓ 高置信度核心网络构建和分析完成，结果已保存。\n")
  } else {
    cat("警告: 在高置信度阈值下，目标蛋白质之间未找到足够的相互作用\n")
  }
} else {
  cat("错误: 未能映射任何靶点到STRING数据库\n")
}

cat("高置信度网络构建分析完成！\n")