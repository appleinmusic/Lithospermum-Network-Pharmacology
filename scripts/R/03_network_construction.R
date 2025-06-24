# 紫草网络药理学分析 - 完整网络构建
# 使用STRING数据库构建蛋白质相互作用网络

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

cat("开始完整网络构建分析...\n")

# ==============================================================================
# 1. 加载通过ADMET筛选的化合物
# 这是关键的修正步骤，确保网络分析基于筛选后的化合物
# ==============================================================================
cat("正在加载通过ADMET筛选的化合物...\n")
admet_filtered_compounds_path <- "results/admet_analysis_data.csv"
if (!file.exists(admet_filtered_compounds_path)) {
  stop("错误: ADMET筛选结果文件 '", admet_filtered_compounds_path, "' 不存在。\n请先运行 05_ADMET_analysis.R 脚本。")
}
admet_filtered_compounds <- read_csv(admet_filtered_compounds_path, show_col_types = FALSE)
filtered_ingredient_ids <- unique(admet_filtered_compounds$np_id)
cat("✓ 加载了 ", length(filtered_ingredient_ids), " 个通过ADMET筛选的化合物。\n")


# ==============================================================================
# 2. 根据筛选后的化合物确定靶点
# ==============================================================================
cat("正在根据筛选后的化合物确定靶点...\n")

# 加载成分-靶点关联数据
ingredient_targets_path <- "zwsjk/CMAUPv2.0_download_Ingredient_Target_Associations_ActivityValues_References.txt"
ingredient_targets <- read_tsv(ingredient_targets_path, show_col_types = FALSE)

# 筛选紫草成分的靶点
zicao_targets_data <- ingredient_targets %>%
  filter(Ingredient_ID %in% filtered_ingredient_ids)

cat("✓ 筛选后化合物的靶点关联数:", nrow(zicao_targets_data), "\n")

# 加载靶点详细信息
targets_path <- "zwsjk/CMAUPv2.0_download_Targets.txt"
targets <- read_tsv(targets_path, show_col_types = FALSE)

# 合并靶点信息
zicao_targets_with_info <- zicao_targets_data %>%
  left_join(targets, by = "Target_ID") %>%
  filter(!is.na(Gene_Symbol))

# 准备最终的靶点列表
all_targets <- unique(zicao_targets_with_info$Gene_Symbol)
all_targets <- all_targets[!is.na(all_targets) & all_targets != ""]

cat("✓ 最终用于网络构建的独特靶点数:", length(all_targets), "\n")


# 3. 加载STRING数据库文件
cat("正在加载STRING数据库...\n")

# 加载蛋白质信息
protein_info <- fread("data/string_db/9606.protein.info.v12.0.txt.gz")
cat("STRING蛋白质信息加载完成，共", nrow(protein_info), "条记录\n")

# 加载蛋白质别名
protein_aliases <- fread("data/string_db/9606.protein.aliases.v12.0.txt.gz")
cat("STRING蛋白质别名加载完成，共", nrow(protein_aliases), "条记录\n")

# 检查是否存在蛋白质相互作用文件
if (file.exists("data/string_db/9606.protein.links.v12.0.txt.gz")) {
  protein_links <- fread("data/string_db/9606.protein.links.v12.0.txt.gz")
  cat("STRING蛋白质相互作用数据加载完成，共", nrow(protein_links), "条记录\n")
  
  # 筛选高质量相互作用 (score >= 400)
  high_quality_links <- protein_links[protein_links$combined_score >= 400, ]
  cat("高质量相互作用 (score>=400):", nrow(high_quality_links), "条\n")
  
} else {
  cat("警告: 未找到蛋白质相互作用数据文件\n")
  cat("请下载 9606.protein.links.v12.0.txt.gz 并放置到 data/string_db/ 目录\n")
  cat("下载链接: https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz\n")
  stop("缺少必需的STRING数据文件")
}

# 4. 靶点名称映射到STRING ID
cat("正在进行靶点名称映射...\n")

# 检查列名
cat("STRING别名文件列名:", colnames(protein_aliases), "\n")

# 创建映射函数
map_to_string_id <- function(target_names, aliases_df, protein_info_df) {
  # 确保使用正确的列名
  protein_id_col <- colnames(aliases_df)[1]  # 第一列是蛋白质ID列
  alias_col <- "alias"
  
  cat("使用的蛋白质ID列名:", protein_id_col, "\n")
  
  mapped_data <- data.frame(
    target_name = character(0),
    string_id = character(0),
    source = character(0),
    stringsAsFactors = FALSE
  )
  
  for (target in target_names) {
    if (is.na(target) || target == "") next
    
    found <- FALSE
    
    # 1. 首先尝试在蛋白质信息文件的preferred_name中查找
    protein_matches <- protein_info_df[protein_info_df$preferred_name == target, ]
    
    if (nrow(protein_matches) > 0) {
      protein_id <- protein_matches[[colnames(protein_info_df)[1]]][1]
      if (length(protein_id) > 0 && !is.na(protein_id) && nchar(as.character(protein_id)) > 0) {
        new_row <- data.frame(
          target_name = target,
          string_id = as.character(protein_id),
          source = "preferred_name",
          stringsAsFactors = FALSE
        )
        mapped_data <- rbind(mapped_data, new_row)
        cat("成功映射:", target, "->", as.character(protein_id), "(来源: preferred_name)\n")
        found <- TRUE
      }
    }
    
    # 2. 如果在preferred_name中没找到，再尝试别名文件
    if (!found) {
      # 精确匹配基因符号
      matches <- aliases_df[aliases_df[[alias_col]] == target, ]
      
      if (nrow(matches) == 0) {
        # 尝试忽略大小写的精确匹配
        matches <- aliases_df[toupper(aliases_df[[alias_col]]) == toupper(target), ]
      }
      
      if (nrow(matches) == 0) {
        # 尝试部分匹配
        matches <- aliases_df[grepl(paste0("^", target, "$"), aliases_df[[alias_col]], ignore.case = TRUE), ]
      }
      
      if (nrow(matches) > 0) {
        # 检查蛋白质ID是否存在且不为NA
        protein_id <- matches[[protein_id_col]][1]
        
        if (length(protein_id) > 0 && !is.na(protein_id) && nchar(as.character(protein_id)) > 0) {
          new_row <- data.frame(
            target_name = target,
            string_id = as.character(protein_id),
            source = paste0("alias:", matches$source[1]),
            stringsAsFactors = FALSE
          )
          mapped_data <- rbind(mapped_data, new_row)
          cat("成功映射:", target, "->", as.character(protein_id), "(来源: alias -", matches$source[1], ")\n")
          found <- TRUE
        } else {
          cat("靶点", target, "在别名中找到匹配但蛋白质ID为空，匹配数:", nrow(matches), "\n")
        }
      }
    }
    
    if (!found) {
      cat("未找到映射:", target, "\n")
    }
  }
  
  return(mapped_data)
}

# 执行映射
target_mapping <- map_to_string_id(all_targets, protein_aliases, protein_info)
cat("成功映射靶点数:", nrow(target_mapping), "/", length(all_targets), "\n")

# --- 关键修正：报告并处理未映射的靶点 ---
unmapped_targets <- setdiff(all_targets, target_mapping$target_name)
if (length(unmapped_targets) > 0) {
    cat("警告：以下", length(unmapped_targets), "个靶点未能映射到STRING ID，将不会包含在网络中:\n")
    for (target in unmapped_targets) {
        cat("  - ", target, "\n")
    }
}

# 显示成功映射的靶点
if (nrow(target_mapping) > 0) {
  cat("成功映射的靶点:\n")
  print(target_mapping)
  
  # 保存映射结果
  write_csv(target_mapping, "results/tables/target_string_mapping.csv")
  cat("✓ 靶点映射结果已保存到 'results/tables/target_string_mapping.csv'\n")
}

# 5. 构建PPI网络
cat("正在构建蛋白质相互作用网络...\n")
string_ids <- unique(target_mapping$string_id)

# 筛选相互作用
ppi_interactions <- high_quality_links %>%
  filter(protein1 %in% string_ids & protein2 %in% string_ids)

cat("目标蛋白质间相互作用数:", nrow(ppi_interactions), "\n")

# 创建igraph对象
ppi_network <- graph_from_data_frame(
  d = ppi_interactions[, c("protein1", "protein2")],
  directed = FALSE
)

# ==============================================================================
# 6. 【科学性修正】提取最大连通分量 (LCC)
#    这是确保网络分析稳健性和聚焦核心功能的标准步骤。
# ==============================================================================
cat("\n正在提取最大连通分量 (LCC)...\n")
components <- decompose(ppi_network)
if (length(components) > 0) {
  # 找到最大的连通分量
  largest_comp_index <- which.max(sapply(components, vcount))
  core_network <- components[[largest_comp_index]]
  
  # 报告提取结果
  original_node_count <- vcount(ppi_network)
  core_node_count <- vcount(core_network)
  
  if (original_node_count > core_node_count) {
    cat("信息：从", original_node_count, "个节点的网络中提取出", core_node_count, "个节点的核心网络。\n")
    
    # 识别并报告被排除的节点
    original_nodes <- V(ppi_network)$name
    core_nodes <- V(core_network)$name
    excluded_nodes_ids <- setdiff(original_nodes, core_nodes)
    
    # 将 STRING ID 转换回基因符号以便阅读
    excluded_nodes_df <- protein_info %>%
      filter(`#string_protein_id` %in% excluded_nodes_ids) %>%
      select(preferred_name)
    
    excluded_node_names <- excluded_nodes_df$preferred_name
    
    if (length(excluded_node_names) > 0) {
      cat("以下", length(excluded_node_names), "个节点是孤立节点，已被排除:\n")
      for (name in excluded_node_names) {
        cat("  - ", name, "\n")
      }
    }
  } else {
    cat("✓ 网络本身已是单一连通分量，无需提取。\n")
  }
} else {
  core_network <- make_empty_graph() # 如果没有组件则创建空图
  cat("警告：初始网络为空或没有连通组件。\n")
}

# 7. 保存核心网络对象和拓扑统计
if (vcount(core_network) > 0) {
  # 计算并保存网络拓扑统计
  network_stats <- data.frame(
    nodes = vcount(core_network),
    edges = ecount(core_network),
    avg_degree = mean(degree(core_network)),
    density = graph.density(core_network),
    clustering_coefficient = transitivity(core_network, type = "global"),
    avg_path_length = mean_distance(core_network),
    diameter = diameter(core_network),
    components = length(decompose(core_network))
  )
  cat("核心网络拓扑统计:\n")
  print(network_stats)
  
  # 保存统计信息
  write_csv(network_stats, "results/tables/network_topology_stats.csv")
  
  # 保存igraph对象
  saveRDS(core_network, file = "results/network/ppi_network.rds")
  
  # 保存边列表
  ppi_edge_list <- as_data_frame(core_network, what = "edges")
  write_csv(ppi_edge_list, "results/tables/ppi_interactions.csv")
  
  cat("✓ 核心网络构建和分析完成，结果已保存。\n")
} else {
  cat("最终网络为空，不保存任何结果。\n")
}

cat("网络构建分析完成！\n") 