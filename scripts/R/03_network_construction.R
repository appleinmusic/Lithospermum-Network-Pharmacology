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

# **关键修正**: 使用已处理的数据文件，确保与02_ADMET_filtering.R输出一致
filtered_targets_path <- "data/processed/lithospermum_targets.tsv"
if (!file.exists(filtered_targets_path)) {
  stop("错误: 已处理的靶点文件 '", filtered_targets_path, "' 不存在。\n请先运行 01_data_preparation.R 和 02_ADMET_filtering.R 脚本。")
}

# 加载已经处理和过滤的目标数据（这些数据已经与ADMET筛选后的化合物对应）
zicao_targets_filtered <- read_tsv(filtered_targets_path, show_col_types = FALSE)

cat("✓ 筛选后化合物的靶点关联数:", nrow(zicao_targets_filtered), "\n")

# 准备最终的靶点列表
all_targets <- unique(zicao_targets_filtered$Gene_Symbol)
all_targets <- all_targets[!is.na(all_targets) & all_targets != ""]

cat("✓ 最终用于网络构建的独特靶点数:", length(all_targets), "\n")

# **重要**: 保存筛选后的目标数据用于验证
filtered_targets_for_network <- zicao_targets_filtered %>%
  select(Ingredient_ID, Gene_Symbol) %>%
  distinct()
write_tsv(filtered_targets_for_network, "data/processed/lithospermum_targets_for_network.tsv")


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
}

# 5. 构建PPI网络
cat("正在构建蛋白质相互作用网络...\n")

if (nrow(target_mapping) > 0) {
  # 获取映射靶点的相互作用
  target_string_ids <- target_mapping$string_id
  
  # 筛选目标蛋白质之间的相互作用
  target_interactions <- high_quality_links[
    high_quality_links$protein1 %in% target_string_ids & 
    high_quality_links$protein2 %in% target_string_ids, 
  ]
  
  cat("目标蛋白质间相互作用数:", nrow(target_interactions), "\n")
  
  if (nrow(target_interactions) > 0) {
    # 创建网络图
    g <- graph_from_data_frame(
      target_interactions[, c("protein1", "protein2", "combined_score")],
      directed = FALSE
    )
    
    # 添加节点属性
    # 映射STRING ID回靶点名称
    id_to_name <- setNames(target_mapping$target_name, target_mapping$string_id)
    
    V(g)$name_display <- ifelse(names(V(g)) %in% names(id_to_name), 
                               id_to_name[names(V(g))], 
                               names(V(g)))
    
    # 标记节点类型
    V(g)$node_type <- "紫草靶点"
    
    # --- 关键修正：提取最大连通分量 (LCC) ---
    cat("\n正在提取最大连通分量 (LCC)...\n")
    components_list <- decompose(g)
    if (length(components_list) > 0) {
        component_sizes <- sapply(components_list, vcount)
        largest_component_index <- which.max(component_sizes)
        g_core_network <- components_list[[largest_component_index]]
        
        core_nodes_string_ids <- V(g_core_network)$name
        all_nodes_string_ids <- V(g)$name
        excluded_nodes_string_ids <- setdiff(all_nodes_string_ids, core_nodes_string_ids)
        
        if (length(excluded_nodes_string_ids) > 0) {
            excluded_node_names <- V(g)$name_display[match(excluded_nodes_string_ids, V(g)$name)]
            cat("信息：从", vcount(g), "个节点的网络中提取出", vcount(g_core_network), "个节点的核心网络。\n")
            cat("以下", length(excluded_node_names), "个节点是孤立节点，已被排除:\n")
            for(node in excluded_node_names) {
                cat("  - ", node, "\n")
            }
        }
    } else {
        g_core_network <- g # 如果网络为空或只有一个分量
    }
    
    # 后续分析应全部基于 g_core_network
    
    # 计算网络拓扑参数 (基于核心网络)
    topo_stats <- data.frame(
      节点数 = vcount(g_core_network),
      边数 = ecount(g_core_network),
      平均度 = mean(degree(g_core_network)),
      网络密度 = edge_density(g_core_network),
      聚集系数 = transitivity(g_core_network),
      平均路径长度 = ifelse(is_connected(g_core_network), mean_distance(g_core_network), NA),
      直径 = ifelse(is_connected(g_core_network), diameter(g_core_network), NA),
      连通分量数 = components(g_core_network)$no
    )
    
    cat("核心网络拓扑统计:\n")
    print(topo_stats)
    
    # 保存网络数据
    write.csv(target_mapping, "results/tables/target_string_mapping.csv", row.names = FALSE)
    write.csv(target_interactions, "results/tables/ppi_interactions.csv", row.names = FALSE)
    write.csv(topo_stats, "results/tables/network_topology_stats.csv", row.names = FALSE)
    
    # 保存网络对象 (保存修正后的核心网络)
    saveRDS(g_core_network, "results/network/ppi_network.rds")
    
    cat("✓ 核心网络构建和分析完成，结果已保存。\n")
    
  } else {
    cat("警告: 在目标蛋白质之间未找到足够的相互作用\n")
  }
  
} else {
  cat("错误: 未能映射任何靶点到STRING数据库\n")
}

cat("网络构建分析完成！\n")