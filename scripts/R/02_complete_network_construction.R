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

# Set working directory for reproducibility
if (!file.exists("data") && file.exists("../../zwsjk")) {
  setwd("../..")
}

# Create output directories (relative paths)
dir.create("results/network", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

cat("Starting complete network construction analysis...\n")
cat("Working directory:", getwd(), "\n")

# 1. Load CMAUP data
cat("Loading CMAUP data from relative path...\n")

# Load plant data (using relative path)
plants <- read_tsv("../../zwsjk/CMAUPv2.0_download_Plants.txt", show_col_types = FALSE)

# 找到紫草
lithospermum <- plants %>%
  filter(str_detect(tolower(Species_Name), "lithospermum") & 
         str_detect(tolower(Species_Name), "erythrorhizon"))

cat("找到紫草，Plant ID:", lithospermum$Plant_ID, "\n")

# Load plant-ingredient associations (using relative path)
plant_ingredients <- read_tsv("../../zwsjk/CMAUPv2.0_download_Plant_Ingredient_Associations_onlyActiveIngredients.txt", 
                             col_names = c("Plant_ID", "Ingredient_ID"), 
                             show_col_types = FALSE)

# Filter ingredients for Lithospermum
zicao_ingredient_ids <- plant_ingredients %>%
  filter(Plant_ID == lithospermum$Plant_ID) %>%
  pull(Ingredient_ID)

cat("Number of Lithospermum-related ingredients:", length(zicao_ingredient_ids), "\n")

# Load ingredient-target associations (using relative path)
ingredient_targets <- read_tsv("../../zwsjk/CMAUPv2.0_download_Ingredient_Target_Associations_ActivityValues_References.txt", 
                              show_col_types = FALSE)

# 筛选紫草成分的靶点
zicao_targets_data <- ingredient_targets %>%
  filter(Ingredient_ID %in% zicao_ingredient_ids)

cat("紫草成分-靶点关联数:", nrow(zicao_targets_data), "\n")

# Load target detailed information (using relative path)
targets <- read_tsv("../../zwsjk/CMAUPv2.0_download_Targets.txt", show_col_types = FALSE)

# 合并靶点信息
zicao_targets_with_info <- zicao_targets_data %>%
  left_join(targets, by = "Target_ID") %>%
  filter(!is.na(Gene_Symbol))

cat("有效紫草靶点数:", length(unique(zicao_targets_with_info$Gene_Symbol)), "\n")

# 准备疾病相关靶点（假设有疾病数据或使用已知银屑病相关基因）
# 这里我们先用紫草靶点作为主要分析对象
all_targets <- unique(zicao_targets_with_info$Gene_Symbol)
all_targets <- all_targets[!is.na(all_targets) & all_targets != ""]

cat("总靶点数:", length(all_targets), "\n")

# 3. Load STRING database files (using relative paths)
cat("Loading STRING database...\n")

# Load protein information (check multiple possible locations)
string_paths <- c("../../data/string_db", "data/string_db")
string_db_path <- NULL
for (path in string_paths) {
  if (file.exists(file.path(path, "9606.protein.info.v12.0.txt.gz"))) {
    string_db_path <- path
    break
  }
}

if (is.null(string_db_path)) {
  stop("STRING database files not found. Please run scripts/download_string_db.sh first.")
}

cat("Using STRING database from:", string_db_path, "\n")

# Load protein information
protein_info <- fread(file.path(string_db_path, "9606.protein.info.v12.0.txt.gz"))
cat("STRING protein info loaded:", nrow(protein_info), "records\n")

# Load protein aliases
protein_aliases <- fread(file.path(string_db_path, "9606.protein.aliases.v12.0.txt.gz"))
cat("STRING protein aliases loaded:", nrow(protein_aliases), "records\n")

# Check if protein interaction file exists
if (file.exists(file.path(string_db_path, "9606.protein.links.v12.0.txt.gz"))) {
  protein_links <- fread(file.path(string_db_path, "9606.protein.links.v12.0.txt.gz"))
  cat("STRING protein interactions loaded:", nrow(protein_links), "records\n")
  
  # Filter high-quality interactions (score >= 400)
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
    
    # 计算网络拓扑参数
    topo_stats <- data.frame(
      节点数 = vcount(g),
      边数 = ecount(g),
      平均度 = mean(degree(g)),
      网络密度 = edge_density(g),
      聚集系数 = transitivity(g),
      平均路径长度 = ifelse(is.connected(g), average.path.length(g), NA),
      直径 = ifelse(is.connected(g), diameter(g), NA),
      连通分量数 = components(g)$no
    )
    
    cat("网络拓扑统计:\n")
    print(topo_stats)
    
    # 保存网络数据
    write.csv(target_mapping, "results/tables/target_string_mapping.csv", row.names = FALSE)
    write.csv(target_interactions, "results/tables/ppi_interactions.csv", row.names = FALSE)
    write.csv(topo_stats, "results/tables/network_topology_stats.csv", row.names = FALSE)
    
    # 保存网络对象
    saveRDS(g, "results/network/ppi_network.rds")
    
    cat("网络构建完成！\n")
    cat("结果已保存到 results/ 目录\n")
    
  } else {
    cat("警告: 在目标蛋白质之间未找到足够的相互作用\n")
  }
  
} else {
  cat("错误: 未能映射任何靶点到STRING数据库\n")
}

cat("网络构建分析完成！\n") 