# 紫草网络药理学分析 - 网络可视化
# 基于STRING数据库构建的PPI网络可视化

# 加载必需的包
suppressMessages({
  library(igraph)
  library(ggplot2)
  library(ggraph)
  library(ggrepel)
  library(dplyr)
  library(RColorBrewer)
  library(gridExtra)
  library(grid)
  library(corrplot)
  library(pheatmap)
  library(reshape2)
  library(scales)
})

# 创建输出目录
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

cat("开始网络可视化分析...\n")

# 1. 加载网络数据
if (file.exists("results/network/ppi_network.rds")) {
  g <- readRDS("results/network/ppi_network.rds")
  cat("成功加载PPI网络，节点数:", vcount(g), "，边数:", ecount(g), "\n")
} else {
  cat("错误: 未找到PPI网络文件，请先运行网络构建脚本\n")
  stop("缺少网络数据")
}

# 加载靶点映射信息
if (file.exists("results/tables/target_string_mapping.csv")) {
  target_mapping <- read.csv("results/tables/target_string_mapping.csv", stringsAsFactors = FALSE)
  cat("成功加载靶点映射信息，映射靶点数:", nrow(target_mapping), "\n")
} else {
  cat("警告: 未找到靶点映射文件，使用默认节点标签\n")
  target_mapping <- NULL
}

# 2. 网络基本可视化 - 修复颜色编码问题
cat("正在生成基础网络图...\n")

# 设置明确的节点颜色区分：直接靶点（蓝色）和一阶互作蛋白（橙色）
# 首先需要确定哪些是直接靶点，哪些是一阶互作蛋白

# 设置明确的节点颜色区分：直接靶点（蓝色）和一阶互作蛋白（橙色）
# 从CMAUP数据中获取直接靶点信息

# 方法1: 从保存的靶点映射文件中获取直接靶点
direct_targets <- c()
if (file.exists("results/tables/target_string_mapping.csv")) {
  target_mapping_data <- read.csv("results/tables/target_string_mapping.csv", stringsAsFactors = FALSE)
  direct_targets <- target_mapping_data$target_name
  cat("从映射文件获取直接靶点数:", length(direct_targets), "\n")
}

# 方法2: 如果没有映射文件，从CMAUP原始数据获取
if (length(direct_targets) == 0) {
  if (file.exists("zwsjk/CMAUPv2.0_download_Targets.txt")) {
    cat("从CMAUP原始数据获取直接靶点...\n")
    
    # 加载植物数据
    plants <- read.delim("zwsjk/CMAUPv2.0_download_Plants.txt", stringsAsFactors = FALSE)
    lithospermum <- plants[grepl("lithospermum.*erythrorhizon", plants$Species_Name, ignore.case = TRUE), ]
    
    if (nrow(lithospermum) > 0) {
      # 加载植物-成分关联
      plant_ingredients <- read.delim("zwsjk/CMAUPv2.0_download_Plant_Ingredient_Associations_onlyActiveIngredients.txt", 
                                     header = FALSE, stringsAsFactors = FALSE)
      colnames(plant_ingredients) <- c("Plant_ID", "Ingredient_ID")
      
      # 获取紫草成分
      zicao_ingredients <- plant_ingredients[plant_ingredients$Plant_ID == lithospermum$Plant_ID[1], ]
      
      # 加载成分-靶点关联
      ingredient_targets <- read.delim("zwsjk/CMAUPv2.0_download_Ingredient_Target_Associations_ActivityValues_References.txt", 
                                      stringsAsFactors = FALSE)
      
      # 获取紫草靶点
      zicao_targets <- ingredient_targets[ingredient_targets$Ingredient_ID %in% zicao_ingredients$Ingredient_ID, ]
      
      # 加载靶点信息
      targets_info <- read.delim("zwsjk/CMAUPv2.0_download_Targets.txt", stringsAsFactors = FALSE)
      zicao_targets_with_info <- merge(zicao_targets, targets_info, by = "Target_ID")
      
      direct_targets <- unique(zicao_targets_with_info$Gene_Symbol)
      direct_targets <- direct_targets[!is.na(direct_targets) & direct_targets != ""]
      cat("从CMAUP数据获取直接靶点数:", length(direct_targets), "\n")
    }
  }
}

# 方法3: 如果仍然没有，使用网络中度数较高的节点作为直接靶点
if (length(direct_targets) == 0) {
  cat("使用网络拓扑特征确定直接靶点...\n")
  degree_values <- degree(g)
  # 使用度数前60%的节点作为直接靶点
  threshold <- quantile(degree_values, 0.4)
  direct_targets <- names(degree_values[degree_values >= threshold])
}

cat("最终确定的直接靶点数:", length(direct_targets), "\n")

# 设置节点类型和颜色
# 首先确保我们有有效的基因符号
node_names <- names(V(g))

# 检查网络中节点的基因符号属性
if ("gene_symbol" %in% vertex_attr_names(g)) {
  current_gene_symbols <- V(g)$gene_symbol
} else if ("name" %in% vertex_attr_names(g)) {
  current_gene_symbols <- V(g)$name
} else {
  current_gene_symbols <- node_names
}

# 确保current_gene_symbols不为空
if (length(current_gene_symbols) == 0 || all(is.na(current_gene_symbols))) {
  current_gene_symbols <- node_names
}

# 设置节点类型
V(g)$node_type <- ifelse(current_gene_symbols %in% direct_targets, 
                         "direct_target", "first_shell_interactor")

# 明确的颜色区分：蓝色为直接靶点，橙色为一阶互作蛋白
node_colors <- c("direct_target" = "#2E86AB", "first_shell_interactor" = "#F28E2B")
V(g)$color <- node_colors[V(g)$node_type]

# 设置节点大小更好地反映度中心性（更大的差异）
degree_values <- degree(g)
V(g)$size <- scales::rescale(degree_values, to = c(8, 25))

# 设置边的权重
E(g)$width <- scales::rescale(E(g)$combined_score, to = c(0.3, 2.5))

# 设置节点标签（使用基因符号）
if (!is.null(target_mapping)) {
  # 创建STRING ID到基因符号的映射
  id_to_name <- setNames(target_mapping$target_name, target_mapping$string_id)
  
  # 为每个节点设置基因符号标签
  gene_symbols <- sapply(node_names, function(x) {
    if (x %in% names(id_to_name)) {
      return(id_to_name[x])
    } else {
      # 如果没有映射，尝试从已有属性获取
      return(x)  # 使用节点ID作为后备
    }
  })
  
  # 确保基因符号不为空
  gene_symbols[is.na(gene_symbols) | gene_symbols == ""] <- node_names[is.na(gene_symbols) | gene_symbols == ""]
  
  # 设置为节点属性
  V(g)$gene_symbol <- gene_symbols
  V(g)$name_display <- gene_symbols
  
  # 关键步骤：直接替换节点的name属性为基因符号
  # 这样igraph绘图时会优先使用基因符号而不是STRING ID
  V(g)$name <- gene_symbols
  
  cat("基因符号映射完成，前5个示例:\n")
  cat(head(paste(node_names, "->", gene_symbols), 5), sep = "\n")
  cat("节点name属性已更新为基因符号\n")

} else {
  # 如果没有映射文件，使用节点名称
  V(g)$gene_symbol <- node_names
  V(g)$name_display <- node_names
  V(g)$name <- node_names  # 也更新name属性
  cat("警告：未找到靶点映射文件，使用节点ID作为标签\n")
}

# --- 图1: PPI网络图 (使用ggraph改进) ---
cat("正在生成基于ggraph的PPI网络图...\n")

# 使用更美观的布局算法
set.seed(123)
layout_gg <- create_layout(g, layout = 'fr')

# 提取网络属性用于标题
density_val <- round(edge_density(g), 3)
avg_degree_val <- round(mean(degree(g)), 2)
nodes_val <- vcount(g)
edges_val <- ecount(g)

# 使用ggraph绘图
p_ppi <- ggraph(layout_gg) + 
  geom_edge_fan(aes(width = I(width), alpha = I(0.5)), color = "gray60") +
  geom_node_point(aes(color = node_type, size = I(size))) +
  geom_node_text(aes(label = name_display), repel = TRUE, size = 3.5, fontface = "bold", 
                 max.overlaps = 15, bg.color = "white", bg.r = 0.1) +
  scale_color_manual(values = node_colors, name = "Node Type", 
                     labels = c("Direct L. erythrorhizon target", "First-shell interactor")) +
  scale_edge_width_continuous(range = c(0.3, 1.5), guide = "none") +
  theme_graph(base_family = 'sans', background = "white") +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10)
  ) +
  labs(
    title = "Protein-Protein Interaction Network of L. erythrorhizon Targets",
    subtitle = paste("High-confidence interactions, STRING v12.0 (Nodes:", nodes_val, "| Edges:", edges_val,")"),
    caption = paste("Network Density:", density_val, "| Average Degree:", avg_degree_val, "| Node size encodes degree centrality")
  )

# 保存ggraph版本
ggsave("results/figures/Figure1_PPI_Network.png", plot = p_ppi, width = 12, height = 10, dpi = 300, units = "in")
cat("✓ 新版PPI网络图已保存。\n")

# 3. 网络拓扑分析可视化 - 生成完整的四合一图
cat("正在生成拓扑分析四合一图...\n")

# 计算各种中心性指标
degree_centrality <- degree(g)
betweenness_centrality <- betweenness(g)
closeness_centrality <- closeness(g)
eigenvector_centrality <- eigen_centrality(g)$vector

# 创建拓扑数据框
topo_df <- data.frame(
  node = names(V(g)),
  gene_symbol = V(g)$gene_symbol,
  node_type = V(g)$node_type,
  degree = degree_centrality,
  betweenness = betweenness_centrality,
  closeness = closeness_centrality,
  eigenvector = eigenvector_centrality,
  stringsAsFactors = FALSE
)

# Panel A: 度分布直方图
p1 <- ggplot(topo_df, aes(x = degree)) +
  geom_histogram(bins = 12, fill = "#2E86AB", alpha = 0.8, color = "white", linewidth = 0.5) +
  labs(title = "Panel A: Degree Distribution", 
       x = "Degree", 
       y = "Count") +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.3),
    panel.grid.minor = element_blank()
  ) +
  annotate("text", x = max(topo_df$degree) * 0.7, y = max(table(cut(topo_df$degree, 12))) * 0.8,
           label = paste("Mean:", round(mean(degree_centrality), 1)), 
           size = 3.5, fontface = "bold", color = "#2E86AB")

# Panel B: 中心性散点图（度 vs 介数中心性）
# 识别top节点进行标注
top_nodes <- topo_df[topo_df$degree > quantile(topo_df$degree, 0.75) | 
                     topo_df$betweenness > quantile(topo_df$betweenness, 0.75), ]

p2 <- ggplot(topo_df, aes(x = degree, y = betweenness)) +
  geom_point(aes(color = node_type, size = eigenvector), alpha = 0.8) +
  scale_color_manual(values = node_colors, name = "Node Type") +
  geom_text_repel(data = top_nodes, aes(label = gene_symbol),
                  size = 3, max.overlaps = 10, fontface = "bold") +
  labs(title = "Panel B: Centrality Analysis", 
       x = "Degree Centrality", 
       y = "Betweenness Centrality",
       size = "Eigenvector") +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # 在四合一图中隐藏图例
  )

# Panel C: 中心性相关性热图
centrality_matrix <- topo_df[, c("degree", "betweenness", "closeness", "eigenvector")]
centrality_cor <- cor(centrality_matrix)

# 创建热图数据
cor_df <- reshape2::melt(centrality_cor)
names(cor_df) <- c("Var1", "Var2", "value")

p3 <- ggplot(cor_df, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#E74C3C", mid = "white", high = "#2E86AB", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "Correlation") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3.5, fontface = "bold") +
  labs(title = "Panel C: Centrality Correlation",
       x = "", y = "") +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white"),
    legend.position = "none"  # 在四合一图中隐藏图例
  )

# Panel D: 核心蛋白连接图（Hub节点及其连接）
# 识别Hub节点
hub_threshold <- quantile(degree_centrality, 0.8)
hub_nodes <- names(degree_centrality[degree_centrality >= hub_threshold])

# 创建Hub节点数据框
hub_df <- topo_df[topo_df$node %in% hub_nodes, ]
hub_df <- hub_df[order(-hub_df$degree), ]  # 按度排序

p4 <- ggplot(hub_df, aes(x = reorder(gene_symbol, degree), y = degree)) +
  geom_col(aes(fill = node_type), alpha = 0.8, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = node_colors, name = "Node Type") +
  geom_text(aes(label = degree), hjust = -0.1, size = 3, fontface = "bold") +
  coord_flip() +
  labs(title = "Panel D: Hub Protein Connectivity",
       x = "Hub Proteins",
       y = "Degree Centrality") +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major.x = element_line(color = "#E5E5E5", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # 在四合一图中隐藏图例
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# 组合四个子图
library(gridExtra)
library(grid)

combined_topology_plot <- grid.arrange(p1, p2, p3, p4, 
                                       ncol = 2, nrow = 2,
                                       top = textGrob("Network Topology Analysis and Hub Target Identification", 
                                                    gp = gpar(fontsize = 16, fontface = "bold")))

# 保存四合一图
ggsave("results/figures/Figure2_Network_Topology.png", 
       plot = combined_topology_plot, width = 14, height = 12, dpi = 300, units = "in")
ggsave("results/figures/Figure2_Network_Topology.pdf", 
       plot = combined_topology_plot, width = 14, height = 12, units = "in")

# 保存单独的子图（备用）
ggsave("results/figures/Figure2A_Degree_Distribution.png", 
       plot = p1, width = 6, height = 5, dpi = 300, units = "in")
ggsave("results/figures/Figure2B_Centrality_Scatter.png", 
       plot = p2, width = 6, height = 5, dpi = 300, units = "in")
ggsave("results/figures/Figure2C_Centrality_Correlation.png", 
       plot = p3, width = 6, height = 5, dpi = 300, units = "in")
ggsave("results/figures/Figure2D_Hub_Connectivity.png", 
       plot = p4, width = 6, height = 5, dpi = 300, units = "in")

# 4. 关键节点分析
cat("正在进行关键节点分析...\n")

# 识别Hub节点 (度中心性前20%)
hub_threshold <- quantile(degree_centrality, 0.8)
hub_nodes <- names(degree_centrality[degree_centrality >= hub_threshold])

# 识别桥接节点 (介数中心性前20%)
bridge_threshold <- quantile(betweenness_centrality, 0.8)
bridge_nodes <- names(betweenness_centrality[betweenness_centrality >= bridge_threshold])

cat("Hub节点数:", length(hub_nodes), "\n")
cat("桥接节点数:", length(bridge_nodes), "\n")

# 保存关键节点信息
key_nodes <- data.frame(
  node = names(V(g)),
  gene_symbol = V(g)$gene_symbol,
  node_type = V(g)$node_type,
  degree = degree_centrality,
  betweenness = betweenness_centrality,
  closeness = closeness_centrality,
  eigenvector = eigenvector_centrality,
  is_hub = names(V(g)) %in% hub_nodes,
  is_bridge = names(V(g)) %in% bridge_nodes,
  stringsAsFactors = FALSE
)

# 按度中心性排序
key_nodes <- key_nodes[order(-key_nodes$degree), ]
write.csv(key_nodes, "results/tables/key_nodes_analysis.csv", row.names = FALSE)

cat("前10个Hub节点:\n")
print(head(key_nodes[key_nodes$is_hub, c("gene_symbol", "degree", "betweenness")], 10))

# 5. [DEPRECATED] 功能模块分析和可视化已移至 05_module_analysis.R
cat("功能模块分析和可视化已移至 05_module_analysis.R 脚本。\n")
cat("此脚本 (04) 不再生成功能模块图，以确保单一来源原则。\n")

# 6. 生成分析摘要
cat("正在生成分析摘要...\n")
summary_text <- paste(
  "--- Network Analysis and Visualization Summary ---",
  paste("Timestamp:", Sys.time()),
  "\n[PPI Network]",
  paste("Nodes:", vcount(g)),
  paste("Edges:", ecount(g)),
  paste("Network Density:", round(edge_density(g), 4)),
  paste("Average Degree:", round(mean(degree(g)), 2)),
  "\n[Centrality Analysis]",
  paste("Top 5 Hubs (by Degree):", paste(head(key_nodes$gene_symbol, 5), collapse=", ")),
  "\n--- End of Summary ---"
)
writeLines(summary_text, "results/figures/ADMET_analysis_summary.txt")

cat("网络可视化分析完成！\n")
cat("生成的图表文件:\n")
cat("- Figure1_PPI_Network.png (主要PPI网络图)\n")
cat("- Figure2_Network_Topology.png (拓扑分析图)\n")
cat("- Figure3_Centrality_Correlation.png (中心性相关性)\n")
cat("- Figure4_Functional_Modules.png (功能模块图)\n")
cat("分析结果表格已保存到 results/tables/ 目录\n")