# Lithospermum erythrorhizon Network Pharmacology Analysis - Network Visualization
# PPI network visualization based on STRING database

# Load required packages
suppressMessages({
  library(igraph)
  library(ggplot2)
  library(ggraph)
  library(ggrepel)
  library(dplyr)
  library(RColorBrewer)
  library(gridExtra)
  library(corrplot)
  library(pheatmap)
})

# Set working directory to project root
if(require(rstudioapi) && rstudioapi::isAvailable()) {
  # Running in RStudio
  project_root <- file.path(dirname(dirname(dirname(rstudioapi::getActiveDocumentContext()$path))))
  setwd(project_root)
} else {
  # Running from command line - assume we're already in project root
  # or set to a known location
  if(basename(getwd()) == "R") {
    setwd("../../")
  } else if(basename(getwd()) == "scripts") {
    setwd("../")
  }
  # If running from project root, do nothing
}

# Create output directories
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

cat("Starting network visualization analysis...\n")

# 1. Load network data
if (file.exists("results/network/ppi_network.rds")) {
  g <- readRDS("results/network/ppi_network.rds")
  cat("Successfully loaded PPI network with", vcount(g), "nodes and", ecount(g), "edges\n")
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

# 2. 网络基本可视化
cat("正在生成基础网络图...\n")

# 设置节点颜色（简化为单一类型，因为都是紫草靶点）
node_colors <- c("紫草靶点" = "#FF6B6B", "其他" = "#95A5A6")

# 确保所有节点都有node_type属性
if (!"node_type" %in% vertex_attr_names(g)) {
  V(g)$node_type <- "紫草靶点"
}

V(g)$color <- node_colors[V(g)$node_type]

# 设置节点大小 (基于度中心性)
V(g)$size <- scales::rescale(degree(g), to = c(5, 20))

# 设置边的权重
E(g)$width <- scales::rescale(E(g)$combined_score, to = c(0.5, 3))

# 设置节点标签（使用基因符号）
if (!is.null(target_mapping)) {
  # 创建STRING ID到基因符号的映射
  id_to_name <- setNames(target_mapping$target_name, target_mapping$string_id)
  
  # 为每个节点设置基因符号标签
  V(g)$gene_symbol <- sapply(names(V(g)), function(x) {
    if (x %in% names(id_to_name)) {
      return(id_to_name[x])
    } else {
      return(V(g)[x]$name_display)
    }
  })
} else {
  V(g)$gene_symbol <- V(g)$name_display
}

# 基础网络布局
set.seed(123)
layout_fr <- layout_with_fr(g)

# 生成网络图
png("results/figures/Figure1_PPI_Network.png", width = 14, height = 12, units = "in", res = 300)
par(mar = c(2, 2, 4, 6))
plot(g, 
     layout = layout_fr,
     vertex.label = V(g)$gene_symbol,
     vertex.label.cex = 0.8,
     vertex.label.color = "black",
     vertex.label.font = 2,
     vertex.frame.color = "white",
     vertex.frame.width = 2,
     edge.color = alpha("gray60", 0.7),
     edge.curved = 0.1,
     main = "Lithospermum erythrorhizon PPI Network\n(Based on STRING Database v12.0)",
     sub = paste("Nodes:", vcount(g), "| Edges:", ecount(g), "| Confidence ≥ 400"),
     cex.main = 1.5,
     cex.sub = 1.2)

# 添加图例
legend("topright", 
       legend = c("Lithospermum Targets", paste("Network Density:", round(edge_density(g), 3))),
       col = c(node_colors[1], "black"),
       pch = c(19, NA),
       cex = 1.0,
       title = "Network Properties",
       bty = "n")
dev.off()

# 3. 网络拓扑分析可视化
cat("正在生成拓扑分析图...\n")

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

# 度分布图
p1 <- ggplot(topo_df, aes(x = degree)) +
  geom_histogram(bins = 15, fill = "#FF6B6B", alpha = 0.7, color = "white") +
  labs(title = "Degree Distribution", 
       x = "Degree", 
       y = "Count",
       subtitle = paste("Mean degree:", round(mean(degree_centrality), 2))) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))

# 中心性散点图
p2 <- ggplot(topo_df, aes(x = degree, y = betweenness, size = eigenvector)) +
  geom_point(color = "#FF6B6B", alpha = 0.8) +
  geom_text_repel(aes(label = ifelse(degree > quantile(degree, 0.8) | 
                                   betweenness > quantile(betweenness, 0.8), 
                                   gene_symbol, "")),
                  size = 3, max.overlaps = 15) +
  labs(title = "Centrality Analysis", 
       x = "Degree Centrality", 
       y = "Betweenness Centrality",
       size = "Eigenvector\nCentrality") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))

# 中心性相关性热图数据
centrality_cor <- cor(topo_df[, c("degree", "betweenness", "closeness", "eigenvector")])

png("results/figures/Figure2_Network_Topology.png", width = 16, height = 8, units = "in", res = 300)
grid.arrange(p1, p2, ncol = 2)
dev.off()

# 中心性相关性热图
png("results/figures/Figure3_Centrality_Correlation.png", width = 10, height = 8, units = "in", res = 300)
corrplot(centrality_cor, method = "color", type = "upper", 
         addCoef.col = "black", tl.col = "black", tl.srt = 45,
         title = "Centrality Measures Correlation",
         mar = c(0,0,2,0))
dev.off()

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

# 5. 功能模块分析
cat("正在进行功能模块分析...\n")

# 社区检测
communities <- cluster_louvain(g)
V(g)$community <- membership(communities)

cat("检测到", max(V(g)$community), "个功能模块\n")

# 模块化系数
modularity_score <- modularity(communities)
cat("模块化系数:", round(modularity_score, 3), "\n")

# 可视化功能模块
set.seed(123)
module_colors <- rainbow(max(V(g)$community))
V(g)$module_color <- module_colors[V(g)$community]

png("results/figures/Figure4_Functional_Modules.png", width = 14, height = 12, units = "in", res = 300)
par(mar = c(2, 2, 4, 6))
plot(g, 
     layout = layout_fr,
     vertex.color = V(g)$module_color,
     vertex.label = V(g)$gene_symbol,
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     vertex.label.font = 2,
     vertex.size = 10,
     vertex.frame.color = "white",
     edge.color = alpha("gray60", 0.6),
     main = paste("Functional Modules\n(Modularity =", round(modularity_score, 3), ")"),
     cex.main = 1.5)

# 添加模块边界
plot(communities, g, layout = layout_fr, add = TRUE, 
     col = alpha(module_colors, 0.3), border = module_colors)
dev.off()

# 6. 生成分析摘要
cat("正在生成分析摘要...\n")

summary_stats <- data.frame(
  指标 = c("节点总数", "边总数", "平均度", "网络密度", "聚集系数", 
           "平均路径长度", "网络直径", "连通分量数", "功能模块数", "模块化系数"),
  值 = c(vcount(g), ecount(g), round(mean(degree(g)), 2), 
         round(edge_density(g), 4), round(transitivity(g), 3),
         ifelse(is_connected(g), round(mean_distance(g), 2), "不连通"),
         ifelse(is_connected(g), diameter(g), "不连通"),
         components(g)$no, max(V(g)$community), round(modularity_score, 3)),
  stringsAsFactors = FALSE
)

write.csv(summary_stats, "results/tables/network_summary.csv", row.names = FALSE)

# 保存模块信息
module_info <- data.frame(
  node = names(V(g)),
  gene_symbol = V(g)$gene_symbol,
  module = V(g)$community,
  stringsAsFactors = FALSE
)

module_info <- module_info[order(module_info$module, -key_nodes$degree), ]
write.csv(module_info, "results/tables/functional_modules.csv", row.names = FALSE)

# 保存可视化图例说明
legends <- c(
  "Figure1_PPI_Network.png: 基于STRING数据库v12.0的蛋白质相互作用网络，节点大小表示度中心性，边宽度表示相互作用置信度",
  "Figure2_Network_Topology.png: 网络拓扑分析，包括度分布和中心性分析散点图",
  "Figure3_Centrality_Correlation.png: 不同中心性指标之间的相关性热图",
  "Figure4_Functional_Modules.png: 基于Louvain算法的功能模块检测结果，不同颜色代表不同功能模块"
)

writeLines(legends, "results/figures/figure_legends.txt")

cat("网络可视化分析完成！\n")
cat("生成的图表文件:\n")
cat("- Figure1_PPI_Network.png (主要PPI网络图)\n")
cat("- Figure2_Network_Topology.png (拓扑分析图)\n")
cat("- Figure3_Centrality_Correlation.png (中心性相关性)\n")
cat("- Figure4_Functional_Modules.png (功能模块图)\n")
cat("分析结果表格已保存到 results/tables/ 目录\n") 