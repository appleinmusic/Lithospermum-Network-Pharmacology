#!/usr/bin/env Rscript
# 功能模块分析脚本 - Lithospermum erythrorhizon
# 基于网络聚类识别功能模块

cat("=== 功能模块分析 ===\n")
cat("开始时间:", as.character(Sys.time()), "\n")

# 加载必需包
suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(gridExtra)
})

# 设置目录路径
if (!file.exists("data") && file.exists("../../zwsjk")) {
  setwd("../..")
}
output_dir <- "results/figures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 1. 加载网络数据
if (file.exists("results/network/ppi_network.rds")) {
  g <- readRDS("results/network/ppi_network.rds")
  cat("✓ 成功加载PPI网络，节点数:", vcount(g), "，边数:", ecount(g), "\n")
} else {
  stop("PPI网络文件不存在，请先运行02_complete_network_construction.R脚本")
}

# 2. 功能模块检测
cat("正在进行功能模块分析...\n")

# 使用Louvain算法进行社区检测
communities <- cluster_louvain(g)
V(g)$community <- membership(communities)

cat("✓ 检测到", max(V(g)$community), "个功能模块\n")

# 计算模块化系数
modularity_score <- modularity(communities)
cat("✓ 模块化系数:", round(modularity_score, 3), "\n")

# 3. 模块统计分析
module_stats <- data.frame(
  module = 1:max(V(g)$community),
  size = as.vector(table(V(g)$community)),
  stringsAsFactors = FALSE
)

module_stats$percentage <- round(module_stats$size / vcount(g) * 100, 1)

cat("\n=== 模块大小分布 ===\n")
for(i in 1:nrow(module_stats)) {
  cat(sprintf("模块 %d: %d 个节点 (%.1f%%)\n", 
              module_stats$module[i], 
              module_stats$size[i], 
              module_stats$percentage[i]))
}

# 4. 可视化功能模块
set.seed(123)
module_colors <- rainbow(max(V(g)$community))
V(g)$module_color <- module_colors[V(g)$community]

# 使用更好的布局算法
layout_coord <- layout_with_fr(g, niter = 500)

# 只显示主要节点的标签，减少密集度
hub_threshold <- quantile(degree(g), 0.7)  
V(g)$show_label <- ifelse(degree(g) >= hub_threshold, V(g)$gene_symbol, "")

png(file.path(output_dir, "Figure4_Functional_Modules.png"), 
    width = 16, height = 14, units = "in", res = 300)
par(mar = c(2, 2, 4, 8), bg = "white")

plot(g, 
     layout = layout_coord,
     vertex.color = V(g)$module_color,
     vertex.label = V(g)$show_label,
     vertex.label.cex = 1.0,
     vertex.label.color = "black",
     vertex.label.font = 2,
     vertex.size = pmax(8, degree(g) * 0.8),
     vertex.frame.color = "white",
     vertex.frame.width = 2,
     edge.color = alpha("gray50", 0.4),
     edge.width = 0.8,
     main = paste("Functional Modules in PPI Network\n(Modularity Score =", round(modularity_score, 3), ")"),
     cex.main = 1.8,
     font.main = 2)

# 添加模块边界
plot(communities, g, layout = layout_coord, add = TRUE, 
     col = alpha(module_colors, 0.2), 
     border = alpha(module_colors, 0.6),
     lwd = 2)

# 添加图例
legend("topright", 
       legend = paste("Module", 1:max(V(g)$community)),
       fill = module_colors[1:max(V(g)$community)],
       cex = 1.2,
       bty = "n",
       title = "Functional Modules",
       title.cex = 1.3)

dev.off()

# 5. 生成模块分布饼图
p_pie <- ggplot(module_stats, aes(x = "", y = size, fill = factor(module))) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = module_colors, name = "Module") +
  labs(
    title = "Distribution of Nodes Across Functional Modules",
    subtitle = paste("Total modules:", max(V(g)$community), "| Modularity:", round(modularity_score, 3)),
    caption = "Based on Louvain clustering algorithm"
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5, color = "grey60"),
    legend.position = "right",
    plot.background = element_rect(fill = "white")
  ) +
  geom_text(aes(label = ifelse(percentage >= 10, paste0(percentage, "%"), "")), 
            position = position_stack(vjust = 0.5), 
            color = "black", fontface = "bold", size = 4)

ggsave(file.path(output_dir, "module_distribution_pie.pdf"), 
       plot = p_pie, width = 10, height = 8, units = "in")
ggsave(file.path(output_dir, "module_distribution_pie.png"), 
       plot = p_pie, width = 10, height = 8, dpi = 300, units = "in")

# 6. 保存模块信息和统计
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# 保存模块成员信息
module_info <- data.frame(
  node = names(V(g)),
  gene_symbol = V(g)$gene_symbol,
  module = V(g)$community,
  degree = degree(g),
  stringsAsFactors = FALSE
)

module_info <- module_info[order(module_info$module, -module_info$degree), ]
write.csv(module_info, "results/tables/functional_modules.csv", row.names = FALSE)

# 保存模块统计
write.csv(module_stats, "results/tables/module_statistics.csv", row.names = FALSE)

# 7. 模块间连接分析
cat("\n=== 模块间连接分析 ===\n")

# 计算模块间边数
module_edges <- data.frame(
  from_module = V(g)$community[get.edgelist(g)[,1]],
  to_module = V(g)$community[get.edgelist(g)[,2]],
  stringsAsFactors = FALSE
)

# 统计模块内和模块间连接
intra_module_edges <- sum(module_edges$from_module == module_edges$to_module)
inter_module_edges <- sum(module_edges$from_module != module_edges$to_module)

cat("模块内连接数:", intra_module_edges, "\n")
cat("模块间连接数:", inter_module_edges, "\n")
cat("模块内连接比例:", round(intra_module_edges / ecount(g) * 100, 1), "%\n")

# 8. 生成分析摘要
summary_info <- list(
  total_modules = max(V(g)$community),
  modularity_score = round(modularity_score, 3),
  largest_module = max(module_stats$size),
  smallest_module = min(module_stats$size),
  intra_module_edges = intra_module_edges,
  inter_module_edges = inter_module_edges,
  intra_module_percentage = round(intra_module_edges / ecount(g) * 100, 1)
)

cat("\n=== 功能模块分析完成 ===\n")
cat("✓ 检测模块数:", summary_info$total_modules, "\n")
cat("✓ 模块化系数:", summary_info$modularity_score, "\n")
cat("✓ 最大模块大小:", summary_info$largest_module, "个节点\n")
cat("✓ 模块内连接比例:", summary_info$intra_module_percentage, "%\n")
cat("✓ 保存图表: Figure4_Functional_Modules.png, module_distribution_pie.png\n")
cat("✓ 保存数据: functional_modules.csv, module_statistics.csv\n")
cat("完成时间:", as.character(Sys.time()), "\n")