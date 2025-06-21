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
  library(ggraph) # 新增ggraph包
})

# 设置目录路径
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

# 确保网络有正确的基因符号属性
if (!"gene_symbol" %in% vertex_attr_names(g)) {
  # 如果没有gene_symbol属性，使用name或name_display
  if ("name_display" %in% vertex_attr_names(g)) {
    V(g)$gene_symbol <- V(g)$name_display
  } else if ("name" %in% vertex_attr_names(g)) {
    V(g)$gene_symbol <- V(g)$name
  } else {
    # 如果都没有，使用节点名称
    V(g)$gene_symbol <- names(V(g))
  }
  cat("✓ 设置基因符号属性\n")
}

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

# 4. 可视化功能模块 - 方案一：分面布局 (Faceted Layout)
cat("正在生成分面布局的功能模块图 (Figure 4)...
")
set.seed(123)

# 使用更清晰、区分度更高的调色板
module_colors <- rainbow(max(V(g)$community), s = 0.8, v = 0.9)
V(g)$module_color <- module_colors[V(g)$community]

# 准备绘图数据
# 动态调整节点大小，基于度数
V(g)$size <- scales::rescale(degree(g), to = c(2.5, 10))
# 确定hub节点用于标签显示 (在每个模块内确定hub)
V(g)$show_label <- ""
for(comm_id in unique(V(g)$community)) {
    module_nodes <- which(V(g)$community == comm_id)
    module_degrees <- degree(g)[module_nodes]
    # 提高阈值，只显示每个模块最重要的节点
    hub_threshold <- quantile(module_degrees, 0.80) 
    # 确保至少显示一个节点，即使模块很小
    if (length(module_degrees) > 0 && all(module_degrees < hub_threshold)) {
        hub_threshold <- min(module_degrees)
    }
    hub_nodes <- module_nodes[module_degrees >= hub_threshold]
    V(g)$show_label[hub_nodes] <- V(g)$gene_symbol[hub_nodes]
}

# 确保模块ID是因子类型，以便分面
V(g)$community_factor <- as.factor(V(g)$community)

# 创建ggraph对象
p_modules_faceted <- ggraph(g, layout = 'fr') + 
  geom_edge_fan(aes(alpha = ..index..), color = 'grey65', width = 0.45) +
  geom_node_point(aes(color = community_factor, size = size), alpha = 0.9) +
  geom_node_text(aes(label = show_label), repel = TRUE, size = 3.5, fontface = 'bold', color = 'black', bg.color = "white", bg.r = 0.1) +
  # 核心：按社区（模块）进行分面
  facet_nodes(~community_factor, scales = 'free', ncol = 4) + # 4列布局
  # 使用美观的调色板
  scale_color_manual(values = module_colors, name = "Module") +
  scale_size_continuous(range = c(2.5, 10), name = "Node Degree") +
  # 主题设置
  theme_graph(base_family = 'sans', background = 'white', plot_margin = margin(10, 10, 10, 10)) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20, margin = margin(b = 5)), # 居中主标题
    plot.subtitle = element_text(hjust = 0.5, size = 14, margin = margin(b = 15)), # 居中副标题
    legend.position = 'bottom',
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(face = "bold", size = 14, margin = margin(t = 5, b = 5)) # 分面标签
  ) +
  labs(
    title = "Functional Modules in the PPI Network",
    subtitle = paste("Modularity Score:", round(modularity_score, 3), "| Nodes are faceted by community"),
    caption = "Layout: Fruchterman-Reingold | Node size corresponds to degree"
  )

# 保存新的分面图
output_filename_faceted <- file.path(output_dir, "Figure4_Functional_Modules_Faceted.png")
ggsave(output_filename_faceted, 
       plot = p_modules_faceted, 
       width = 18, height = 16, units = "in", dpi = 600, bg = "white")

cat(paste("✓ 成功保存分面模块图:", output_filename_faceted, "
"))

# 清理旧的Figure4文件，避免混淆
old_fig_path <- file.path(output_dir, "Figure4_Functional_Modules_Clean.png")
if (file.exists(old_fig_path)) {
  file.remove(old_fig_path)
  cat("✓ 已删除旧的Figure4文件:", old_fig_path, "
")
}


# 5. 生成模块分布可视化 - 改进的条形图和饼图
# 条形图版本
p_module_bar <- ggplot(module_stats, aes(x = reorder(factor(module), -size), y = size)) +
  geom_col(aes(fill = factor(module)), color = "white", linewidth = 0.8, alpha = 0.9) +
  geom_text(aes(label = paste0(size, " nodes\n(", percentage, "%)")), 
            vjust = -0.5, color = "black", fontface = "bold", size = 4.5) +
  scale_fill_brewer(name = "Module", palette = "Set2") +
  labs(
    title = "Distribution of Nodes Across Functional Modules",
    subtitle = paste("Total modules:", max(V(g)$community), "| Modularity:", round(modularity_score, 3)),
    caption = "Based on Louvain clustering algorithm",
    x = "Module (ordered by size)",
    y = "Number of Nodes"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(b = 15)),
    plot.caption = element_text(size = 11, hjust = 0.5, color = "grey60", margin = margin(t = 15)),
    legend.position = "bottom",
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# 改进的饼图版本 - 使用更好的颜色和更大的标签
p_module_pie <- ggplot(module_stats, aes(x = "", y = size, fill = factor(module))) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 1.2) +
  coord_polar("y", start = 0) +
  scale_fill_brewer(name = "Module", palette = "Set2") +
  labs(
    title = "Module Distribution (Pie Chart)",
    subtitle = paste("Total modules:", max(V(g)$community), "| Modularity:", round(modularity_score, 3)),
    caption = "Based on Louvain clustering algorithm"
  ) +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(b = 15)),
    plot.caption = element_text(size = 11, hjust = 0.5, color = "grey60", margin = margin(t = 15)),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  geom_text(aes(label = paste0(percentage, "%")), 
            position = position_stack(vjust = 0.5), 
            color = "white", fontface = "bold", size = 6, family = "sans")

ggsave(file.path(output_dir, "module_distribution_pie.pdf"), 
       plot = p_module_pie, width = 12, height = 10, units = "in")
ggsave(file.path(output_dir, "module_distribution_pie.png"), 
       plot = p_module_pie, width = 12, height = 10, dpi = 300, units = "in")

# 保存条形图版本
ggsave(file.path(output_dir, "module_distribution_bar.pdf"), 
       plot = p_module_bar, width = 12, height = 8, units = "in")
ggsave(file.path(output_dir, "module_distribution_bar.png"), 
       plot = p_module_bar, width = 12, height = 8, dpi = 300, units = "in")

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
edge_list <- as_edgelist(g, names = FALSE)  # 使用新的函数
module_edges <- data.frame(
  from_module = V(g)$community[edge_list[,1]],
  to_module = V(g)$community[edge_list[,2]],
  stringsAsFactors = FALSE
)

# 统计模块内和模块间连接
intra_module_edges <- sum(module_edges$from_module == module_edges$to_module)
inter_module_edges <- sum(module_edges$from_module != module_edges$to_module)

cat("模块内连接数:", intra_module_edges, "\n")
cat("模块间连接数:", inter_module_edges, "\n")
cat("模块内连接比例:", round(intra_module_edges / ecount(g) * 100, 1), "%\n")

# 8. 模块连接性分析和可视化
cat("\n=== 生成模块连接性图 ===\n")

# 计算模块间连接矩阵
module_connectivity_matrix <- matrix(0, nrow = max(V(g)$community), ncol = max(V(g)$community))
rownames(module_connectivity_matrix) <- paste("Module", 1:max(V(g)$community))
colnames(module_connectivity_matrix) <- paste("Module", 1:max(V(g)$community))

# 填充连接矩阵
for(i in 1:nrow(module_edges)) {
  from_mod <- module_edges$from_module[i]
  to_mod <- module_edges$to_module[i]
  module_connectivity_matrix[from_mod, to_mod] <- module_connectivity_matrix[from_mod, to_mod] + 1
  if(from_mod != to_mod) {
    module_connectivity_matrix[to_mod, from_mod] <- module_connectivity_matrix[to_mod, from_mod] + 1
  }
}

# 生成模块连接性热图
library(reshape2)
connectivity_df <- melt(module_connectivity_matrix)
colnames(connectivity_df) <- c("Module1", "Module2", "Connections")

p_connectivity <- ggplot(connectivity_df, aes(x = Module1, y = Module2, fill = Connections)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = Connections), color = "white", size = 6, fontface = "bold") +
  scale_fill_gradient(name = "Connections", 
                      low = "#3498DB", high = "#E74C3C",
                      guide = guide_colorbar(title.position = "top")) +
  labs(
    title = "Module Connectivity Matrix",
    subtitle = paste("Inter- and intra-module connections (Total modules:", max(V(g)$community), ")"),
    x = "Module",
    y = "Module",
    caption = "Values represent the number of edges between modules"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    plot.caption = element_text(size = 11, hjust = 0.5, color = "grey60"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white")
  ) +
  coord_fixed()

# 保存模块连接性图
ggsave(file.path(output_dir, "module_connectivity_plot.pdf"), 
       plot = p_connectivity, width = 10, height = 8, units = "in")
ggsave(file.path(output_dir, "module_connectivity_plot.png"), 
       plot = p_connectivity, width = 10, height = 8, dpi = 300, units = "in")

# 9. 模块组成分析图
cat("正在生成模块组成分析图...\n")

# 计算每个模块的连接度统计
module_composition_stats <- data.frame(
  module = 1:max(V(g)$community),
  size = as.vector(table(V(g)$community)),
  avg_degree = sapply(1:max(V(g)$community), function(m) {
    nodes_in_module <- which(V(g)$community == m)
    mean(degree(g)[nodes_in_module])
  }),
  max_degree = sapply(1:max(V(g)$community), function(m) {
    nodes_in_module <- which(V(g)$community == m)
    max(degree(g)[nodes_in_module])
  }),
  internal_edges = sapply(1:max(V(g)$community), function(m) {
    sum(module_edges$from_module == m & module_edges$to_module == m)
  }),
  external_edges = sapply(1:max(V(g)$community), function(m) {
    sum((module_edges$from_module == m & module_edges$to_module != m) |
        (module_edges$to_module == m & module_edges$from_module != m))
  })
)

module_composition_stats$density <- with(module_composition_stats, 
                                        2 * internal_edges / (size * (size - 1)))
module_composition_stats$external_ratio <- with(module_composition_stats, 
                                               external_edges / (internal_edges + external_edges))

# 生成组合图表
p1_comp <- ggplot(module_composition_stats, aes(x = factor(module), y = size)) +
  geom_col(fill = "#3498DB", alpha = 0.8, color = "black") +
  geom_text(aes(label = size), vjust = -0.5, size = 4, fontface = "bold") +
  labs(title = "Module Size Distribution", x = "Module", y = "Number of Nodes") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold")
  )

p2_comp <- ggplot(module_composition_stats, aes(x = factor(module), y = avg_degree)) +
  geom_col(fill = "#E74C3C", alpha = 0.8, color = "black") +
  geom_text(aes(label = round(avg_degree, 1)), vjust = -0.5, size = 4, fontface = "bold") +
  labs(title = "Average Degree per Module", x = "Module", y = "Average Degree") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold")
  )

p3_comp <- ggplot(module_composition_stats, aes(x = factor(module), y = density)) +
  geom_col(fill = "#2ECC71", alpha = 0.8, color = "black") +
  geom_text(aes(label = round(density, 2)), vjust = -0.5, size = 4, fontface = "bold") +
  labs(title = "Module Density", x = "Module", y = "Internal Density") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold")
  )

p4_comp <- ggplot(module_composition_stats, aes(x = factor(module), y = external_ratio)) +
  geom_col(fill = "#F39C12", alpha = 0.8, color = "black") +
  geom_text(aes(label = round(external_ratio, 2)), vjust = -0.5, size = 4, fontface = "bold") +
  labs(title = "External Connection Ratio", x = "Module", y = "External/Total Ratio") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold")
  )

# 组合图表
p_composition <- grid.arrange(p1_comp, p2_comp, p3_comp, p4_comp, 
                             ncol = 2, nrow = 2,
                             top = "Module Composition Analysis")

# 保存模块组成图
ggsave(file.path(output_dir, "module_composition_plot.pdf"), 
       plot = p_composition, width = 12, height = 10, units = "in")
ggsave(file.path(output_dir, "module_composition_plot.png"), 
       plot = p_composition, width = 12, height = 10, dpi = 300, units = "in")

# 保存模块统计数据
write.csv(module_composition_stats, "results/tables/module_composition_stats.csv", row.names = FALSE)

# 10. 生成分析摘要
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
cat("✓ 保存图表: Figure4_Functional_Modules_Faceted.png, module_distribution_pie.png\n")
cat("✓ 保存新图表: module_connectivity_plot.png, module_composition_plot.png\n")
cat("✓ 保存数据: functional_modules.csv, module_statistics.csv, module_composition_stats.csv\n")
cat("完成时间:", as.character(Sys.time()), "\n")