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
  library(ggforce) # 新增，用于geom_mark_hull
  library(concaveman) # geom_mark_hull的依赖包
  library(RColorBrewer)
  library(tidygraph) # 【结构性修正】引入tidygraph包，确保数据处理的稳健性
})

# 设置目录路径
# setwd("/Users/lgmoon/Desktop/zdhky")  # 相对路径版本，注释掉硬编码路径
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

# 4. 可视化功能模块 - 最终方案：单一、清晰的ggraph图
cat("正在生成最终版的功能模块图 (Figure4_Functional_Modules.png)...\n")
set.seed(42) # 保证布局可重复

# 使用更适合发表的调色板
module_colors <- RColorBrewer::brewer.pal(n = max(4, max(V(g)$community)), name = "Set2")
V(g)$module_color <- module_colors[V(g)$community]

# 准备绘图数据
V(g)$size <- scales::rescale(degree(g), to = c(3, 12))
V(g)$degree <- degree(g) # 将度中心性存为属性

# --- 智能标签逻辑 ---
# 目标：只为每个模块中度中心性最高的2个节点添加标签
# 显式创建数据框以确保所有列都存在
node_df <- data.frame(
  name = V(g)$name,
  gene_symbol = V(g)$gene_symbol,
  community = V(g)$community,
  degree = V(g)$degree
)

label_nodes <- node_df %>%
  group_by(community) %>%
  top_n(2, degree) %>%
  pull(name)

V(g)$show_label <- ifelse(V(g)$name %in% label_nodes, V(g)$gene_symbol, "")

# 创建布局
layout_fr <- create_layout(g, layout = 'fr')

# --- ggraph绘图 ---
p_modules_final <- ggraph(layout_fr) +
  # 绘制模块的凸包背景
  geom_mark_hull(
    aes(x, y, group = community, fill = factor(community), label = paste("Module", community)),
    concavity = 4,
    expand = unit(3, 'mm'),
    alpha = 0.2,
    show.legend = FALSE,
    label.buffer = unit(10, 'mm'),
    label.fontsize = 10
  ) +
  # 绘制边
  geom_edge_fan(aes(alpha = ..index..), color = 'grey50', width = 0.5) +
  # 绘制节点
  geom_node_point(aes(color = factor(community), size = size), alpha = 0.8) +
  # 绘制智能标签
  geom_node_text(aes(label = show_label), repel = TRUE, size = 3.5, fontface = 'bold', 
                 bg.color = "white", bg.r = 0.1, max.overlaps = 15) +
  
  # 设置颜色、大小和主题
  scale_fill_manual(values = module_colors) +
  scale_color_manual(values = module_colors, name = "Functional Module") +
  scale_size_continuous(range = c(3, 12), name = "Node Degree") +
  theme_graph(base_family = 'sans', background = 'white') +
  guides(color = guide_legend(override.aes = list(size=5))) + # 增大图例中的点
  labs(
    title = "Functional Modules in the PPI Network",
    subtitle = paste("Identified by Louvain clustering (Modularity:", round(modularity_score, 3), ")"),
    caption = "Node size corresponds to degree. Only top 2 hub nodes per module are labeled."
  )

# 保存最终的图表
output_filename_final <- file.path(output_dir, "Figure4_Functional_Modules.png")
ggsave(output_filename_final,
       plot = p_modules_final,
       width = 14, height = 12, units = "in", dpi = 300, bg = "white")

cat(paste("✓ 成功保存最终版功能模块图:", output_filename_final, "\n"))

# 清理旧的、可能混淆的图表文件
other_fig_paths <- file.path(output_dir, c("Figure4_Functional_Modules_Faceted.png", "Figure4_Functional_Modules_Enhanced.png"))
for (f_path in other_fig_paths) {
    if (file.exists(f_path)) {
        file.remove(f_path)
        cat("✓ 已删除旧图:", f_path, "\n")
    }
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

# 7. 模块间连接分析 和 Top 3 模块可视化 (最终修复版)
cat("\n=== Top 3 功能模块深度分析 ===\n")

# 找出最大的3个模块的ID
top3_modules <- module_stats %>%
  top_n(3, size) %>%
  pull(module)

cat("--> 已确定Top 3模块 (ID):", paste(top3_modules, collapse=", "), "\n")

# 【结构性修正】使用 tidygraph 进行稳健的子图操作
# 从主图 g 创建 tidygraph 对象，并加入社区信息
g_tidy <- as_tbl_graph(g, directed = FALSE) %>%
  activate(nodes) %>%
  mutate(community = V(g)$community)

# 筛选出只包含 Top 3 模块的节点，并创建子图
g_top3_tidy <- g_tidy %>%
  filter(community %in% top3_modules)

# 获取子图的节点和边数
v_count_top3 <- g_top3_tidy %>% activate(nodes) %>% as_tibble() %>% nrow()
e_count_top3 <- g_top3_tidy %>% activate(edges) %>% as_tibble() %>% nrow()

# 【错误修复 V3】正确的模块化得分计算方法
# 在子图上重新进行社区检测
g_top3_igraph <- as.igraph(g_top3_tidy)
communities_top3 <- cluster_louvain(g_top3_igraph)

# 使用 igraph 的 modularity 函数计算模块化得分
modularity_top3 <- modularity(communities_top3)

# 将社区信息添加回tidygraph对象
g_top3_tidy <- g_top3_tidy %>%
    mutate(community = as.factor(membership(communities_top3)))

cat("✓ 已创建Top 3模块子图，包含", v_count_top3, "个节点和", e_count_top3, "条边。新模块化得分:", round(modularity_top3, 3), "\n")


# --- 可视化Top 3模块网络 ---
cat("--> 正在生成Top 3模块网络图 (Figure4_Top3_Modules.png)...\n")
set.seed(123) # 使用不同的种子以获得可能更优的布局

# 检查子图是否为空
if (v_count_top3 == 0) {
  cat("! Top 3 模块子图为空，跳过可视化。\n")
} else {
  # 创建布局
  layout_top3 <- create_layout(g_top3_tidy, layout = 'fr')

  # 准备智能标签
  node_df_top3 <- as_tibble(activate(g_top3_tidy, "nodes"))
  label_nodes_top3 <- node_df_top3 %>%
    group_by(community) %>%
    top_n(3, degree) %>%
    pull(name)
  
  g_top3_tidy <- g_top3_tidy %>%
    mutate(show_label = ifelse(name %in% label_nodes_top3, gene_symbol, ""))

  p_top3_modules <- ggraph(g_top3_tidy, layout = layout_top3) +
    geom_mark_hull(
      aes(x, y, group = community, fill = factor(community)),
      concavity = 4, expand = unit(3, 'mm'), alpha = 0.2, show.legend = FALSE
    ) +
    # 【错误修复】不再引用不存在的 'weight' 属性。使用固定的边样式，与主图一致。
    geom_edge_fan(color = 'grey60', width = 0.6, alpha = 0.7) +
    geom_node_point(aes(color = factor(community), size = degree), alpha = 0.9) +
    geom_node_text(aes(label = show_label), repel = TRUE, size = 4, fontface = 'bold',
                   bg.color = "white", bg.r = 0.1, max.overlaps = 20) +
    scale_fill_manual(values = module_colors) +
    scale_color_manual(values = module_colors, name = "Functional Module") +
    scale_size_continuous(name = "Node Degree", range=c(4,15)) +
    theme_graph(base_family = 'sans', background = 'white') +
    labs(
      title = "Analysis of Top 3 Functional Modules",
      subtitle = paste("Network of", v_count_top3, "nodes and", e_count_top3, "edges."),
      caption = "Node size corresponds to degree. Top 3 hub nodes per module are labeled."
    )

  output_filename_top3 <- file.path(output_dir, "Figure4_Top3_Modules.png")
  ggsave(output_filename_top3,
         plot = p_top3_modules,
         width = 14, height = 12, units = "in", dpi = 300, bg = "white")

  cat("✓ 成功保存Top 3模块网络图:", output_filename_top3, "\n")
}


cat("\n=== 功能模块分析完成 ===\n")
cat("✓ 检测模块数:", max(V(g)$community), "\n")
cat("✓ 模块化系数:", round(modularity_score, 3), "\n")
cat("✓ 最大模块大小:", max(module_stats$size), "个节点\n")
cat("✓ 保存图表: Figure4_Functional_Modules.png, Figure4_Top3_Modules.png, module_distribution_*.png\n")
cat("✓ 保存数据: functional_modules.csv, module_statistics.csv\n")
cat("完成时间:", as.character(Sys.time()), "\n")