#!/usr/bin/env Rscript
# Functional Module Analysis Script for Lithospermum erythrorhizon
# Network clustering-based functional module identification

# Set random seeds for reproducibility
set.seed(42)  # For Louvain clustering algorithm

cat("=== Functional Module Analysis ===\n")
cat("Start time:", as.character(Sys.time()), "\n")

# Load required packages
suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(gridExtra)
})

# Set working directory to project root for relative paths
if(require(rstudioapi) && rstudioapi::isAvailable()) {
  project_root <- file.path(dirname(dirname(dirname(rstudioapi::getActiveDocumentContext()$path))))
  setwd(project_root)
}
output_dir <- "results/figures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 1. Load network data
if (file.exists("results/network/ppi_network.rds")) {
  g <- readRDS("results/network/ppi_network.rds")
  cat("✓ Successfully loaded PPI network, nodes:", vcount(g), ", edges:", ecount(g), "\n")
} else {
  stop("PPI network file does not exist, please run 02_complete_network_construction.R first")
}

# 2. Functional module detection
cat("Performing functional module analysis...\n")

# Ensure the network has correct gene symbol attributes
if (!"gene_symbol" %in% vertex_attr_names(g)) {
  # If no gene_symbol attribute, use name or name_display
  if ("name_display" %in% vertex_attr_names(g)) {
    V(g)$gene_symbol <- V(g)$name_display
  } else if ("name" %in% vertex_attr_names(g)) {
    V(g)$gene_symbol <- V(g)$name
  } else {
    # If neither exists, use vertex names
    V(g)$gene_symbol <- names(V(g))
  }
  cat("✓ Set gene symbol attributes\n")
}

# Use Louvain algorithm for community detection
communities <- cluster_louvain(g)
V(g)$community <- membership(communities)

cat("✓ Detected", max(V(g)$community), "functional modules\n")

# Calculate modularity coefficient
modularity_score <- modularity(communities)
cat("✓ Modularity score:", round(modularity_score, 3), "\n")

# 3. Module statistical analysis
module_stats <- data.frame(
  module = 1:max(V(g)$community),
  size = as.vector(table(V(g)$community)),
  stringsAsFactors = FALSE
)

module_stats$percentage <- round(module_stats$size / vcount(g) * 100, 1)

cat("\n=== Module Size Distribution ===\n")
for(i in 1:nrow(module_stats)) {
  cat(sprintf("Module %d: %d nodes (%.1f%%)\n", 
              module_stats$module[i], 
              module_stats$size[i], 
              module_stats$percentage[i]))
}

# 4. Visualize functional modules
set.seed(123)
module_colors <- rainbow(max(V(g)$community))
V(g)$module_color <- module_colors[V(g)$community]

# Use better layout algorithm
layout_coord <- layout_with_fr(g, niter = 500)

# Only show labels for major nodes to reduce density
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

# Add module boundaries
plot(communities, g, layout = layout_coord, add = TRUE, 
     col = alpha(module_colors, 0.2), 
     border = alpha(module_colors, 0.6),
     lwd = 2)

# Add legend
legend("topright", 
       legend = paste("Module", 1:max(V(g)$community)),
       fill = module_colors[1:max(V(g)$community)],
       cex = 1.2,
       bty = "n",
       title = "Functional Modules",
       title.cex = 1.3)

dev.off()

# 5. Generate module distribution visualization - improved bar chart and pie chart
# Bar chart version
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

# Improved pie chart version - with better colors and larger labels
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

# Save bar chart version
ggsave(file.path(output_dir, "module_distribution_bar.pdf"), 
       plot = p_module_bar, width = 12, height = 8, units = "in")
ggsave(file.path(output_dir, "module_distribution_bar.png"), 
       plot = p_module_bar, width = 12, height = 8, dpi = 300, units = "in")

# 6. Save module information and statistics
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Save module member information
module_info <- data.frame(
  node = names(V(g)),
  gene_symbol = V(g)$gene_symbol,
  module = V(g)$community,
  degree = degree(g),
  stringsAsFactors = FALSE
)

module_info <- module_info[order(module_info$module, -module_info$degree), ]
write.csv(module_info, "results/tables/functional_modules.csv", row.names = FALSE)

# Save module statistics
write.csv(module_stats, "results/tables/module_statistics.csv", row.names = FALSE)

# 7. Inter-module connectivity analysis
cat("\n=== Inter-module Connectivity Analysis ===\n")

# Calculate inter-module edges
edge_list <- as_edgelist(g, names = FALSE)  # Use new function
module_edges <- data.frame(
  from_module = V(g)$community[edge_list[,1]],
  to_module = V(g)$community[edge_list[,2]],
  stringsAsFactors = FALSE
)

# Count intra-module and inter-module connections
intra_module_edges <- sum(module_edges$from_module == module_edges$to_module)
inter_module_edges <- sum(module_edges$from_module != module_edges$to_module)

cat("Intra-module connections:", intra_module_edges, "\n")
cat("Inter-module connections:", inter_module_edges, "\n")
cat("Intra-module connection ratio:", round(intra_module_edges / ecount(g) * 100, 1), "%\n")

# 8. Module connectivity analysis and visualization
cat("\n=== Generating Module Connectivity Graph ===\n")

# Calculate inter-module connection matrix
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

# Combine plots
p_composition <- grid.arrange(p1_comp, p2_comp, p3_comp, p4_comp, 
                             ncol = 2, nrow = 2,
                             top = "Module Composition Analysis")

# Save module composition plot
ggsave(file.path(output_dir, "module_composition_plot.pdf"), 
       plot = p_composition, width = 12, height = 10, units = "in")
ggsave(file.path(output_dir, "module_composition_plot.png"), 
       plot = p_composition, width = 12, height = 10, dpi = 300, units = "in")

# Save module statistics data
write.csv(module_composition_stats, "results/tables/module_composition_stats.csv", row.names = FALSE)

# 10. Generate analysis summary
summary_info <- list(
  total_modules = max(V(g)$community),
  modularity_score = round(modularity_score, 3),
  largest_module = max(module_stats$size),
  smallest_module = min(module_stats$size),
  intra_module_edges = intra_module_edges,
  inter_module_edges = inter_module_edges,
  intra_module_percentage = round(intra_module_edges / ecount(g) * 100, 1)
)

cat("\n=== Functional Module Analysis Complete ===\n")
cat("✓ Detected modules:", summary_info$total_modules, "\n")
cat("✓ Modularity score:", summary_info$modularity_score, "\n")
cat("✓ Largest module size:", summary_info$largest_module, "nodes\n")
cat("✓ Intra-module connection ratio:", summary_info$intra_module_percentage, "%\n")
cat("✓ Saved plots: Figure4_Functional_Modules.png, module_distribution_pie.png\n")
cat("✓ Saved new plots: module_connectivity_plot.png, module_composition_plot.png\n")
cat("✓ Saved data: functional_modules.csv, module_statistics.csv, module_composition_stats.csv\n")
cat("Completion time:", as.character(Sys.time()), "\n")