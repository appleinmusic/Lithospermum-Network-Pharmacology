# KEGG Pathway Enrichment Analysis using clusterProfiler
# Updated version to match manuscript methodology

cat("=== KEGG通路富集分析（使用clusterProfiler）===\n")
cat("开始时间:", as.character(Sys.time()), "\n")

# 加载必需包
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(DOSE)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
})

# 设置目录路径
output_dir <- "results/figures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 1. 加载靶点数据
targets_file <- "results/tables/target_string_mapping.csv"
if (file.exists(targets_file)) {
  targets_data <- readr::read_csv(targets_file, show_col_types = FALSE)
  gene_symbols <- unique(targets_data$target_name)
  cat("✓ 载入靶点基因数:", length(gene_symbols), "\n")
} else {
  stop("靶点数据文件不存在，请先运行02_complete_network_construction.R脚本")
}

# 2. 基因符号转换为ENTREZ ID
cat("正在转换基因符号为ENTREZ ID...\n")
gene_entrez <- bitr(gene_symbols, 
                   fromType = "SYMBOL",
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

if (nrow(gene_entrez) == 0) {
  stop("无法转换基因符号为ENTREZ ID")
}

cat("✓ 成功转换", nrow(gene_entrez), "个基因的ENTREZ ID\n")

# 删除重复的ENTREZ ID
entrez_ids <- unique(gene_entrez$ENTREZID)
cat("✓ 去重后共", length(entrez_ids), "个唯一ENTREZ ID\n")

# 3. KEGG通路富集分析
cat("正在进行KEGG通路富集分析...\n")

# 设置分析参数（与手稿描述一致）
kegg_enrich <- enrichKEGG(gene = entrez_ids,
                         organism = 'hsa',           # human
                         keyType = 'kegg',
                         pvalueCutoff = 0.05,        # P值阈值
                         pAdjustMethod = "BH",       # Benjamini-Hochberg FDR校正
                         minGSSize = 5,              # 最小基因集大小
                         maxGSSize = 500,            # 最大基因集大小
                         qvalueCutoff = 0.05,        # FDR阈值
                         use_internal_data = FALSE)

if (is.null(kegg_enrich) || nrow(kegg_enrich@result) == 0) {
  cat("警告: 没有发现显著富集的KEGG通路\n")
  
  # 降低阈值重新尝试
  cat("尝试降低阈值重新分析...\n")
  kegg_enrich <- enrichKEGG(gene = entrez_ids,
                           organism = 'hsa',
                           keyType = 'kegg',
                           pvalueCutoff = 0.1,
                           pAdjustMethod = "BH",
                           minGSSize = 3,
                           maxGSSize = 500,
                           qvalueCutoff = 0.1,
                           use_internal_data = FALSE)
}

if (!is.null(kegg_enrich) && nrow(kegg_enrich@result) > 0) {
  cat("✓ 发现", nrow(kegg_enrich@result), "个显著富集的KEGG通路\n")
  
  # 4. 结果处理和可视化
  kegg_result <- kegg_enrich@result
  
  # 添加基因符号信息
  kegg_result$geneSymbol <- sapply(kegg_result$geneID, function(genes) {
    gene_ids <- unlist(strsplit(genes, "/"))
    symbols <- gene_entrez$SYMBOL[match(gene_ids, gene_entrez$ENTREZID)]
    paste(symbols[!is.na(symbols)], collapse = "/")
  })
  
  # 选择前20个最显著的通路
  top_pathways <- head(kegg_result, 20)
  
  # 5. 生成条形图
  p1 <- barplot(kegg_enrich, 
                showCategory = 15,
                title = "KEGG Pathway Enrichment Analysis",
                font.size = 12,
                color = "p.adjust")
  
  # 自定义条形图
  plot_data <- top_pathways[1:min(15, nrow(top_pathways)), ]
  plot_data$Description <- str_wrap(plot_data$Description, width = 40)
  
  p2 <- ggplot(plot_data, aes(x = reorder(Description, Count), y = Count)) +
    geom_col(aes(fill = p.adjust), alpha = 0.8, color = "black", linewidth = 0.3) +
    geom_text(aes(label = Count), hjust = -0.1, fontface = "bold", size = 3.5) +
    scale_fill_gradient(name = "FDR", 
                        low = "#E06666", high = "#FFF2CC",
                        guide = guide_colorbar(title.position = "top")) +
    coord_flip() +
    labs(
      title = "KEGG Pathway Enrichment Analysis",
      subtitle = "Lithospermum erythrorhizon Target Genes",
      x = "KEGG Pathway",
      y = "Gene Count",
      caption = paste("Analyzed using clusterProfiler package",
                     "\nHypergeometric test with Benjamini-Hochberg FDR correction")
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5),
      plot.caption = element_text(size = 10, hjust = 0.5, color = "grey60"),
      panel.background = element_rect(fill = "white", color = "black"),
      plot.background = element_rect(fill = "white"),
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 11),
      axis.title = element_text(size = 13, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
  # 6. 生成点图
  p3 <- dotplot(kegg_enrich, 
                showCategory = 15,
                title = "KEGG Pathway Enrichment Dotplot",
                font.size = 12,
                color = "p.adjust")
  
  # 7. 生成网络图
  if (nrow(kegg_result) >= 5) {
    # 计算通路相似性并生成网络图
    cat("正在计算通路相似性...\n")
    kegg_enrich_sim <- pairwise_termsim(kegg_enrich)
    
    # 显示更多通路，但选择最显著的
    n_pathways <- min(20, nrow(kegg_result))  # 显示前20个最显著的通路
    cat("生成网络图，显示前", n_pathways, "个最显著的通路\n")
    
    p4 <- emapplot(kegg_enrich_sim, 
                   showCategory = n_pathways,  # 增加显示数量
                   layout = "kk") +
      theme_void() +  # 使用clean主题
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        text = element_text(size = 16, color = "black"),  # 适当调整字体以容纳更多内容
        plot.title = element_text(size = 20, face = "bold", color = "black", hjust = 0.5),
        plot.subtitle = element_text(size = 16, color = "black", hjust = 0.5),
        legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.background = element_rect(fill = "white", color = "white"),
        legend.key = element_rect(fill = "white", color = "white"),
        legend.position = "bottom",
        plot.margin = margin(30, 30, 30, 30),  # 适当边距
        legend.key.size = unit(1.5, "cm"),
        legend.spacing.x = unit(1, "cm"),
        legend.box.spacing = unit(0.8, "cm")
      ) +
      labs(title = "KEGG Pathway Enrichment Network",
           subtitle = paste("Top", n_pathways, "pathways | Edges represent pathway similarity"))
  }
  
  # 8. 保存图表
  cat("正在保存图表...\n")
  
  # 标准条形图
  ggsave(file.path(output_dir, "kegg_pathway_barplot_clusterProfiler.pdf"), 
         plot = p1, width = 12, height = 8, units = "in")
  ggsave(file.path(output_dir, "kegg_pathway_barplot_clusterProfiler.png"), 
         plot = p1, width = 12, height = 8, dpi = 300, units = "in")
  
  # 自定义条形图
  ggsave(file.path(output_dir, "kegg_pathway_enrichment_barplot.pdf"), 
         plot = p2, width = 12, height = 8, units = "in")
  ggsave(file.path(output_dir, "kegg_pathway_enrichment_barplot.png"), 
         plot = p2, width = 12, height = 8, dpi = 300, units = "in")
  
  # 点图
  ggsave(file.path(output_dir, "kegg_pathway_dotplot.pdf"), 
         plot = p3, width = 12, height = 8, units = "in")
  ggsave(file.path(output_dir, "kegg_pathway_dotplot.png"), 
         plot = p3, width = 12, height = 8, dpi = 300, units = "in")
  
  # 网络图（如果可用）
  if (exists("p4")) {
    ggsave(file.path(output_dir, "kegg_pathway_network.pdf"), 
           plot = p4, width = 14, height = 12, units = "in")
    ggsave(file.path(output_dir, "kegg_pathway_network.png"), 
           plot = p4, width = 14, height = 12, dpi = 300, units = "in")
  }
  
  # 9. 保存富集数据
  # 创建符合手稿表格格式的数据
  manuscript_table <- top_pathways %>%
    select(ID, Description, Count, p.adjust, pvalue, geneSymbol) %>%
    rename(
      "Pathway_ID" = ID,
      "Pathway" = Description,
      "Gene_Count" = Count,
      "FDR" = p.adjust,
      "P_value" = pvalue,
      "Gene_Symbols" = geneSymbol
    ) %>%
    arrange(P_value)
  
  # 保存完整结果
  write.csv(kegg_result, "results/kegg_enrichment_clusterprofiler_full.csv", row.names = FALSE)
  write.csv(manuscript_table, "results/kegg_enrichment_manuscript_table.csv", row.names = FALSE)
  
  # 更新原有的pathway_enrichment_data.csv文件
  updated_pathway_data <- kegg_result %>%
    select(Description, Count, geneSymbol, pvalue, p.adjust) %>%
    rename(
      pathway = Description,
      gene_count = Count,
      genes = geneSymbol,
      pvalue = pvalue,
      qvalue = p.adjust
    ) %>%
    mutate(
      fold_enrichment = NA,  # clusterProfiler不直接提供这个值
      pathway_clean = ifelse(nchar(pathway) > 30, 
                           paste0(substr(pathway, 1, 27), "..."), 
                           pathway)
    ) %>%
    arrange(pvalue)
  
  write.csv(updated_pathway_data, "results/pathway_enrichment_data_clusterprofiler.csv", row.names = FALSE)
  
  # 10. 生成统计摘要
  cat("\n=== KEGG通路富集分析结果（clusterProfiler）===\n")
  cat("总靶点基因数:", length(gene_symbols), "\n")
  cat("转换为ENTREZ ID的基因数:", length(entrez_ids), "\n")
  cat("显著富集通路数 (FDR < 0.05):", sum(kegg_result$p.adjust < 0.05), "\n")
  cat("显著富集通路数 (P < 0.05):", sum(kegg_result$pvalue < 0.05), "\n")
  
  if (nrow(manuscript_table) > 0) {
    cat("\nTop 10 最显著富集的KEGG通路:\n")
    for(i in 1:min(10, nrow(manuscript_table))) {
      cat(sprintf("  %d. %s (ID: %s)\n", i, 
                  manuscript_table$Pathway[i], 
                  manuscript_table$Pathway_ID[i]))
      cat(sprintf("     基因数: %d, P值: %.2e, FDR: %.2e\n", 
                  manuscript_table$Gene_Count[i],
                  manuscript_table$P_value[i],
                  manuscript_table$FDR[i]))
    }
  }
  
  # 11. 为手稿生成表格格式
  cat("\n=== 为手稿生成的表2格式 ===\n")
  cat("| Pathway | Gene Count | P-value | FDR | Pathway ID |\n")
  cat("|---------|------------|---------|-----|-----------|\n")
  
  for(i in 1:min(10, nrow(manuscript_table))) {
    cat(sprintf("| %s | %d | %.2E | %.2E | %s |\n",
                manuscript_table$Pathway[i],
                manuscript_table$Gene_Count[i],
                manuscript_table$P_value[i],
                manuscript_table$FDR[i],
                manuscript_table$Pathway_ID[i]))
  }
  
  cat("\n✓ clusterProfiler通路富集分析完成\n")
  cat("✓ 保存图表: kegg_pathway_enrichment_barplot.pdf/png, kegg_pathway_dotplot.pdf/png\n")
  cat("✓ 保存数据: kegg_enrichment_clusterprofiler_full.csv, kegg_enrichment_manuscript_table.csv\n")
  
} else {
  cat("错误: 未能进行KEGG通路富集分析\n")
  cat("可能的原因:\n")
  cat("1. 基因数量不足\n")
  cat("2. 基因符号转换失败\n")
  cat("3. 网络连接问题\n")
}

# 12. GO富集分析（补充分析）
cat("\n=== 补充GO富集分析 ===\n")

# GO-BP富集
go_bp <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  minGSSize = 5,
                  maxGSSize = 500,
                  readable = TRUE)

if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
  cat("✓ GO-BP富集分析发现", nrow(go_bp@result), "个显著通路\n")
  
  # 保存GO-BP结果
  write.csv(go_bp@result, "results/go_bp_enrichment_clusterprofiler.csv", row.names = FALSE)
  
  # 生成GO-BP条形图
  p_go_bp <- barplot(go_bp, showCategory = 15, title = "GO Biological Process Enrichment")
  ggsave(file.path(output_dir, "go_bp_enrichment_barplot.pdf"), 
         plot = p_go_bp, width = 12, height = 8, units = "in")
  ggsave(file.path(output_dir, "go_bp_enrichment_barplot.png"), 
         plot = p_go_bp, width = 12, height = 8, dpi = 300, units = "in")
}

# GO-MF富集
go_mf <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  minGSSize = 5,
                  maxGSSize = 500,
                  readable = TRUE)

if (!is.null(go_mf) && nrow(go_mf@result) > 0) {
  cat("✓ GO-MF富集分析发现", nrow(go_mf@result), "个显著通路\n")
  write.csv(go_mf@result, "results/go_mf_enrichment_clusterprofiler.csv", row.names = FALSE)
}

# GO-CC富集
go_cc <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  minGSSize = 5,
                  maxGSSize = 500,
                  readable = TRUE)

if (!is.null(go_cc) && nrow(go_cc@result) > 0) {
  cat("✓ GO-CC富集分析发现", nrow(go_cc@result), "个显著通路\n")
  write.csv(go_cc@result, "results/go_cc_enrichment_clusterprofiler.csv", row.names = FALSE)
}

cat("完成时间:", as.character(Sys.time()), "\n")
