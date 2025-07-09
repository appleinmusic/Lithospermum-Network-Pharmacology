#!/usr/bin/env Rscript
# ==============================================================================
# Pathway Enrichment Analysis for High-Confidence Network (V3 - Strict Statistics)
# 职责：对高置信度网络 (score >= 700) 的靶点进行KEGG和GO富集分析。
# ==============================================================================

cat("=== High-Confidence Pathway Enrichment Analysis (V3) ===\n")
cat("开始时间:", as.character(Sys.time()), "\n")

# 加载必需包
suppressPackageStartupMessages({
  library(igraph) # 修正：添加igraph包
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# --- 1. 设置路径和加载高置信度靶点数据 ---
output_dir <- "results/figures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# **核心修改**: 读取高置信度网络对象
network_file <- "results/network/ppi_network_high_confidence.rds"
if (!file.exists(network_file)) {
  stop("错误：高置信度网络文件 '", network_file, "' 不存在。\n请先运行 03a_network_construction_high_confidence.R 脚本。")
}

g_core_network_high_confidence <- readRDS(network_file)

# 从网络对象中提取基因符号
gene_symbols <- V(g_core_network_high_confidence)$name_display

cat("✓ 成功载入", length(gene_symbols), "个高置信度核心网络中的靶点基因进行富集分析。\n")

# --- 2. 基因ID转换 ---
cat("--> 正在转换基因符号为ENTREZ ID...\n")
gene_entrez <- bitr(gene_symbols, 
                   fromType = "SYMBOL",
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

if (nrow(gene_entrez) == 0) {
  stop("无法将任何基因符号转换为ENTREZ ID，富集分析中止。")
}

cat("✓ 成功转换", nrow(gene_entrez), "个基因的ENTREZ ID。\n")
entrez_ids <- unique(gene_entrez$ENTREZID)

# --- 3. KEGG通路富集分析 (严格标准) ---
cat("--> 正在进行KEGG通路富集分析 (FDR < 0.05).\n")
FDR_CUTOFF <- 0.05
P_VALUE_CUTOFF <- 0.05

kegg_enrich <- enrichKEGG(gene = entrez_ids,
                         organism = 'hsa',
                         keyType = 'kegg',
                         pvalueCutoff = P_VALUE_CUTOFF,
                         pAdjustMethod = "BH",
                         minGSSize = 3, # 放宽minGSSize以适应更小的基因集
                         maxGSSize = 500,
                         qvalueCutoff = FDR_CUTOFF)

# --- 4. 结果处理和报告 ---
if (!is.null(kegg_enrich) && nrow(kegg_enrich@result) > 0) {
  cat("✓ 成功! 在FDR <", FDR_CUTOFF, "的阈值下，发现", nrow(kegg_enrich@result), "个显著富集的KEGG通路。\n")
  kegg_result <- as.data.frame(kegg_enrich)
  kegg_result$geneSymbol <- sapply(kegg_result$geneID, function(genes) {
    gene_ids <- unlist(strsplit(genes, "/"))
    symbols <- gene_entrez$SYMBOL[match(gene_ids, gene_entrez$ENTREZID)]
    paste(symbols[!is.na(symbols)], collapse = "/")
  })
  
  # 保存高置信度结果
  write.csv(kegg_result, "results/kegg_enrichment_clusterprofiler_full_high_confidence.csv", row.names = FALSE)
  cat("✓ 完整的高置信度KEGG富集结果已保存。\n")

  p_dot <- dotplot(kegg_enrich, showCategory = 15, title = "High-Confidence KEGG Pathway Enrichment (FDR < 0.05)", font.size = 12, color = "p.adjust")
  ggsave(file.path(output_dir, "kegg_pathway_dotplot_high_confidence.png"), plot = p_dot, width = 10, height = 8, dpi = 300, units = "in")
  cat("✓ 高置信度KEGG富集点图已保存。\n")

} else {
  cat("--- 在高置信度网络下，未发现显著富集的KEGG通路。---\n")
  file.create("results/kegg_enrichment_high_confidence_no_results.txt")
}

# --- 5. GO富集分析 (同样使用严格标准) ---
cat("--> 正在进行GO富集分析 (FDR < 0.05).\n")
go_enrich <- enrichGO(gene = entrez_ids,
                    OrgDb = org.Hs.eg.db,
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = P_VALUE_CUTOFF,
                    qvalueCutoff = FDR_CUTOFF,
                    readable = TRUE)

if (!is.null(go_enrich) && nrow(go_enrich@result) > 0) {
  cat("✓ 成功! 在FDR <", FDR_CUTOFF, "的阈值下，发现", nrow(go_enrich@result), "个显著富集的GO terms。\n")
  write.csv(as.data.frame(go_enrich), "results/go_enrichment_clusterprofiler_full_high_confidence.csv", row.names = FALSE)
  cat("✓ 完整的高置信度GO富集结果已保存。\n")

  go_enrich_simplified <- simplify(go_enrich, cutoff=0.7, by="p.adjust", select_fun=min)
  p_go_dot <- dotplot(go_enrich_simplified, showCategory = 20, title = "High-Confidence GO Enrichment (FDR < 0.05)", font.size = 10, color = "p.adjust")
  ggsave(file.path(output_dir, "Figure6_GO_Enrichment_Dotplot_high_confidence.png"), plot = p_go_dot, width = 12, height = 8, dpi = 300, units = "in")
  cat("✓ 高置信度GO富集点图已保存。\n")

} else {
  cat("--- 在高置信度网络下，未发现显著富集的GO terms。---\n")
  file.create("results/go_enrichment_high_confidence_no_results.txt")
}

cat("\n=== 高置信度富集分析完成 (V3) ===\n")
cat("完成时间:", as.character(Sys.time()), "\n")
