#!/usr/bin/env Rscript
# ==============================================================================
# Pathway Enrichment Analysis (V3 - Strict Statistics)
# 职责：严格按照预设的统计阈值(FDR < 0.05)进行KEGG和GO富集分析，
# 不进行任何自动化的阈值调整。
# ==============================================================================

cat("=== Pathway Enrichment Analysis (V3) ===
")
cat("开始时间:", as.character(Sys.time()), "
")

# 加载必需包
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
})

# --- 1. 设置路径和加载靶点数据 ---
output_dir <- "results/figures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# **核心修改**: 靶点列表应来源于最终的核心网络，而不是所有靶点
# 读取网络分析后得到的、在核心网络中的节点列表
network_nodes_file <- "results/tables/key_nodes_analysis.csv"
if (!file.exists(network_nodes_file)) {
  stop("错误：核心网络节点文件 '", network_nodes_file, "' 不存在。
请先运行 04_network_analysis_and_viz.R 脚本。")
}
targets_data <- readr::read_csv(network_nodes_file, show_col_types = FALSE)
# 使用 gene_symbol 列，因为它是最准确的
gene_symbols <- unique(targets_data$gene_symbol)
cat("✓ 成功载入", length(gene_symbols), "个核心网络中的靶点基因进行富集分析。
")

# --- 2. 基因ID转换 ---
cat("--> 正在转换基因符号为ENTREZ ID...
")
gene_entrez <- bitr(gene_symbols, 
                   fromType = "SYMBOL",
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

if (nrow(gene_entrez) == 0) {
  stop("无法将任何基因符号转换为ENTREZ ID，富集分析中止。")
}

cat("✓ 成功转换", nrow(gene_entrez), "个基因的ENTREZ ID。
")
entrez_ids <- unique(gene_entrez$ENTREZID)

# --- 3. KEGG通路富集分析 (严格标准) ---
cat("--> 正在进行KEGG通路富集分析 (FDR < 0.05)...
")

# **核心修改**: 移除自动放宽阈值的逻辑，只执行一次严格分析
FDR_CUTOFF <- 0.05
P_VALUE_CUTOFF <- 0.05

kegg_enrich <- enrichKEGG(gene = entrez_ids,
                         organism = 'hsa',
                         keyType = 'kegg',
                         pvalueCutoff = P_VALUE_CUTOFF,
                         pAdjustMethod = "BH",
                         minGSSize = 5,
                         maxGSSize = 500,
                         qvalueCutoff = FDR_CUTOFF)

# --- 4. 结果处理和报告 ---
if (!is.null(kegg_enrich) && nrow(kegg_enrich@result) > 0) {
  cat("✓ 成功! 在FDR <", FDR_CUTOFF, "的阈值下，发现", nrow(kegg_enrich@result), "个显著富集的KEGG通路。
")
  
  kegg_result <- as.data.frame(kegg_enrich)
  
  # ... (后续的可视化和保存代码与您原脚本类似，此处为保证逻辑正确性而重写关键部分)
  
  # 添加基因符号以便阅读
  kegg_result$geneSymbol <- sapply(kegg_result$geneID, function(genes) {
    gene_ids <- unlist(strsplit(genes, "/"))
    symbols <- gene_entrez$SYMBOL[match(gene_ids, gene_entrez$ENTREZID)]
    paste(symbols[!is.na(symbols)], collapse = "/")
  })
  
  # 保存完整结果
  write.csv(kegg_result, "results/kegg_enrichment_clusterprofiler_full.csv", row.names = FALSE)
  cat("✓ 完整的KEGG富集结果已保存。
")

  # 生成并保存关键图表 (以点图为例)
  p_dot <- dotplot(kegg_enrich, 
                   showCategory = 15, # 显示前15个
                   title = "KEGG Pathway Enrichment (FDR < 0.05)",
                   font.size = 12,
                   color = "p.adjust")
  
  ggsave(file.path(output_dir, "kegg_pathway_dotplot.png"), 
         plot = p_dot, width = 10, height = 8, dpi = 300, units = "in")
  cat("✓ KEGG富集点图已保存。
")

} else {
  cat("--- 未发现显著富集的KEGG通路 ---
")
  cat("在严格的FDR <", FDR_CUTOFF, "阈值下，没有KEGG通路达到统计学显著性。
")
  cat("这是一个有效的阴性结果，将在结果文件中记录。
")
  
  # 创建一个空的占位符文件，表明分析已运行但无结果
  file.create("results/kegg_enrichment_no_significant_results.txt")
}

# --- 5. GO富集分析 (同样使用严格标准) ---
cat("--> 正在进行GO富集分析 (FDR < 0.05)...
")

go_enrich <- enrichGO(gene = entrez_ids,
                    OrgDb = org.Hs.eg.db,
                    ont = "ALL", # 一次性分析BP, CC, MF
                    pAdjustMethod = "BH",
                    pvalueCutoff = P_VALUE_CUTOFF,
                    qvalueCutoff = FDR_CUTOFF,
                    readable = TRUE) # 将ENTREZ ID直接转换为基因符号

if (!is.null(go_enrich) && nrow(go_enrich@result) > 0) {
  cat("✓ 成功! 在FDR <", FDR_CUTOFF, "的阈值下，发现", nrow(go_enrich@result), "个显著富集的GO terms。
")
  
  # 保存GO结果
  write.csv(as.data.frame(go_enrich), "results/go_enrichment_clusterprofiler_full.csv", row.names = FALSE)
  cat("✓ 完整的GO富集结果已保存。
")

  # --- 改进GO点图：先去冗余，再绘图 ---
  cat("--> 正在使用simplify()函数去除GO term的冗余...\n")
  go_enrich_simplified <- simplify(go_enrich, cutoff=0.7, by="p.adjust", select_fun=min)
  
  cat("✓ GO term去冗余完成。原始terms:", nrow(go_enrich@result), "，简化后terms:", nrow(go_enrich_simplified@result), "\n")

  # --- 生成并保存GO富集网络图 (稳健方案) ---
  # 检查简化后的结果是否有足够的terms用于网络分析
  if (nrow(go_enrich_simplified@result) >= 5) {
    cat("--> 正在生成GO富集网络图...\n")
    
    # 使用原始的go_enrich对象进行网络图绘制，避免setReadable冲突
    tryCatch({
      p_go_emap <- emapplot(go_enrich_simplified, showCategory = min(30, nrow(go_enrich_simplified@result)), layout = "fr") + 
        labs(title = "GO Enrichment Map (FDR < 0.05)",
             subtitle = "Similar GO terms are clustered together") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5))
              
      ggsave(file.path(output_dir, "Figure6_GO_Enrichment_Map.png"), 
             plot = p_go_emap, width = 12, height = 10, dpi = 300, units = "in")
      
      cat("✓ 已保存GO富集网络图 (emapplot)。\n")
    }, error = function(e) {
      cat("! emapplot生成失败，改用dotplot作为替代可视化方案。\n")
      cat("  错误信息:", e$message, "\n")
      
      # 使用稳健的dotplot替代方案
      p_go_dot <- dotplot(go_enrich_simplified, 
                         showCategory = min(20, nrow(go_enrich_simplified@result)),
                         title = "GO Enrichment Analysis (FDR < 0.05)",
                         font.size = 10,
                         color = "p.adjust") +
        theme(axis.text.y = element_text(size = 8))
      
      ggsave(file.path(output_dir, "Figure6_GO_Enrichment_Dotplot.png"), 
             plot = p_go_dot, width = 12, height = 8, dpi = 300, units = "in")
      
      cat("✓ 已保存GO富集点图作为替代方案。\n")
    })
  } else {
    cat("! GO terms数量过少(", nrow(go_enrich_simplified@result), ")，无法生成网络图。\n")
    
    # 生成简单的条形图
    p_go_bar <- barplot(go_enrich_simplified, 
                       showCategory = nrow(go_enrich_simplified@result),
                       title = "GO Enrichment Analysis (FDR < 0.05)")
    
    ggsave(file.path(output_dir, "Figure6_GO_Enrichment_Barplot.png"), 
           plot = p_go_bar, width = 10, height = 6, dpi = 300, units = "in")
    
    cat("✓ 已保存GO富集条形图。\n")
  }

} else {
  cat("--- 未发现显著富集的GO terms ---
")
  cat("在严格的FDR <", FDR_CUTOFF, "阈值下，没有GO term达到统计学显著性。
")
  file.create("results/go_enrichment_no_significant_results.txt")
}

cat("
=== 富集分析完成 (V3) ===
")
cat("✓ 所有分析均严格遵守FDR <", FDR_CUTOFF, "的统计标准。
")
cat("✓ 结果和图表已更新。
")
cat("完成时间:", as.character(Sys.time()), "
")