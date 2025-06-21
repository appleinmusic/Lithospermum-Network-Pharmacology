#!/usr/bin/env Rscript
# ADMET Properties Visualization for Lithospermum erythrorhizon

# 设置更长的超时时间以处理大文件下载
options(timeout = 300)

# 加载必需包
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
  library(gridExtra)
  library(grid)
  library(viridis)
  # 新增 rcdk 和 rJava 用于分子描述符计算
  # 请确保已安装Java环境
  if (!requireNamespace("rcdk", quietly = TRUE)) install.packages("rcdk", repos = "https://cloud.r-project.org")
  if (!requireNamespace("rJava", quietly = TRUE)) install.packages("rJava", repos = "https://cloud.r-project.org")
  library(rJava)
  library(rcdk)
})

# --- 专家级解决方案：正确的rcdk初始化 ---
cat("--> [专家解决方案] 使用rcdk推荐的初始化方法...\n")

# 方法1：使用rcdk内置的初始化
tryCatch({
  # 卸载并重新加载rJava以确保干净状态
  if("rJava" %in% loadedNamespaces()) {
    unloadNamespace("rJava")
  }
  if("rcdk" %in% loadedNamespaces()) {
    unloadNamespace("rcdk")
  }
  
  # 重新加载包
  library(rJava, quietly = TRUE)
  library(rcdk, quietly = TRUE)
  
  cat("✓ 方法1：标准rcdk初始化成功\n")
}, error = function(e) {
  cat("✗ 方法1失败:", e$message, "\n")
  
  # 方法2：手动但更智能的初始化
  tryCatch({
    cat("--> 尝试方法2：手动智能初始化...\n")
    
    # 强制初始化Java
    .jinit(force.init = TRUE)
    
    # 获取正确的jar路径
    rcdk_path <- system.file(package = "rcdk")
    rcdklibs_path <- system.file(package = "rcdklibs")
    
    cat("rcdk路径:", rcdk_path, "\n")
    cat("rcdklibs路径:", rcdklibs_path, "\n")
    
    # 添加rcdklibs的jar文件
    rcdklibs_jars <- list.files(file.path(rcdklibs_path, "cont"), 
                               pattern = "\\.jar$", full.names = TRUE)
    if (length(rcdklibs_jars) > 0) {
      for (jar in rcdklibs_jars) {
        .jaddClassPath(jar)
      }
      cat("✓ 添加了", length(rcdklibs_jars), "个rcdklibs jar文件\n")
    }
    
    # 添加rcdk的jar文件
    rcdk_jars <- list.files(file.path(rcdk_path, "cont"), 
                           pattern = "\\.jar$", full.names = TRUE)
    if (length(rcdk_jars) > 0) {
      for (jar in rcdk_jars) {
        .jaddClassPath(jar)
      }
      cat("✓ 添加了", length(rcdk_jars), "个rcdk jar文件\n")
    }
    
    cat("✓ 方法2：手动智能初始化完成\n")
    
  }, error = function(e2) {
    cat("✗ 方法2也失败:", e2$message, "\n")
    cat("将切换到备用方案...\n")
  })
})

# 显示当前classpath（用于调试）
cat("当前Java Classpath包含", length(.jclassPath()), "个条目\n")
# ---

cat("=== ADMET性质可视化分析 (V2 - 严谨预测版) ===\n")
cat("开始时间:", as.character(Sys.time()), "\n")

# 设置目录路径
if (!file.exists("data") || !file.exists("results")) {
  stop("请确保在项目根目录运行此脚本（包含data/和results/文件夹的目录）")
}
output_dir <- "results/figures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 检查并读取数据文件
ingredients_file <- "data/processed/lithospermum_ingredients_filtered.tsv"
if (!file.exists(ingredients_file)) {
  stop("成分数据文件不存在，请先运行01_complete_data_loading.R脚本")
}

# 读取过滤后的成分数据
ingredients_data <- readr::read_tsv(ingredients_file)


cat("✓ 载入成分数据:", nrow(ingredients_data), "条记录\n")

# --- 方法学重大更新：使用RCDK进行真实的ADMET性质预测 ---
cat("--> 开始进行基于RCDK的分子描述符计算...\n")

# 确保SMILES列存在
if (!"SMILES_ingredient" %in% colnames(ingredients_data)) {
  stop("错误：输入文件中缺少 'SMILES_ingredient' 列，无法进行描述符计算。")
}

# 解析SMILES为rcdk分子对象
# 添加错误处理，跳过无法解析的SMILES
mols <- rcdk::parse.smiles(ingredients_data$SMILES_ingredient)
valid_indices <- !sapply(mols, is.null)
if(sum(!valid_indices) > 0) {
    cat("警告：有", sum(!valid_indices), "个SMILES字符串无法被解析，将被跳过。\n")
    mols <- mols[valid_indices]
    ingredients_data_valid <- ingredients_data[valid_indices, ]
} else {
    ingredients_data_valid <- ingredients_data
}


# 计算描述符 - 使用正确的描述符类名
# 注意：部分分子的某些描述符可能计算失败，返回NA
descriptors <- rcdk::eval.desc(mols, 
                                which.desc = c("org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor",
                                              "org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor",
                                              "org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor",
                                              "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor",
                                              "org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor"),
                                verbose = FALSE)

cat("✓ RCDK分子描述符计算完成\n")

# 将计算结果与原始数据合并
# 检查实际的列名并进行正确的重命名
cat("描述符计算结果列名:", colnames(descriptors), "\n")
predicted_admet <- as.data.frame(descriptors) %>% 
  rename(LogP_pred = ALogP, nHA_pred = nHBAcc, nHD_pred = nHBDon, TPSA_pred = TopoPSA, MW_pred = MW)

active_ingredients_pred <- bind_cols(ingredients_data_valid, predicted_admet)

# --- 基于真实预测值进行ADMET筛选 ---
# 移除旧的、有缺陷的OB和DL筛选，采用更可靠的Lipinski规则和TPSA
active_ingredients <- active_ingredients_pred %>%
  mutate(
    # 处理计算失败的NA值，将其视为不通过筛选
    MW_pred = ifelse(is.na(MW_pred), 9999, MW_pred),
    LogP_pred = ifelse(is.na(LogP_pred), 99, LogP_pred),
    nHD_pred = ifelse(is.na(nHD_pred), 99, nHD_pred),
    nHA_pred = ifelse(is.na(nHA_pred), 99, nHA_pred),
    TPSA_pred = ifelse(is.na(TPSA_pred), 999, TPSA_pred),

    # 基于RCDK的真实预测值进行筛选
    passes_lipinski_mw = MW_pred <= 500,
    passes_lipinski_logp = LogP_pred <= 5,
    passes_lipinski_hbd = nHD_pred <= 5,
    passes_lipinski_hba = nHA_pred <= 10,
    passes_tpsa = TPSA_pred <= 140,
    
    # 总体ADMET合规性 (移除OB和DL)
    passes_ADMET = passes_lipinski_mw & passes_lipinski_logp & 
                   passes_lipinski_hbd & passes_lipinski_hba & 
                   passes_tpsa,
    
    # 化合物名称简化
    compound = pref_name
  )

# 筛选通过所有标准的核心活性成分
final_active_compounds <- active_ingredients %>% 
  filter(passes_ADMET == TRUE)

cat("--> ADMET筛选完成。有效化合物:", nrow(active_ingredients), ", 通过筛选:", nrow(final_active_compounds), "\n")

# 保存包含所有预测值和筛选结果的完整数据
write.csv(active_ingredients, "results/admet_analysis_full_data_predicted.csv", row.names = FALSE)

# 仅保存通过筛选的化合物，用于下游分析
# 这是替换旧的 admet_analysis_data.csv 的关键文件
write.csv(final_active_compounds, "results/admet_analysis_data.csv", row.names = FALSE)

# --- 可视化部分更新 ---
# 后续的ggplot代码将使用更新后的 `active_ingredients` 数据框和新列名

# 1. ADMET筛选标准可视化 (更新标准)
admet_standards <- data.frame(
  Property = c("MW (Da)", "LogP", "nHA", "nHD", "TPSA (Å²)"),
  Threshold = c(500, 5, 10, 5, 140),
  Type = c("<=", "<=", "<=", "<=", "<=")
)

# 2. 四合一ADMET图的各个子图

# Panel A: MW vs LogP散点图 
p1 <- ggplot(active_ingredients, aes(x = MW_pred, y = LogP_pred)) +
  geom_point(aes(color = passes_ADMET), size = 3, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "#E74C3C", "TRUE" = "#2ECC71"), 
                     labels = c("FALSE" = "Failed", "TRUE" = "Passed")) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "#E74C3C", linewidth = 0.8) +
  geom_vline(xintercept = 500, linetype = "dashed", color = "#E74C3C", linewidth = 0.8) +
  labs(
    title = "Panel A: MW vs LogP (Predicted)",
    x = "Molecular Weight (Da)",
    y = "ALOGP",
    color = "ADMET Status"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    legend.position = "none"
  )

# Panel B: TPSA分布直方图
p2 <- ggplot(active_ingredients, aes(x = TPSA_pred)) +
  geom_histogram(bins = 8, fill = "#3498DB", alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 140, linetype = "dashed", color = "#E74C3C", linewidth = 0.8) +
  labs(
    title = "Panel B: TPSA Distribution (Predicted)",
    x = "TPSA (Å²)",
    y = "Count"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
  annotate("text", x = 145, y = max(table(cut(active_ingredients$TPSA_pred, 8))) * 0.8, 
           label = "TPSA ≤ 140", color = "#E74C3C", size = 3, angle = 90, fontface = "bold")

# Panel C: Lipinski规则符合率
lipinski_compliance <- data.frame(
  Rule = c("MW ≤ 500", "LogP ≤ 5", "HBD ≤ 5", "HBA ≤ 10", "All Rules Pass"),
  Pass_Rate = c(
    sum(active_ingredients$passes_lipinski_mw, na.rm = TRUE) / nrow(active_ingredients) * 100,
    sum(active_ingredients$passes_lipinski_logp, na.rm = TRUE) / nrow(active_ingredients) * 100,
    sum(active_ingredients$passes_lipinski_hbd, na.rm = TRUE) / nrow(active_ingredients) * 100,
    sum(active_ingredients$passes_lipinski_hba, na.rm = TRUE) / nrow(active_ingredients) * 100,
    sum(active_ingredients$passes_ADMET, na.rm = TRUE) / nrow(active_ingredients) * 100
  )
)

p3 <- ggplot(lipinski_compliance, aes(x = reorder(Rule, Pass_Rate), y = Pass_Rate)) +
  geom_col(fill = "#E67E22", alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_text(aes(label = paste0(round(Pass_Rate, 1), "%")), 
            hjust = -0.1, size = 3, fontface = "bold") +
  coord_flip() +
  labs(
    title = "Panel C: Lipinski Rule Compliance (Predicted)",
    x = "Lipinski Rules",
    y = "Compliance Rate (%)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major.x = element_line(color = "#E5E5E5", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
  scale_y_continuous(limits = c(0, 105), expand = c(0, 0))

# Panel D: 替换为更有意义的图，例如 nHA vs nHD
p4 <- ggplot(active_ingredients, aes(x = nHD_pred, y = nHA_pred)) +
  geom_point(aes(color = passes_ADMET), size = 3, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "#E74C3C", "TRUE" = "#2ECC71"), 
                     labels = c("FALSE" = "Failed", "TRUE" = "Passed")) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "#E74C3C", linewidth = 0.8) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "#E74C3C", linewidth = 0.8) +
  labs(
    title = "Panel D: H-Bond Donors vs Acceptors",
    x = "Hydrogen Bond Donors (nHD)",
    y = "Hydrogen Bond Acceptors (nHA)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    legend.position = "bottom"
  )

# 组合四个子图为一个四合一图
combined_plot <- grid.arrange(p1, p2, p3, p4, 
                              ncol = 2, nrow = 2,
                              top = textGrob("ADMET Properties Analysis of L. erythrorhizon Compounds (Predicted)", 
                                           gp = gpar(fontsize = 16, fontface = "bold")))

# 保存其他独立分析图（保留原有的分析图）
if (nrow(final_active_compounds) > 0) {
    passed_compounds <- final_active_compounds %>%
      arrange(desc(MW_pred)) %>%
      slice_head(n = 15) %>%
      mutate(compound = ifelse(nchar(compound) > 20, 
                              paste0(substr(compound, 1, 17), "..."), 
                              compound))

    p_mw_dist <- ggplot(passed_compounds, aes(x = reorder(compound, MW_pred), y = MW_pred)) +
      geom_col(fill = "#3498DB", alpha = 0.8, color = "black", linewidth = 0.2) +
      coord_flip() +
      labs(
        title = "Molecular Weight Distribution",
        subtitle = "Top 15 ADMET-Compliant Compounds by Molecular Weight",
        x = "Active Compounds",
        y = "Molecular Weight (Da)"
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        panel.background = element_rect(fill = "white", color = "black"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major.x = element_line(color = "#E5E5E5", linewidth = 0.3),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text.y = element_text(size = 9)
      )

    # 保存独立的分子量分布图
    ggsave(file.path(output_dir, "Figure6_Lipinski_Analysis.png"), 
           plot = p_mw_dist, width = 8, height = 6, dpi = 600, units = "in")
    ggsave(file.path(output_dir, "Figure6_Lipinski_Analysis.pdf"), 
           plot = p_mw_dist, width = 8, height = 6, units = "in")
}


# 按TPSA的top化合物图
if (nrow(final_active_compounds) > 0) {
  top_compounds <- final_active_compounds %>%
    arrange(desc(TPSA_pred)) %>%
    slice_head(n = 10) %>%
    mutate(compound = ifelse(nchar(compound) > 25, 
                            paste0(substr(compound, 1, 22), "..."), 
                            compound))
  
  p_tpsa <- ggplot(top_compounds, aes(x = reorder(compound, TPSA_pred), y = TPSA_pred)) +
    geom_col(fill = "#E67E22", alpha = 0.8, color = "black", linewidth = 0.2) +
    geom_text(aes(label = round(TPSA_pred, 1)), hjust = -0.1, size = 3, fontface = "bold") +
    coord_flip() +
    labs(
      title = "Top 10 Compounds by Polar Surface Area",
      subtitle = "TPSA Values for ADMET-Compliant Compounds",
      x = "Compound Name",
      y = "Topological Polar Surface Area (A^2)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      panel.background = element_rect(fill = "white", color = "black"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major.x = element_line(color = "#E5E5E5", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.text.y = element_text(size = 9)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
  ggsave(file.path(output_dir, "Figure7_Top_Compounds.png"), 
         plot = p_tpsa, width = 8, height = 6, dpi = 600, units = "in")
  ggsave(file.path(output_dir, "Figure7_Top_Compounds.pdf"), 
         plot = p_tpsa, width = 8, height = 6, units = "in")
}

# 保存四合一图
ggsave(file.path(output_dir, "Figure5_ADMET_Properties.png"), 
       plot = combined_plot, width = 12, height = 10, dpi = 600, units = "in")
ggsave(file.path(output_dir, "Figure5_ADMET_Properties.pdf"), 
       plot = combined_plot, width = 12, height = 10, units = "in")

# 保存单独的子图（备用）
ggsave(file.path(output_dir, "Figure5A_MW_LogP.png"), 
       plot = p1, width = 6, height = 5, dpi = 600, units = "in")
ggsave(file.path(output_dir, "Figure5B_TPSA_Distribution.png"), 
       plot = p2, width = 6, height = 5, dpi = 600, units = "in")
ggsave(file.path(output_dir, "Figure5C_Lipinski_Compliance.png"), 
       plot = p3, width = 6, height = 5, dpi = 600, units = "in")
ggsave(file.path(output_dir, "Figure5D_Hbond_Plot.png"), 
       plot = p4, width = 6, height = 5, dpi = 600, units = "in")

# 6. 生成ADMET统计摘要
admet_summary <- list(
  total_compounds_screened = nrow(active_ingredients),
  passed_admet = sum(active_ingredients$passes_ADMET, na.rm = TRUE),
  admet_success_rate = paste0(round(sum(active_ingredients$passes_ADMET, na.rm = TRUE) / nrow(active_ingredients) * 100, 1), "%"),
  mw_mean = round(mean(active_ingredients$MW_pred, na.rm = TRUE), 2),
  mw_range = paste(round(min(active_ingredients$MW_pred, na.rm = TRUE), 2), 
                   round(max(active_ingredients$MW_pred, na.rm = TRUE), 2), sep = "-"),
  logp_mean = round(mean(active_ingredients$LogP_pred, na.rm = TRUE), 3),
  logp_range = paste(round(min(active_ingredients$LogP_pred, na.rm = TRUE), 3), 
                     round(max(active_ingredients$LogP_pred, na.rm = TRUE), 3), sep = "-"),
  tpsa_mean = round(mean(active_ingredients$TPSA_pred, na.rm = TRUE), 2),
  tpsa_range = paste(round(min(active_ingredients$TPSA_pred, na.rm = TRUE), 2), 
                     round(max(active_ingredients$TPSA_pred, na.rm = TRUE), 2), sep = "-")
)

# 保存统计摘要
writeLines(
  c("=== ADMET Properties Summary (Predicted) ===",
    paste("Total Compounds Screened:", admet_summary$total_compounds_screened),
    paste("ADMET Compliant Compounds:", admet_summary$passed_admet),
    paste("ADMET Success Rate:", admet_summary$admet_success_rate),
    paste("MW Mean (Da):", admet_summary$mw_mean),
    paste("MW Range (Da):", admet_summary$mw_range),
    paste("LogP Mean:", admet_summary$logp_mean),
    paste("LogP Range:", admet_summary$logp_range),
    paste("TPSA Mean (A^2):", admet_summary$tpsa_mean),
    paste("TPSA Range (A^2):", admet_summary$tpsa_range),
    "",
    "=== Four-Panel Figure Components ===",
    "Panel A: MW vs LogP scatter plot (Predicted)",
    "Panel B: TPSA distribution histogram (Predicted)", 
    "Panel C: Lipinski rules compliance rates (Predicted)",
    "Panel D: H-Bond Donors vs Acceptors plot"),
  file.path(output_dir, "ADMET_analysis_summary.txt")
)

cat("=== ADMET性质分析完成 (V2) ===\n")
cat("✓ 生成主要四合一图: Figure5_ADMET_Properties.png\n")
cat("✓ 保存路径:", output_dir, "\n")
cat("完成时间:", as.character(Sys.time()), "\n")