#!/usr/bin/env Rscript
# ==============================================================================
# ADMET-Based Compound Screening and Visualization (V3 - Unified Workflow)
# 职责：对所有初始化合物执行一次性的、基于rcdk的de novo预测和严格的
# Lipinski规则筛选，并生成可视化图表。
# ==============================================================================

# 设置和加载包
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
  library(gridExtra)
  library(grid)
  library(viridis)
  # 确保rcdk和rJava已安装
  if (!requireNamespace("rcdk", quietly = TRUE)) install.packages("rcdk", repos = "https://cloud.r-project.org")
  if (!requireNamespace("rJava", quietly = TRUE)) install.packages("rJava", repos = "https://cloud.r-project.org")
  library(rJava)
  library(rcdk)
})

# 初始化Java和rcdk
# ... (保留了您脚本中稳健的rcdk初始化代码)

cat("=== ADMET Screening and Visualization (V3) ===
")
cat("开始时间:", as.character(Sys.time()), "
")

# --- 1. 设置路径和加载数据 ---
output_dir <- "results/figures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# **核心修改**: 读取上一步生成的、未经任何筛选的初始化合物文件
ingredients_file <- "data/processed/lithospermum_ingredients_curated.tsv"
if (!file.exists(ingredients_file)) {
  stop("错误：初始化合物文件 '", ingredients_file, "' 不存在。
请先运行 01_data_preparation.R 脚本。")
}

ingredients_initial <- readr::read_tsv(ingredients_file, show_col_types = FALSE)
cat("✓ 成功载入", nrow(ingredients_initial), "个初始化合物进行筛选。
")

# --- 2. De Novo ADMET性质预测 (基于RCDK) ---
cat("--> 正在对所有初始化合物进行de novo ADMET性质预测...
")

# 确保SMILES列存在 (列名在CMAUP数据中为SMILES)
if (!"SMILES" %in% colnames(ingredients_initial)) {
  stop("错误：输入文件中缺少 'SMILES' 列，无法进行描述符计算。")
}

# 解析SMILES为rcdk分子对象，处理可能存在的错误
mols <- rcdk::parse.smiles(ingredients_initial$SMILES)
valid_indices <- !sapply(mols, is.null)
if(sum(!valid_indices) > 0) {
    cat("警告：有", sum(!valid_indices), "个SMILES字符串无法被解析，将被跳过。
")
    mols <- mols[valid_indices]
    ingredients_valid_smiles <- ingredients_initial[valid_indices, ]
} else {
    ingredients_valid_smiles <- ingredients_initial
}

# 计算描述符
descriptors <- rcdk::eval.desc(mols, 
                                which.desc = c("org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor",
                                              "org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor",
                                              "org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor",
                                              "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor",
                                              "org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor"),
                                verbose = FALSE)

cat("✓ RCDK分子描述符计算完成。
")

# 合并预测结果
predicted_admet <- as.data.frame(descriptors) %>% 
  rename(LogP_pred = ALogP, nHA_pred = nHBAcc, nHD_pred = nHBDon, TPSA_pred = TopoPSA, MW_pred = MW)

all_compounds_with_pred <- bind_cols(ingredients_valid_smiles, predicted_admet)

# --- 3. 基于预测值进行唯一一次的严格筛选 ---
cat("--> 正在应用Lipinski's Rule of Five进行筛选...
")

# 定义筛选标准
LIPINSKI_MW_MAX <- 500
LIPINSKI_LOGP_MAX <- 5
LIPINSKI_HBD_MAX <- 5
LIPINSKI_HBA_MAX <- 10
TPSA_MAX <- 140

# 应用筛选
compounds_screened <- all_compounds_with_pred %>%
  mutate(
    # 处理计算失败的NA值，将其视为不通过筛选
    MW_pred = ifelse(is.na(MW_pred), 9999, MW_pred),
    LogP_pred = ifelse(is.na(LogP_pred), 99, LogP_pred),
    nHD_pred = ifelse(is.na(nHD_pred), 99, nHD_pred),
    nHA_pred = ifelse(is.na(nHA_pred), 99, nHA_pred),
    TPSA_pred = ifelse(is.na(TPSA_pred), 999, TPSA_pred),

    # 判断是否通过各项标准
    passes_lipinski_mw = MW_pred <= LIPINSKI_MW_MAX,
    passes_lipinski_logp = LogP_pred <= LIPINSKI_LOGP_MAX,
    passes_lipinski_hbd = nHD_pred <= LIPINSKI_HBD_MAX,
    passes_lipinski_hba = nHA_pred <= LIPINSKI_HBA_MAX,
    passes_tpsa = TPSA_pred <= TPSA_MAX,
    
    # 最终是否通过所有筛选
    passes_ADMET = passes_lipinski_mw & passes_lipinski_logp & 
                   passes_lipinski_hbd & passes_lipinski_hba & 
                   passes_tpsa,
    
    compound = pref_name
  )

# 筛选出最终的核心活性成分
final_core_compounds <- compounds_screened %>% 
  filter(passes_ADMET == TRUE)

cat("--> ADMET筛选完成。
")
cat("    初始化合物总数:", nrow(ingredients_initial), "
")
cat("    通过筛选的核心化合物数:", nrow(final_core_compounds), "
")

# --- 4. 保存筛选结果 ---
# 保存包含所有预测值和筛选过程的完整数据，用于透明度
write.csv(compounds_screened, "results/admet_analysis_full_data_predicted.csv", row.names = FALSE)

# **核心输出**：仅保存通过筛选的化合物，作为所有下游分析的唯一输入源
write.csv(final_core_compounds, "results/admet_analysis_data.csv", row.names = FALSE)
cat("✓ 核心化合物数据已保存至: results/admet_analysis_data.csv
")

# **关键新增**：根据筛选后的化合物更新靶点数据
cat("--> 正在根据筛选后的化合物更新靶点数据...
")

# 加载初始靶点数据
targets_initial_file <- "data/processed/lithospermum_targets_initial.tsv"
if (!file.exists(targets_initial_file)) {
  stop("错误：初始靶点文件 '", targets_initial_file, "' 不存在。
请先运行 01_data_preparation.R 脚本。")
}

targets_initial <- readr::read_tsv(targets_initial_file, show_col_types = FALSE)

# 获取通过ADMET筛选的化合物ID列表
filtered_ingredient_ids <- final_core_compounds$np_id

# 筛选对应的靶点数据
targets_filtered <- targets_initial %>%
  filter(Ingredient_ID %in% filtered_ingredient_ids)

# 保存筛选后的靶点数据
write_tsv(targets_filtered, "data/processed/lithospermum_targets.tsv")
cat("✓ 筛选后的靶点数据已保存至: data/processed/lithospermum_targets.tsv
")
cat("    筛选后靶点关联数:", nrow(targets_filtered), "
")
cat("    筛选后独特靶点数:", length(unique(targets_filtered$Gene_Symbol)), "
")

# 同时保存筛选后的化合物数据为TSV格式（与其他处理文件格式一致）
write_tsv(final_core_compounds, "data/processed/lithospermum_ingredients_filtered.tsv")
cat("✓ 筛选后的化合物数据也已保存为TSV格式: data/processed/lithospermum_ingredients_filtered.tsv
")

# --- 5. 可视化分析 ---
# (这部分的可视化代码与您原脚本一致，因为它已经是基于预测值进行的，
#  现在它将作用于未经预筛选的、更完整的数据集上，使其分析更具代表性)

# Panel A: MW vs LogP散点图
p1 <- ggplot(compounds_screened, aes(x = MW_pred, y = LogP_pred)) +
  geom_point(aes(color = passes_ADMET), size = 3, alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "#E74C3C", "TRUE" = "#2ECC71"), name = "Passes Filter") +
  geom_hline(yintercept = LIPINSKI_LOGP_MAX, linetype = "dashed", color = "#E74C3C") +
  geom_vline(xintercept = LIPINSKI_MW_MAX, linetype = "dashed", color = "#E74C3C") +
  labs(title = "A: Molecular Weight vs. LogP", x = "Molecular Weight (Da)", y = "ALOGP") +
  theme_bw() + theme(legend.position = "none")

# Panel B: TPSA Distribution
p2 <- ggplot(compounds_screened, aes(x = TPSA_pred)) +
  geom_histogram(aes(fill = passes_ADMET), bins = 20, alpha = 0.8, position = "identity") +
  scale_fill_manual(values = c("FALSE" = "#E74C3C", "TRUE" = "#2ECC71"), name = "Passes Filter") +
  geom_vline(xintercept = TPSA_MAX, linetype = "dashed", color = "#E74C3C") +
  labs(title = "B: Topological Polar Surface Area (TPSA)", x = "TPSA (Å²)", y = "Count") +
  theme_bw() + theme(legend.position = "none")

# Panel C: Compliance with Lipinski's Rule of Five
lipinski_summary <- compounds_screened %>%
  select(passes_lipinski_mw, passes_lipinski_logp, passes_lipinski_hbd, passes_lipinski_hba) %>%
  summarise_all(sum) %>%
  tidyr::gather(key = "Rule", value = "Count") %>%
  mutate(Rule = factor(Rule, levels = c("passes_lipinski_mw", "passes_lipinski_logp", "passes_lipinski_hbd", "passes_lipinski_hba"),
                       labels = c("MW <= 500", "LogP <= 5", "HBD <= 5", "HBA <= 10")))

p3 <- ggplot(lipinski_summary, aes(x = Rule, y = Count)) +
  geom_col(aes(fill = Rule), alpha = 0.9) +
  geom_text(aes(label = Count), vjust = -0.5, fontface = "bold") +
  scale_fill_viridis_d(option = "C") +
  labs(title = "C: Compliance with Lipinski's Rules", x = "Lipinski's Rule of Five", y = "Number of Passing Compounds") +
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))

# Panel D: Hydrogen Bond Donors vs. Acceptors
p4 <- ggplot(compounds_screened, aes(x = nHD_pred, y = nHA_pred)) +
  geom_point(aes(color = passes_ADMET), size = 3, alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "#E74C3C", "TRUE" = "#2ECC71"), name = "Passes Filter") +
  geom_hline(yintercept = LIPINSKI_HBA_MAX, linetype = "dashed", color = "#E74C3C") +
  geom_vline(xintercept = LIPINSKI_HBD_MAX, linetype = "dashed", color = "#E74C3C") +
  labs(title = "D: H-Bond Donors vs. Acceptors", x = "Hydrogen Bond Donors (HBD)", y = "Hydrogen Bond Acceptors (HBA)") +
  theme_bw() + theme(legend.position = "right")

# Combine all plots into a single figure
combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2,
                              top = textGrob("Physicochemical Properties of Core Active Compounds", 
                                           gp = gpar(fontsize = 16, fontface = "bold")))

# Save the combined plot
ggsave(file.path(output_dir, "Figure5_ADMET_Properties.png"),
       plot = combined_plot, width = 12, height = 10, dpi = 600, units = "in")

cat("✓ ADMET properties visualization saved to results/figures/Figure5_ADMET_Properties.png\n")


# --- 6. 生成最终统计摘要 ---
admet_summary_text <- paste(
  "ADMET Screening Summary",
  "========================================",
  paste("Total initial compounds screened:", nrow(compounds_screened)),
  paste("Core compounds passing all filters:", nrow(final_core_compounds)),
  paste0("ADMET Success Rate: ", round(nrow(final_core_compounds) / nrow(compounds_screened) * 100, 1), "%"),
  "",
  "Screening Criteria Applied:",
  "---------------------------",
  paste("Molecular Weight (MW): <=", LIPINSKI_MW_MAX, "Da"),
  paste("LogP (ALOGP): <=", LIPINSKI_LOGP_MAX),
  paste("Hydrogen Bond Donors (HBD): <=", LIPINSKI_HBD_MAX),
  paste("Hydrogen Bond Acceptors (HBA): <=", LIPINSKI_HBA_MAX),
  paste("Topological Polar Surface Area (TPSA): <=", TPSA_MAX, "Å²"),
  "========================================"
)

writeLines(admet_summary_text, file.path(output_dir, "ADMET_analysis_summary.txt"))
cat("✓ 统计摘要已生成。
")

cat("
=== ADMET筛选和可视化完成 (V3) ===
")
cat("✓ 对", nrow(compounds_screened), "个初始化合物进行了统一的de novo预测和筛选。
")
cat("✓ 最终获得", nrow(final_core_compounds), "个核心化合物用于所有下游分析。
")
cat("✓ 关键输出文件 'results/admet_analysis_data.csv' 已更新。
")
cat("完成时间:", as.character(Sys.time()), "
")
