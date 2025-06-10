#!/usr/bin/env Rscript
# ADMET Properties Visualization for Lithospermum erythrorhizon

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
  library(gridExtra)
  library(viridis)
})

cat("Start time:", as.character(Sys.time()), "\n")

# Set directory paths
setwd("/Users/lgmoon/Desktop/zdhky")
output_dir <- "results/figures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Load filtered active ingredients
active_ingredients <- readr::read_tsv("results/tables/lithospermum_active_ingredients_admet_lipinski.tsv")
cat("Loaded active ingredients:", nrow(active_ingredients), "\n")

# 1. ADMET screening standards visualization
admet_standards <- data.frame(
  Property = c("OB (%)", "DL", "MW (Da)", "LogP", "nHA", "nHD"),
  Threshold = c(30, 0.18, 500, 5, 10, 5),
  Unit = c("%", "", "Da", "", "", ""),
  Type = c("≥", "≥", "≤", "≤", "≤", "≤")
)

# 2. Active ingredient ADMET property distribution
p1 <- ggplot(active_ingredients, aes(x = OB_篩選用, y = DL_篩選用)) +
  geom_point(color = "#2E86AB", size = 3, alpha = 0.8) +
  geom_hline(yintercept = 0.18, linetype = "dashed", color = "#E74C3C", linewidth = 0.8) +
  geom_vline(xintercept = 30, linetype = "dashed", color = "#E74C3C", linewidth = 0.8) +
  labs(
    title = "ADMET Screening: OB vs DL Distribution",
    subtitle = "Oral Bioavailability and Drug-likeness of Active Compounds",
    x = "Oral Bioavailability (OB, %)",
    y = "Drug-likeness (DL)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    legend.position = "none"
  ) +
  annotate("text", x = 35, y = 0.16, label = "Screening\nThreshold", 
           color = "#E74C3C", size = 3, hjust = 0, fontface = "bold")

# 3. Lipinski's Rule of Five 评估
if (all(c("MW.x", "LogP", "nHA", "nHD") %in% names(active_ingredients))) {
  lipinski_data <- active_ingredients %>%
    mutate(
      MW_pass = MW.x <= 500,
      LogP_pass = LogP <= 5,
      nHA_pass = nHA <= 10,
      nHD_pass = nHD <= 5,
      lipinski_violations = 4 - (MW_pass + LogP_pass + nHA_pass + nHD_pass)
    )
  
  p2 <- ggplot(lipinski_data, aes(x = factor(lipinski_violations))) +
    geom_bar(fill = "#3498DB", alpha = 0.8, width = 0.6, color = "black", linewidth = 0.3) +
    labs(
      title = "Lipinski's Rule of Five Compliance",
      subtitle = "Distribution of Rule Violations Among Active Compounds",
      x = "Number of Lipinski Violations",
      y = "Number of Compounds"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      panel.background = element_rect(fill = "white", color = "black"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major.y = element_line(color = "#E5E5E5", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5)
    ) +
    scale_x_discrete(labels = c("0", "1", "2", "3", "4")[1:length(unique(lipinski_data$lipinski_violations))])
} else {
  # 如果没有完整的Lipinski数据，创建简化版本
  top_compounds_for_plot <- active_ingredients %>%
    arrange(desc(OB_篩選用)) %>%
    slice_head(n = 15) %>%
    mutate(pref_name = ifelse(nchar(pref_name) > 20, 
                              paste0(substr(pref_name, 1, 17), "..."), 
                              pref_name))
  
  p2 <- ggplot(top_compounds_for_plot, aes(x = reorder(pref_name, OB_篩選用), y = OB_篩選用)) +
    geom_col(fill = "#3498DB", alpha = 0.8, color = "black", linewidth = 0.2) +
    coord_flip() +
    labs(
      title = "Oral Bioavailability of Top Active Compounds",
      subtitle = "OB Values for Top 15 Screened Active Ingredients",
      x = "Active Compounds",
      y = "Oral Bioavailability (%)"
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
}

# 4. Drug property radar chart data preparation
if (nrow(active_ingredients) > 0) {
  # Select representative compounds for radar chart display
  top_compounds <- active_ingredients %>%
    arrange(desc(OB_篩選用)) %>%
    slice_head(n = 10) %>%
    mutate(pref_name = ifelse(nchar(pref_name) > 25, 
                              paste0(substr(pref_name, 1, 22), "..."), 
                              pref_name))
  
  p3 <- ggplot(top_compounds, aes(x = reorder(pref_name, OB_篩選用), y = OB_篩選用)) +
    geom_col(fill = "#E67E22", alpha = 0.8, color = "black", linewidth = 0.2) +
    geom_text(aes(label = round(OB_篩選用, 1)), hjust = -0.1, size = 3, fontface = "bold") +
    coord_flip() +
    labs(
      title = "Top 10 Compounds by Oral Bioavailability",
      subtitle = "Ranking of Most Promising Active Ingredients",
      x = "Compound Name",
      y = "Oral Bioavailability (%)"
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
}

# 5. Save charts
ggsave(file.path(output_dir, "Figure5_ADMET_Properties.png"), 
       plot = p1, width = 8, height = 6, dpi = 600, units = "in")
ggsave(file.path(output_dir, "Figure5_ADMET_Properties.pdf"), 
       plot = p1, width = 8, height = 6, units = "in")

ggsave(file.path(output_dir, "Figure6_Lipinski_Analysis.png"), 
       plot = p2, width = 8, height = 6, dpi = 600, units = "in")
ggsave(file.path(output_dir, "Figure6_Lipinski_Analysis.pdf"), 
       plot = p2, width = 8, height = 6, units = "in")

if (exists("p3")) {
  ggsave(file.path(output_dir, "Figure7_Top_Compounds.png"), 
         plot = p3, width = 8, height = 6, dpi = 600, units = "in")
  ggsave(file.path(output_dir, "Figure7_Top_Compounds.pdf"), 
         plot = p3, width = 8, height = 6, units = "in")
}

# 6. Generate ADMET statistical summary
admet_summary <- list(
  total_compounds = nrow(active_ingredients),
  ob_mean = round(mean(active_ingredients$OB_篩選用, na.rm = TRUE), 2),
  ob_range = paste(round(min(active_ingredients$OB_篩選用, na.rm = TRUE), 2), 
                   round(max(active_ingredients$OB_篩選用, na.rm = TRUE), 2), sep = "-"),
  dl_mean = round(mean(active_ingredients$DL_篩選用, na.rm = TRUE), 3),
  dl_range = paste(round(min(active_ingredients$DL_篩選用, na.rm = TRUE), 3), 
                   round(max(active_ingredients$DL_篩選用, na.rm = TRUE), 3), sep = "-"),
  screening_success_rate = paste0(round(nrow(active_ingredients) / 26 * 100, 1), "%")
)

# Save statistical summary
writeLines(
    paste("Total Active Compounds:", admet_summary$total_compounds),
    paste("OB Mean (%):", admet_summary$ob_mean),
    paste("OB Range (%):", admet_summary$ob_range),
    paste("DL Mean:", admet_summary$dl_mean),
    paste("DL Range:", admet_summary$dl_range),
    paste("Screening Success Rate:", admet_summary$screening_success_rate)),
  file.path(output_dir, "ADMET_analysis_summary.txt")
)

cat("✓ Generated charts:", ifelse(exists("p3"), 3, 2), "个\n")
cat("✓ 保存路径:", output_dir, "\n")
cat("End time:", as.character(Sys.time()), "\n") 