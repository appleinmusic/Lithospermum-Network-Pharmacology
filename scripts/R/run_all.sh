#!/bin/bash
# ==============================================================================
# 一键执行脚本：运行完整的紫草网络药理学分析流程
#
# 功能:
#   - 自动定位项目根目录，支持在任何位置执行此脚本。
#   - 按正确顺序 (01 to 08) 依次调用所有R分析脚本。
#   - 提供一个统一的入口点，简化整个分析流程的复现。
#   - 记录每个脚本的开始和结束时间，并报告总耗时。
#
# 使用方法:
#   1. 确保所有 R 脚本 (01-08) 和数据都在正确的目录结构下。
#   2. 在终端中，直接运行此脚本，无需关心当前路径。
#   3. 首次运行前赋予脚本执行权限: chmod +x scripts/R/run_all.sh
#   4. 运行命令: ./scripts/R/run_all.sh
# ==============================================================================

# Function to print a timestamped message
log_msg() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

# --- 1. Environment Sanity Check & Package Installation ---
log_msg "--- [环境检查] 正在验证并安装必需的R包... ---"
R -e '
options(repos = c(CRAN = "https://cran.rstudio.com/"))
required_packages <- c(
    "igraph", "ggplot2", "ggraph", "ggrepel", "dplyr", "RColorBrewer",
    "gridExtra", "grid", "corrplot", "pheatmap", "reshape2", "scales",
    "clusterProfiler", "org.Hs.eg.db", "enrichplot", "readr", "stringr",
    "data.table", "rcdk", "ggforce", "concaveman"
)
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages) > 0) {
    cat("发现未安装的包，正在安装:", paste(new_packages, collapse=", "), "\n")
    install.packages(new_packages)
} else {
    cat("✓ 所有必需的R包均已安装。\n")
}
'

# Check if the R command was successful
if [ $? -ne 0 ]; then
    log_msg "!!! [错误] R包安装或检查失败。请检查您的R环境和网络连接。流程中断。 !!!"
    exit 1
fi
log_msg "--- [环境检查] 完成 ---"


# --- 2. Change to Project Root Directory ---
# Find the directory where this script is located
# SCRIPT_DIR=$(dirname "$(realpath "$0")")
# PROJECT_ROOT=$(realpath "$SCRIPT_DIR/../../")
# cd "$PROJECT_ROOT"

# Fallback for environments without realpath
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/../../" &> /dev/null && pwd)

log_msg "项目根目录: $PROJECT_ROOT"
cd "$PROJECT_ROOT" || { log_msg "无法切换到项目根目录，退出。"; exit 1; }

# --- 3. Execute Analysis Scripts in Order ---
log_msg "=== 开始执行完整的分析流程 ==="
START_TIME=$SECONDS

# List of scripts to run
SCRIPTS=(
    "scripts/R/01_data_preparation.R"
    "scripts/R/02_ADMET_filtering.R"
    "scripts/R/03_network_construction.R"
    "scripts/R/04_network_analysis_and_viz.R"
    "scripts/R/05_module_analysis.R"
    "scripts/R/06_enrichment_analysis.R"
    "scripts/R/07_compound_target_viz.R"
    "scripts/R/08_docking_validation_viz.R"
)

for SCRIPT in "${SCRIPTS[@]}"; do
    if [ -f "$SCRIPT" ]; then
        log_msg "--- [开始] 正在运行: $SCRIPT ---"
        script_start_time=$SECONDS
        Rscript "$SCRIPT"
        if [ $? -ne 0 ]; then
            log_msg "!!! [错误] $SCRIPT 执行失败。流程中断。 !!!"
            exit 1
        fi
        script_end_time=$SECONDS
        log_msg "--- [成功] 完成: $SCRIPT (耗时: $((script_end_time - script_start_time))s) ---"
        echo "" # Add a blank line for readability
    else
        log_msg "!!! [警告] 脚本未找到: $SCRIPT ，已跳过。 !!!"
    fi
done

END_TIME=$SECONDS
log_msg "=== 所有脚本已成功执行完毕 ==="
log_msg "总耗时: $((END_TIME - START_TIME))s"
log_msg "========================================"

exit 0
