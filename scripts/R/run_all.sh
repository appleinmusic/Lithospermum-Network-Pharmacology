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

# --- 自动定位并切换到项目根目录 ---
# 获取脚本自身的绝对路径
SCRIPT_PATH=$(realpath "$0")
# 获取脚本所在的目录
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
# 计算项目根目录 (假定此脚本位于 'scripts/R/' 子目录中)
PROJECT_ROOT=$(realpath "$SCRIPT_DIR/../..")

# 切换到项目根目录
cd "$PROJECT_ROOT" || exit 1
# --- 完成目录切换 ---


# 函数：打印带有时间戳的日志信息
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# 确保在项目根目录运行
if [ ! -d "scripts" ] || [ ! -d "data" ]; then
    log "错误：请在项目根目录（包含 scripts/ 和 data/ 的目录）下运行此脚本。"
    exit 1
fi

# 定义要执行的脚本列表
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

# 记录总开始时间
TOTAL_START_TIME=$(date +%s)

log "=== 开始执行完整的分析流程 ==="
echo "========================================"

# 循环执行所有脚本
for script in "${SCRIPTS[@]}"; do
    log "--- [开始] 正在运行: $script ---"
    START_TIME=$(date +%s)
    
    # 运行R脚本
    Rscript "$script"
    
    # 检查R脚本是否成功执行
    if [ $? -ne 0 ]; then
        log "!!! [错误] $script 执行失败。流程中断。 !!!"
        exit 1
    fi
    
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    log "--- [成功] 完成: $script (耗时: ${DURATION}s) ---"
    echo ""
done

TOTAL_END_TIME=$(date +%s)
TOTAL_DURATION=$((TOTAL_END_TIME - TOTAL_START_TIME))

echo "========================================"
log "=== 所有脚本已成功执行完毕 ==="
log "总耗时: ${TOTAL_DURATION}s"
echo "========================================"
