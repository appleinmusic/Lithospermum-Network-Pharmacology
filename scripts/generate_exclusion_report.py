import pandas as pd
import os

# --- 配置 ---
INITIAL_INGREDIENTS_PATH = "data/processed/lithospermum_ingredients_initial.tsv"
EXCLUDED_COMPOUNDS_OUTPUT_PATH = "results/tables/supplementary_excluded_compounds.csv"

# 定义被排除的化合物ID及其原因
# 污染物
contaminant_ids = {
    "NPC159589": "Data Contaminant (Non-plant origin)",
    "NPC234053": "Data Contaminant (Non-plant origin)"
}
# ADMET筛选失败
admet_failed_ids = {
    "NPC113428": "Failed ADMET Screening",
    "NPC144557": "Failed ADMET Screening",
    "NPC216630": "Failed ADMET Screening",
    "NPC308799": "Failed ADMET Screening"
}

all_excluded_ids = {**contaminant_ids, **admet_failed_ids}

# --- 脚本执行 ---
def generate_exclusion_report():
    """
    生成一个详细说明被排除化合物及其原因的报告。
    """
    print(f"开始生成化合物排除报告...")

    # 确保输出目录存在
    os.makedirs(os.path.dirname(EXCLUDED_COMPOUNDS_OUTPUT_PATH), exist_ok=True)

    # 读取初始的26个化合物数据
    try:
        initial_df = pd.read_csv(INITIAL_INGREDIENTS_PATH, sep='\\t', engine='python')
        print(f"成功读取初始成分文件: {INITIAL_INGREDIENTS_PATH}")
    except FileNotFoundError:
        print(f"错误: 初始成分文件未找到 at {INITIAL_INGREDIENTS_PATH}")
        return

    # 筛选出被排除的化合物
    excluded_df = initial_df[initial_df['np_id'].isin(all_excluded_ids.keys())].copy()
    
    # 添加排除原因列
    excluded_df['exclusion_reason'] = excluded_df['np_id'].map(all_excluded_ids)
    
    # 重新排列A了列的顺序，以提高可读性
    cols_to_move = ['np_id', 'pref_name', 'exclusion_reason', 'MW', 'LogP', 'TPSA']
    other_cols = [col for col in excluded_df.columns if col not in cols_to_move]
    excluded_df = excluded_df[cols_to_move + other_cols]

    # 保存报告
    excluded_df.to_csv(EXCLUDED_COMPOUNDS_OUTPUT_PATH, index=False)
    
    print(f"报告已成功保存至: {EXCLUDED_COMPOUNDS_OUTPUT_PATH}")
    print(f"报告中包含 {len(excluded_df)} 个化合物的详细信息。")

if __name__ == "__main__":
    generate_exclusion_report() 