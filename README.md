# 🌿 Lithospermum Network Pharmacology Analysis - Reproducibility Guide (Validated)

## 📋 Project Overview

This repository contains the complete computational pipeline for the manuscript: **"Repositioning *Lithospermum erythrorhizon* as a Neuro-Immune Modulator: A Network Pharmacology-Driven Hypothesis Centered on the Serotonergic Synapse"**.

**Repository**: https://github.com/appleinmusic/Lithospermum-Network-Pharmacology  
**Status**: Manuscript ready for submission  
**Last Updated**: July 9, 2025
**Analysis Version**: 4.0 (High-Confidence Sensitivity Analysis Update)

## 🆕 **Latest Update (v4.0) - High-Confidence Validation**
**Date**: July 9, 2025
**Type**: Methodological Enhancement & Hypothesis Validation

### **Key Changes:**
✅ **Added Sensitivity Analysis**: Implemented a two-tiered network analysis to ensure the robustness of our findings. We now generate and analyze both a medium-confidence (score ≥ 400) and a high-confidence (score ≥ 700) network.
✅ **Validated Scientific Conclusion**: The study's central finding—the significant enrichment of the **"Serotonergic synapse" pathway**—is confirmed to be the top hit in **both** networks, providing powerful evidence for our hypothesis.
✅ **Added New Analysis Scripts**: `03a_...` and `06a_...` have been added to the pipeline to perform the high-confidence analysis.
✅ **Updated All Scripts & Results**: The `run_all.sh` script has been updated to execute the complete, validated pipeline.

> **🎯 研究总结 (V4.0 - 经高置信度验证):**
> 本计算流程提供了一个透明且可重复的分析，旨在揭示中药紫草（Lithospermum erythrorhizon）的分子作用机制。本研究的核心是通过一项双层网络敏感性分析，验证了其核心靶点与“血清素能突触”通路的紧密联系。
>
> **方法流程与关键发现:**
> 1.  **数据审查与清理:** 从CMAUP v2.0数据库中确认了 **24** 种真实的紫草来源化合物。
> 2.  **ADMET筛选:** 筛选出 **20** 种核心活性化合物，对应 **32** 个独特的人类蛋白质靶点。
> 3.  **网络构建与验证:** 
>     - **探索性网络 (Score ≥ 400):** 构建了一个包含29个节点和174条边的PPI网络。
>     - **高置信度网络 (Score ≥ 700):** 构建了一个更严格的、包含24个节点和90条边的验证网络。
> 4.  **通路富集与假说确认:** 对两个网络分别进行通路富集分析。结果显示，**“血清素能突触 (Serotonergic synapse)”** 信号通路在两个网络中均排名第一 (高置信度FDR = 1.88e-05)，从而为紫草通过神经-免疫轴发挥作用的核心假说提供了强有力的、可信的证据。

---

## 🔬 Analysis Pipeline Overview

The pipeline is executed by the `scripts/run_all.sh` master script, which runs the R scripts in the correct order. The key analysis stages are now:

1.  **Data Preparation (`01_...`)**: Curation of compounds and targets.
2.  **ADMET Filtering (`02_...`)**: Selection of 20 core active compounds.
3.  **Medium-Confidence Analysis (`03_...`, `06_...`)**: Construction and enrichment analysis of the exploratory network (score ≥ 400).
4.  **High-Confidence Analysis (`03a_...`, `06a_...`)**: Construction and enrichment analysis of the validation network (score ≥ 700).
5.  **Visualization (`04_...`, `05_...`, `07_...`, `08_...`)**: Generation of all figures and network visualizations based on the analysis outputs.

## 🚀 How to Reproduce the Analysis

1.  **Prerequisites**: Ensure you have `R` and `bash` installed.
2.  **Download Data**: The necessary STRING database files are large and should be downloaded separately. Run the `scripts/download_string_db.sh` script to automatically download them into the `data/string_db/` directory.
3.  **Run Pipeline**: Execute the master script from the project root directory:
    ```bash
    bash scripts/run_all.sh
    ```
    This script will first check for and install required R packages, then run the complete, validated analysis pipeline. All results, including those from the high-confidence analysis, will be generated in the `/results` directory.
