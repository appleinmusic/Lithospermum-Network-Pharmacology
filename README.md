# 🌿 Lithospermum Network Pharmacology Analysis - Complete Reproducibility Guide

## 📋 Project Overview

This repository contains the complete computational pipeline for scientific research on "Network Pharmacology Analysis Reveals the Multi-Target Mechanisms of Lithospermum erythrorhizon Against Inflammatory Diseases: A Comprehensive Study Based on CMAUP v2.0 and STRING Database Integration".

**Repository**: https://github.com/appleinmusic/Lithospermum-Network-Pharmacology  
**Research Topic**: Network Pharmacology Analysis Reveals the Multi-Target Mechanisms of Lithospermum erythrorhizon Against Inflammatory Diseases: A Comprehensive Study Based on CMAUP v2.0 and STRING Database Integration  
**Status**: Research manuscript  
**Last Updated**: June 21, 2025

> **🎯 Research Summary**: This computational pipeline investigates how traditional Chinese medicine Lithospermum erythrorhizon (紫草) works against inflammatory diseases through a sophisticated multi-target network involving TP53, PPARG, EGFR, and PTGS2 hub proteins. Out of 21 bioactive compounds, 20 (95.2%) showed excellent drug-like properties and passed the ADMET screen. Shikonofuran C was excluded due to a predicted LogP value greater than 5.

---

## 🚀 Quick Start (30 seconds)

```bash
# 1. Clone and navigate
git clone https://github.com/appleinmusic/Lithospermum-Network-Pharmacology.git
cd Lithospermum_Network_Pharmacology_Reproducibility

# 2. Download required databases (see setup section for details)
# - CMAUP database files → zwsjk/ directory
# - STRING database files → data/string_db/ directory

# 3. Run complete analysis pipeline
# This single script executes all R scripts (01 to 08) in order.
# Ensure it has execute permissions (chmod +x scripts/R/run_all.sh)
./scripts/R/run_all.sh

# 4. View key results
ls results/figures/                                  # Research figures
ls results/tables/                                   # Data tables
```

> **⚠️ Note**: Both CMAUP and STRING database files are required but not included due to size. See setup instructions below.

---

## 🛠️ System Requirements

### Operating System
- **Primary**: macOS (tested on macOS 14.0+)
- **Alternative**: Linux (Ubuntu 20.04+) or Windows 10+ with WSL2

### R Environment
- **R Version**: 4.3.0 or higher
- **RStudio**: 2023.06.0 or higher (optional but recommended)
- **Memory**: Minimum 8GB RAM, 16GB recommended for network analysis
- **Storage**: 2GB free space for data and results

---

## 📦 Required R Packages

### Core Analysis Packages
```r
# Install required packages (run this once)
install.packages(c(
  "dplyr", "tidyr", "readr", "ggplot2",           # Data manipulation & visualization
  "igraph", "visNetwork", "ggraph",               # Network analysis & visualization  
  "clusterProfiler", "org.Hs.eg.db",             # Pathway enrichment analysis
  "stringr", "purrr", "lubridate",               # Utility packages
  "RColorBrewer", "viridis", "scales",           # Visualization enhancement
  "corrplot", "pheatmap", "VennDiagram",         # Statistical visualization
  "rstudioapi",                                   # RStudio integration
  "rcdk", "rJava"                                 # ADMET analysis
))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c(
  "clusterProfiler",
  "org.Hs.eg.db", 
  "GOSemSim",
  "DOSE"
))
```

> **Java Dependency**: The `rJava` and `rcdk` packages require a properly configured Java Development Kit (JDK). Please ensure Java is installed and `JAVA_HOME` is set correctly before installing these packages. You can check with `R CMD javareconf` in your terminal.

### Package Versions Used in Research
```r
# Verified package versions for exact reproducibility
sessionInfo()
# R version 4.3.0 (2023-04-21)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Ventura 13.4

# Key package versions:
# - dplyr_1.1.3
# - ggplot2_3.4.4  
# - igraph_1.5.1
# - clusterProfiler_4.8.3
# - org.Hs.eg.db_3.17.0
# - stringr_1.5.0
# - visNetwork_2.1.2
# - ggraph_2.1.0
# - rcdk_3.5.0
```

---

## 📂 Repository Structure

```
Lithospermum_Network_Pharmacology_Reproducibility/
├── README.md                     # This file
├── LICENSE                       # MIT License
├── 
├── data/                        # Input data directory
│   ├── processed_data/          # Intermediate processed files
│   └── string_db/              # STRING database files
│
├── zwsjk/                      # CMAUP v2.0 database files
│   ├── CMAUPv2.0_download_Plants.txt
│   ├── CMAUPv2.0_download_Ingredients_onlyActive.txt
│   ├── CMAUPv2.0_download_Targets.txt
│   ├── CMAUPv2.0_download_Plant_Ingredient_Associations_onlyActiveIngredients.txt
│   ├── CMAUPv2.0_download_Ingredient_Target_Associations_ActivityValues_References.txt
│   └── CMAUPv2.0_download_Human_Oral_Bioavailability_information_of_Ingredients_All.txt
│
├── scripts/                    # Analysis scripts
│   └── R/                     # R analysis scripts (01-08)
│       ├── 01_data_preparation.R
│       ├── 02_ADMET_filtering.R
│       ├── 03_network_construction.R
│       ├── 04_network_analysis_and_viz.R
│       ├── 05_module_analysis.R
│       ├── 06_enrichment_analysis.R
│       ├── 07_compound_target_viz.R
│       ├── 08_docking_validation_viz.R
│       └── run_all.sh
│
└── results/                   # Generated outputs
    ├── figures/              # High-quality research figures
    ├── tables/               # Data tables (CSV format)
    └── network/              # Network data files
```

---

## 🚀 Quick Start Guide

### Step 1: Clone the Repository
```bash
git clone https://github.com/appleinmusic/Lithospermum-Network-Pharmacology.git
cd Lithospermum_Network_Pharmacology_Reproducibility
```

### Step 2: Set Up CMAUP Database Files
**Note**: CMAUP database files are not included in the Git repository due to size constraints.

**Option A: Download from CMAUP website**
1. Visit https://www.bidd.group/CMAUP/download.html (official download page)
2. Download the required files listed in the repository structure
3. Place them in the `zwsjk/` directory

**Citation requirements**: If you use CMAUP database in your research, please cite both publications:
- Hou et al. (2024) for the 2024 update: PMID: 37897343; DOI: 10.1093/nar/gkad921
- Zeng et al. (2019) for the original database: PMID: 30357356; DOI: 10.1093/nar/gky965

**Option B: Contact authors**
- Email the corresponding author for pre-processed CMAUP files
- Files will be provided in a compressed format

### Step 3: Set Up STRING Database Files
**Note**: STRING database files are required for protein-protein interaction network construction.

**Download STRING v12.0 data:**
1. Visit https://string-db.org/cgi/download (official STRING download page)
2. Download the following human protein files (Homo sapiens - taxon ID: 9606):
   ```
   https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz
   https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz
   https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz
   https://stringdb-downloads.org/download/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.txt.gz
   ```
3. Place these files in the `data/string_db/` directory (keep them compressed as .gz files)

**Quick download commands:**
```bash
# Create directory
mkdir -p data/string_db

# Download STRING database files
cd data/string_db
wget https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz
wget https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz  
wget https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz
wget https://stringdb-downloads.org/download/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.txt.gz
cd ../..
```

**Citation requirement**: If you use STRING database, please cite:
- Szklarczyk et al. (2023) STRING v12.0. Nucleic Acids Res. DOI: 10.1093/nar/gkac1000

### Step 4: Install R Dependencies
```r
# Run in R console
source("install_dependencies.R")  # If available, or install manually as shown above
```

### Step 5: Execute Analysis Pipeline
```bash
# Navigate to project root
cd /path/to/Lithospermum_Network_Pharmacology_Reproducibility

# Run complete analysis pipeline (single command)
./scripts/R/run_all.sh
```

**Alternative: RStudio Execution**
1. Open the project in RStudio
2. Open and run each script in order (01-08)
3. Results will be generated in the `results/` directory

---

## 📊 Expected Outputs

### Key Research Findings
- **Compound Screening**: From an initial 21 compounds, 20 passed ADMET screening. This set was further filtered to **12 core active compounds** that have experimentally validated protein targets in the CMAUP database, ensuring all subsequent analyses are based on established interactions.
- **32 unique protein targets** identified for the 12 core compounds.
- **High-confidence PPI network**: Constructed with **32 nodes and 156 edges** (network density: 0.321), indicating a highly interconnected functional cluster.
- **7 hub proteins identified**: EGFR, PTGS2, CYP2C9, CYP3A4, ALDH1A1, CYP1A2, and HIF1A, representing the most influential nodes in the network.
- **3 functional modules** identified via Louvain clustering, highlighting distinct clusters for metabolism, signaling, and inflammation.
- **Key Pathway Enrichment**: The "Arachidonic acid metabolism" pathway was identified as the most significant, directly linking the herb's targets to core inflammatory processes.

### Generated Files
After successful execution, the following files should be created:

**Main Results Tables**
- `results/admet_analysis_data.csv` - ADMET properties of filtered compounds
- `results/compound_target_network_summary.txt` - Network summary statistics
- `results/kegg_enrichment_manuscript_table.csv` - Key pathway enrichment results
- `results/molecular_docking_results.csv` - Docking analysis results

**Network Data**
- `results/network/ppi_network.rds` - R object containing the PPI network
- `results/tables/network_topology_stats.csv` - Network topology metrics
- `results/tables/functional_modules.csv` - Module detection results

**Pathway Enrichment**
- `results/go_bp_enrichment_clusterprofiler.csv` - GO Biological Process enrichment
- `results/go_cc_enrichment_clusterprofiler.csv` - GO Cellular Component enrichment  
- `results/go_mf_enrichment_clusterprofiler.csv` - GO Molecular Function enrichment
- `results/kegg_enrichment_clusterprofiler_full.csv` - Complete KEGG pathway results

**Research Figures**
- `results/figures/Figure1_PPI_Network.png` - Main PPI network visualization
- `results/figures/Figure2_Network_Topology.png` - Network topology analysis (4-panel figure)
- `results/figures/Figure2C_Centrality_Correlation.png` - Centrality measures correlation
- `results/figures/Figure4_Functional_Modules_Faceted.png` - Functional module visualization
- `results/figures/Figure5_ADMET_Properties.png` - ADMET properties analysis
- `results/figures/Figure6_Lipinski_Analysis.png` - Drug-likeness analysis
- `results/figures/Figure7_Top_Compounds.png` - Top compounds analysis
- `results/figures/Figure7_Compound_Target_Network.png` - Compound-target network

---

## 🔧 Troubleshooting

### Common Issues

**1. CMAUP Database Files Missing**
```
Error: cannot open file 'zwsjk/CMAUPv2.0_download_Plants.txt'
```
**Solution**: Ensure all CMAUP files are placed in the `zwsjk/` directory

**2. STRING Database Files Missing**
```
Error: cannot open file 'data/string_db/9606.protein.info.v12.0.txt.gz'
```
**Solution**: Download STRING database files as described in Step 3 above and place them in `data/string_db/` directory

**3. R Package Installation Issues**
```
Error: package 'clusterProfiler' is not available
```
**Solution**: Install Bioconductor packages using BiocManager as shown above

**4. Memory Issues During Network Analysis**
```
Error: cannot allocate vector of size X Gb
```
**Solution**: 
- Increase R memory limit: `memory.limit(size=16000)` (Windows)
- Close other applications to free up RAM
- Consider running on a machine with more memory

**5. Working Directory Issues**
```
Error: cannot find project root directory
```
**Solution**: The scripts include automatic working directory detection. If issues persist:
```r
# Set working directory manually
setwd("/path/to/Lithospermum_Network_Pharmacology_Reproducibility")
```

### Performance Optimization

**Parallel Processing** (Optional)
```r
# Enable parallel processing for faster pathway enrichment
library(parallel)
options(mc.cores = detectCores() - 1)
```

**Memory Management**
```r
# Clear workspace between scripts if memory is limited
rm(list = ls())
gc()
```

---

## 📈 Data Processing Details

### Random Seeds
All random processes use fixed seeds for reproducibility:
- Network clustering: `set.seed(42)`
- Permutation tests: `set.seed(123)`
- Bootstrap sampling: `set.seed(456)`

### File Formats
- **Input**: TSV/TXT files from CMAUP database
- **Intermediate**: CSV files for processed data
- **Network**: RDS format for R objects
- **Figures**: PNG (high-resolution) and PDF formats

### Quality Control
Each script includes quality control checks:
- Data completeness verification
- File existence validation  
- Result consistency checks
- Error handling with informative messages

---

## 📞 Support and Contact

### Issues and Bug Reports
- **GitHub Issues**: https://github.com/appleinmusic/Lithospermum-Network-Pharmacology/issues
- **Primary Contact**: [研究团队邮箱]
- **Technical Support**: [技术支持邮箱]

### Citation
If you use this code or data in your research, please cite:
```
[Research manuscript citation information]
Network Pharmacology Analysis Reveals the Multi-Target Mechanisms of Lithospermum erythrorhizon 
Against Inflammatory Diseases: A Comprehensive Study Based on CMAUP v2.0 and STRING Database Integration.
[Journal Name]. [Year]; [Volume]([Issue]): [Pages]. DOI: [DOI]
```

### Supplementary Data Citation
Please also cite the underlying databases:
```
CMAUP Database:
- Hou et al. (2024) CMAUP v2.0. Nucleic Acids Res. DOI: 10.1093/nar/gkad921
- Zeng et al. (2019) CMAUP. Nucleic Acids Res. DOI: 10.1093/nar/gky965

STRING Database:
- Szklarczyk et al. (2023) STRING v12.0. Nucleic Acids Res. DOI: 10.1093/nar/gkac1000
```

### License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 📝 Version History

- **v1.0.0** (June 21, 2025): Complete analysis pipeline with validated results and full reproducibility

---

## 🙏 Acknowledgments

- CMAUP database developers for providing comprehensive TCM data
- STRING database team for protein interaction data
- Bioconductor community for excellent pathway analysis tools
- R community for statistical computing resources

---

**Last Updated**: June 21, 2025  
**Maintained by**: Lithospermum Network Pharmacology Research Team  
**Repository**: https://github.com/appleinmusic/Lithospermum-Network-Pharmacology  
**Status**: ✅ Research Ready - All analysis scripts validated and results verified
