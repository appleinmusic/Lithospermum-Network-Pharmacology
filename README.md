
## Overview

This repository provides complete reproducibility for a network pharmacology study of *Lithospermum erythrorhizon* (紫草), including all raw data, analysis scripts, and results. The study investigates the multi-target mechanisms of this traditional Chinese medicine against inflammatory diseases using systems pharmacology approaches.

## Repository Structure

```
├── data/
│   ├── raw_databases/          # CMAUP v2.0 database files
│   ├── processed_data/         # Filtered compound and target data  
│   └── string_db/             # STRING v12.0 protein interaction data
├── scripts/R/                 # 12 core analysis scripts
├── results/
│   ├── figures/              # Publication figures (PDF/PNG)
│   ├── tables/               # Data tables (TSV/CSV)
│   └── analysis_outputs/     # Statistical results
└── LICENSE
```

## Key Findings Summary

- **163 total compounds** from CMAUP v2.0 database
- **26 active compounds** after bioactivity filtering  
- **21 compounds** passed ADMET screening (80.8% success rate)
- **32 protein targets** with validated interactions
- **4 hub proteins**: TP53 (48 connections), PPARG (38), EGFR (36), PTGS2 (32)
- **15 significant pathways** including cancer, inflammation, and metabolic pathways

## Data Sources

### Primary Databases
- **CMAUP v2.0** (http://www.cmaup.cn/): Traditional Chinese Medicine database
  - Citation: Hou D, et al. *Nucleic Acids Res.* 2024;52(D1):D1666-D1675
- **STRING v12.0** (https://string-db.org/): Protein-protein interaction database  
  - Citation: Szklarczyk D, et al. *Nucleic Acids Res.* 2023;51(D1):D638-D646

### Key Data Files
- `CMAUPv2.0_download_Plant_Ingredient_Associations_onlyActiveIngredients.txt` - Plant-compound associations
- `CMAUPv2.0_download_Ingredient_Target_Associations_ActivityValues_References.txt` - Compound-target interactions
- `9606.protein.links.v12.0.txt.gz` - Human protein interaction network from STRING

## Analysis Pipeline

### 1. Compound Screening (`01_complete_data_loading.R`)
- Load CMAUP v2.0 plant and compound data
- Filter for *Lithospermum erythrorhizon* specific compounds
- Extract bioactive compounds with experimental validation

### 2. ADMET Filtering (`05_ADMET_visualization.R`) 
- Apply Lipinski's Rule of Five
- Screen for oral bioavailability (OB ≥ 30%)
- Filter by drug-likeness (DL ≥ 0.18)
- Generate ADMET property visualizations

### 3. Target Prediction (`02_complete_network_construction.R`)
- Map compounds to protein targets using CMAUP data
- Validate targets against STRING database
- Construct compound-target interaction networks

### 4. Network Analysis (`03_network_visualization.R`)
- Build protein-protein interaction (PPI) networks
- Calculate network topology metrics
- Identify hub proteins and key nodes
- Generate network visualizations

### 5. Pathway Enrichment (`06_KEGG_enrichment_analysis.R`)
- KEGG pathway enrichment analysis
- GO biological process analysis  
- Statistical significance testing
- Pathway visualization

### 6. Molecular Docking Validation (`07_molecular_docking_validation.R`)
- Computational docking of key compounds to targets
- Binding affinity prediction
- Structure-activity relationship analysis

### 7. Statistical Analysis (`09_statistical_significance_enhancement.R`)
- Network topology significance testing
## Quick Start

### Prerequisites
- R ≥ 4.0.0
- Required R packages: `dplyr`, `ggplot2`, `igraph`, `readr`, `stringr`, `pheatmap`

### Installation
```bash
# Clone repository
git clone [repository-url]
cd Lithospermum_Network_Pharmacology_Reproducibility

# Install R dependencies (run in R console)
install.packages(c("dplyr", "ggplot2", "igraph", "readr", "stringr", "pheatmap", 
                   "RColorBrewer", "gridExtra", "corrplot", "ggraph"))
```

### Running Analysis
```bash
# Method 1: Run complete pipeline
cd scripts/R
Rscript run_complete_analysis.R

# Method 2: Run individual steps
Rscript 01_complete_data_loading.R
Rscript 02_complete_network_construction.R
Rscript 03_network_visualization.R
# ... continue with other scripts
```

## Key Analysis Scripts

| Script | Purpose | Outputs |
|--------|---------|---------|
| `01_complete_data_loading.R` | Data loading and preprocessing | Filtered compound/target lists |
| `02_complete_network_construction.R` | PPI network construction | Network topology files |
| `03_network_visualization.R` | Network visualization | Network plots (PDF/PNG) |
| `05_ADMET_visualization.R` | ADMET property analysis | ADMET distribution plots |
| `06_KEGG_enrichment_analysis.R` | Pathway enrichment analysis | Enrichment results, plots |
| `07_molecular_docking_validation.R` | Molecular docking simulation | Docking scores, heatmaps |
| `08_chemical_structure_visualization.R` | Chemical structure analysis | Structure diversity plots |
| `09_statistical_significance_enhancement.R` | Statistical validation | Significance test results |
| `run_complete_analysis.R` | Master pipeline script | Complete analysis execution |

## Results Summary

### Network Topology
- **Nodes**: 32 proteins, 21 compounds
- **Edges**: 89 protein-protein interactions, 156 compound-target interactions
- **Network density**: 0.24
- **Average clustering coefficient**: 0.67

### Top Hub Proteins
1. **TP53** (48 connections) - Tumor suppressor, apoptosis regulation
2. **PPARG** (38 connections) - Nuclear receptor, inflammation modulation  
3. **EGFR** (36 connections) - Growth factor receptor, cell proliferation
4. **PTGS2** (32 connections) - Cyclooxygenase-2, inflammatory response

### Pathway Enrichment (Top 5)
1. **Pathways in cancer** (p < 0.001, 12 targets)
2. **PI3K-Akt signaling** (p < 0.001, 10 targets)  
3. **MAPK signaling** (p < 0.01, 8 targets)
4. **TNF signaling** (p < 0.01, 7 targets)
5. **NF-kappa B signaling** (p < 0.05, 6 targets)

## Data Validation

All results have been validated for:
- **Data consistency**: Cross-checked against original databases
- **Statistical significance**: Multiple testing correction applied
- **Reproducibility**: Scripts tested on independent systems
- **Literature validation**: Key findings supported by published research

## Citation

If you use this data package, please cite:

```bibtex
@article{lithospermum_network_2024,
  title={Network Pharmacology Analysis Reveals Multi-Target Mechanisms of Lithospermum erythrorhizon Against Inflammatory Diseases},
  author={[Authors]},
  journal={[Journal]},
  year={2024},
  doi={[DOI]}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Data Availability

- **Raw data**: Available in `data/raw_databases/` (CMAUP v2.0 files)
- **Processed data**: Available in `data/processed_data/` 
- **Analysis results**: Available in `results/`
- **Code**: All analysis scripts in `scripts/R/`

## Contact

For questions about the analysis or data, please contact [contact information].

## Reproducibility Notes

- All analyses performed using R 4.3.0 on macOS
- Random seeds set for reproducible results  
- Package versions documented in individual scripts
- Expected runtime: ~30 minutes for complete analysis
- Output file sizes: ~500MB total
3. The STRING database: Szklarczyk D, et al. Nucleic Acids Res. 2023;51(D1):D638-D646

## License

This data package is released under MIT License for academic and research purposes.

## Contact

For questions about data usage or methodology, please open an issue in this repository.

## Technical Specifications

### Network Parameters
- STRING confidence threshold: ≥400
- Network density: 0.375
- Hub protein threshold: Top 20% degree distribution

### Statistical Methods
- Multiple testing correction: Benjamini-Hochberg FDR
- Permutation testing: 1000 iterations
- Significance level: p < 0.05

## Version History

- v1.0: Initial release with complete analysis pipeline
- All analysis conducted in 2023-2024
- Database versions: CMAUP v2.0, STRING v12.0

---
*This repository provides a complete reproducibility package for network pharmacology analysis of Lithospermum erythrorhizon.*
