#!/bin/bash
# Master script to run the entire network pharmacology analysis pipeline
# V5 - Robust, self-contained R dependency check

# --- Configuration ---
# Set the base directory to the project root where this script is located's parent
cd "$(dirname "$0")/.."

LOG_DIR="logs"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/run_all_${TIMESTAMP}.log"

# R scripts to be executed in order, with correct relative paths from project root
R_SCRIPTS=(
  "scripts/R/01_data_preparation.R"
  "scripts/R/02_ADMET_filtering.R"
  "scripts/R/03_network_construction.R"
  "scripts/R/04_network_analysis_and_viz.R"
  "scripts/R/05_module_analysis.R"
  "scripts/R/06_enrichment_analysis.R"
  "scripts/R/07_compound_target_viz.R"
  "scripts/R/08_docking_validation_viz.R"
  "scripts/R/09_interactive_network_viz.R"
)

# --- Setup ---
mkdir -p "$LOG_DIR"
exec &> >(tee -a "$LOG_FILE") # Redirect stdout and stderr to log file and console
echo "=== Log started at $(date) ==="
echo "Project directory: $(pwd)"


# --- R Package Dependency Check ---
echo "--- [INFO] Checking for required R packages... ---"
R -e '
  # Required packages are defined directly in R for robustness
  required_packages <- c(
    "igraph", "ggplot2", "dplyr", "readr", "gridExtra", "ggraph", "RColorBrewer",
    "ggrepel", "tidygraph", "ggforce", "concaveman", "clusterProfiler", "org.Hs.eg.db",
    "enrichplot", "ggplotify", "cowplot", "visNetwork", "htmlwidgets", "rcdk"
  );

  # Identify Bioconductor packages
  bio_pkg_list <- intersect(required_packages, c("clusterProfiler", "org.Hs.eg.db", "enrichplot"));
  
  installed_pkgs <- rownames(installed.packages());
  
  # Check CRAN packages first
  cran_pkgs_to_check <- setdiff(required_packages, bio_pkg_list)
  missing_cran_pkgs <- setdiff(cran_pkgs_to_check, installed_pkgs)
  
  if (length(missing_cran_pkgs) > 0) {
    cat("--- [INFO] The following required CRAN packages are missing:", paste(missing_cran_pkgs, collapse=", "), "\n");
    cat("--- [INFO] Attempting to install missing CRAN packages...\n");
    install.packages(missing_cran_pkgs, repos="https://cloud.r-project.org/");
  } else {
    cat("--- [SUCCESS] All required CRAN packages are already installed.\n");
  }

  # Check Bioconductor packages
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos="https://cloud.r-project.org/")
  }
  
  # Re-check installed packages after potential installations
  installed_pkgs <- rownames(installed.packages());
  missing_bio_pkgs <- setdiff(bio_pkg_list, installed_pkgs)
  if (length(missing_bio_pkgs) > 0) {
      cat("--- [INFO] The following required Bioconductor packages are missing:", paste(missing_bio_pkgs, collapse=", "), "\n");
      cat("--- [INFO] Attempting to install missing Bioconductor packages...\n");
      BiocManager::install(missing_bio_pkgs, update=FALSE, ask=FALSE);
  } else {
      cat("--- [SUCCESS] All required Bioconductor packages are already installed.\n");
  }

  # Final check to ensure all installations were successful
  installed_after <- rownames(installed.packages())
  final_missing <- setdiff(required_packages, installed_after)
  if (length(final_missing) > 0) {
      cat("--- [ERROR] Failed to install the following packages:", paste(final_missing, collapse=", "), ". Please install them manually and restart the script.\n");
      q(status=1)
  } else {
      cat("--- [SUCCESS] All R package dependencies are satisfied.\n");
  }
'

if [ $? -ne 0 ]; then
  echo "--- [FATAL] R package installation failed. Aborting pipeline. ---"
  exit 1
fi


# --- Main Execution Loop ---
echo -e "\n\n=== [INFO] Starting complete analysis pipeline at $(date) ==="
START_TIME=$SECONDS

for script in "${R_SCRIPTS[@]}"; do
  if [ -f "$script" ]; then
    echo -e "\n--- [START] Running: $script ---"
    SCRIPT_START_TIME=$SECONDS
    
    Rscript "$script"
    
    if [ $? -eq 0 ]; then
      SCRIPT_END_TIME=$SECONDS
      ELAPSED_TIME=$((SCRIPT_END_TIME - SCRIPT_START_TIME))
      echo "--- [SUCCESS] Finished: $script (Duration: ${ELAPSED_TIME}s) ---"
    else
      echo -e "\n\n!!! [ERROR] Execution of $script failed. Pipeline aborted. !!!"
      echo "Check the output above and the log file for details: $LOG_FILE"
      exit 1
    fi
  else
    echo "--- [WARNING] Script not found, skipping: $script ---"
  fi
done

END_TIME=$SECONDS
TOTAL_DURATION=$((END_TIME - START_TIME))

echo -e "\n\n=== [SUCCESS] Entire analysis pipeline completed successfully. ==="
echo "Total execution time: ${TOTAL_DURATION} seconds."
echo "All results are in the 'results/' directory."
echo "Full log available at: $LOG_FILE"
echo "=== Log ended at $(date) ===" 