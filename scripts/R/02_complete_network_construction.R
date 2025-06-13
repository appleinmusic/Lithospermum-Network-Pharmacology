# Lithospermum erythrorhizon Network Pharmacology Analysis - Network Construction
# Protein-protein interaction network construction using STRING database

# Set random seeds for reproducibility
set.seed(42)  # For network clustering and other random processes

# Load required packages
suppressMessages({
  library(dplyr)
  library(igraph)
  library(data.table)
  library(VennDiagram)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(corrplot)
  library(gridExtra)
  library(stringr)
  library(readr)
})

# Set working directory to project root
if(require(rstudioapi) && rstudioapi::isAvailable()) {
  # Running in RStudio
  project_root <- file.path(dirname(dirname(dirname(rstudioapi::getActiveDocumentContext()$path))))
  setwd(project_root)
} else {
  # Running from command line - assume we're already in project root
  # or set to a known location
  if(basename(getwd()) == "R") {
    setwd("../../")
  } else if(basename(getwd()) == "scripts") {
    setwd("../")
  }
  # If running from project root, do nothing
}
  setwd(project_root)
}

# Create output directories
dir.create("results/network", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

cat("Starting comprehensive network construction analysis...\n")

# 1. Load CMAUP data
cat("Loading CMAUP data from zwsjk directory...\n")

# Load plant data
plants <- read_tsv("zwsjk/CMAUPv2.0_download_Plants.txt", show_col_types = FALSE)

# Find Lithospermum
lithospermum <- plants %>%
  filter(str_detect(tolower(Species_Name), "lithospermum") & 
         str_detect(tolower(Species_Name), "erythrorhizon"))

cat("Found Lithospermum, Plant ID:", lithospermum$Plant_ID, "\n")

# Load plant-ingredient association data
plant_ingredients <- read_tsv("zwsjk/CMAUPv2.0_download_Plant_Ingredient_Associations_onlyActiveIngredients.txt", 
                             col_names = c("Plant_ID", "Ingredient_ID"), 
                             show_col_types = FALSE)

# Filter Lithospermum ingredients
zicao_ingredient_ids <- plant_ingredients %>%
  filter(Plant_ID == lithospermum$Plant_ID) %>%
  pull(Ingredient_ID)

cat("Lithospermum-related ingredients:", length(zicao_ingredient_ids), "\n")

# Load ingredient-target association data
ingredient_targets <- read_tsv("zwsjk/CMAUPv2.0_download_Ingredient_Target_Associations_ActivityValues_References.txt", 
                              show_col_types = FALSE)

# Filter targets of Lithospermum ingredients
zicao_targets_data <- ingredient_targets %>%
  filter(Ingredient_ID %in% zicao_ingredient_ids)

cat("Lithospermum ingredient-target associations:", nrow(zicao_targets_data), "\n")

# Load target detailed information
targets <- read_tsv("zwsjk/CMAUPv2.0_download_Targets.txt", show_col_types = FALSE)

# Merge target information
zicao_targets_with_info <- zicao_targets_data %>%
  left_join(targets, by = "Target_ID") %>%
  filter(!is.na(Gene_Symbol))

cat("Valid Lithospermum targets:", length(unique(zicao_targets_with_info$Gene_Symbol)), "\n")

# Prepare disease-related targets (assuming disease data or using known psoriasis-related genes)
# Here we first use Lithospermum targets as the main analysis object
all_targets <- unique(zicao_targets_with_info$Gene_Symbol)
all_targets <- all_targets[!is.na(all_targets) & all_targets != ""]

cat("Total targets:", length(all_targets), "\n")

# 3. Load STRING database files
cat("Loading STRING database...\n")

# Load protein information
protein_info <- fread("data/string_db/9606.protein.info.v12.0.txt.gz")
cat("STRING protein information loaded successfully, total", nrow(protein_info), "records\n")

# Load protein aliases
protein_aliases <- fread("data/string_db/9606.protein.aliases.v12.0.txt.gz")
cat("STRING protein aliases loaded successfully, total", nrow(protein_aliases), "records\n")

# Check if protein interaction file exists
if (file.exists("data/string_db/9606.protein.links.v12.0.txt.gz")) {
  protein_links <- fread("data/string_db/9606.protein.links.v12.0.txt.gz")
  cat("STRING protein interaction data loaded successfully, total", nrow(protein_links), "records\n")
  
  # Filter high-quality interactions (score >= 400)
  high_quality_links <- protein_links[protein_links$combined_score >= 400, ]
  cat("High-quality interactions (score>=400):", nrow(high_quality_links), "records\n")
  
} else {
  cat("Warning: Protein interaction data file not found\n")
  cat("Please download 9606.protein.links.v12.0.txt.gz and place it in data/string_db/ directory\n")
  cat("Download link: https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz\n")
  stop("Missing required STRING data file")
}

# 4. Map target names to STRING ID
cat("Performing target name mapping...\n")

# Check column names
cat("STRING alias file column names:", colnames(protein_aliases), "\n")

# Create mapping function
map_to_string_id <- function(target_names, aliases_df, protein_info_df) {
  # Ensure using correct column names
  protein_id_col <- colnames(aliases_df)[1]  # First column is protein ID column
  alias_col <- "alias"
  
  cat("Using protein ID column name:", protein_id_col, "\n")
  
  mapped_data <- data.frame(
    target_name = character(0),
    string_id = character(0),
    source = character(0),
    stringsAsFactors = FALSE
  )
  
  for (target in target_names) {
    if (is.na(target) || target == "") next
    
    found <- FALSE
    
    # 1. First try to find in protein info file's preferred_name
    protein_matches <- protein_info_df[protein_info_df$preferred_name == target, ]
    
    if (nrow(protein_matches) > 0) {
      protein_id <- protein_matches[[colnames(protein_info_df)[1]]][1]
      if (length(protein_id) > 0 && !is.na(protein_id) && nchar(as.character(protein_id)) > 0) {
        new_row <- data.frame(
          target_name = target,
          string_id = as.character(protein_id),
          source = "preferred_name",
          stringsAsFactors = FALSE
        )
        mapped_data <- rbind(mapped_data, new_row)
        cat("Successfully mapped:", target, "->", as.character(protein_id), "(source: preferred_name)\n")
        found <- TRUE
      }
    }
    
    # 2. If not found in preferred_name, try alias file
    if (!found) {
      # Exact match gene symbol
      matches <- aliases_df[aliases_df[[alias_col]] == target, ]
      
      if (nrow(matches) == 0) {
        # Try case-insensitive exact match
        matches <- aliases_df[toupper(aliases_df[[alias_col]]) == toupper(target), ]
      }
      
      if (nrow(matches) == 0) {
        # Try partial match
        matches <- aliases_df[grepl(paste0("^", target, "$"), aliases_df[[alias_col]], ignore.case = TRUE), ]
      }
      
      if (nrow(matches) > 0) {
        # Check if protein ID exists and is not NA
        protein_id <- matches[[protein_id_col]][1]
        
        if (length(protein_id) > 0 && !is.na(protein_id) && nchar(as.character(protein_id)) > 0) {
          new_row <- data.frame(
            target_name = target,
            string_id = as.character(protein_id),
            source = paste0("alias:", matches$source[1]),
            stringsAsFactors = FALSE
          )
          mapped_data <- rbind(mapped_data, new_row)
          cat("Successfully mapped:", target, "->", as.character(protein_id), "(source: alias -", matches$source[1], ")\n")
          found <- TRUE
        } else {
          cat("Target found in alias but protein ID is empty, matches:", nrow(matches), "\n")
        }
      }
    }
    
    if (!found) {
      cat("Mapping not found:", target, "\n")
    }
  }
  
  return(mapped_data)
}

# Execute mapping
target_mapping <- map_to_string_id(all_targets, protein_aliases, protein_info)
cat("Successfully mapped targets:", nrow(target_mapping), "/", length(all_targets), "\n")

# Display successfully mapped targets
if (nrow(target_mapping) > 0) {
  cat("Successfully mapped targets:\n")
  print(target_mapping)
}

# 5. Build PPI Network
cat("Building protein-protein interaction network...\n")

if (nrow(target_mapping) > 0) {
  # Get interactions for mapped targets
  target_string_ids <- target_mapping$string_id
  
  # Filter interactions between target proteins
  target_interactions <- high_quality_links[
    high_quality_links$protein1 %in% target_string_ids & 
    high_quality_links$protein2 %in% target_string_ids, 
  ]
  
  cat("Target protein interactions:", nrow(target_interactions), "\n")
  
  if (nrow(target_interactions) > 0) {
    # Create network graph
    g <- graph_from_data_frame(
      target_interactions[, c("protein1", "protein2", "combined_score")],
      directed = FALSE
    )
    
    # Add node attributes
    # Map STRING ID back to target names
    id_to_name <- setNames(target_mapping$target_name, target_mapping$string_id)
    
    V(g)$name_display <- ifelse(names(V(g)) %in% names(id_to_name), 
                               id_to_name[names(V(g))], 
                               names(V(g)))
    
    # Mark node types
    V(g)$node_type <- "Lithospermum_Target"
    
    # 计算网络拓扑参数
    topo_stats <- data.frame(
      节点数 = vcount(g),
      边数 = ecount(g),
      平均度 = mean(degree(g)),
      网络密度 = edge_density(g),
      聚集系数 = transitivity(g),
      平均路径长度 = ifelse(is_connected(g), mean_distance(g), NA),
      直径 = ifelse(is_connected(g), diameter(g), NA),
      连通分量数 = components(g)$no
    )
    
    cat("网络拓扑统计:\n")
    print(topo_stats)
    
    # 保存网络数据
    write.csv(target_mapping, "results/tables/target_string_mapping.csv", row.names = FALSE)
    write.csv(target_interactions, "results/tables/ppi_interactions.csv", row.names = FALSE)
    write.csv(topo_stats, "results/tables/network_topology_stats.csv", row.names = FALSE)
    
    # 保存网络对象
    saveRDS(g, "results/network/ppi_network.rds")
    
    cat("网络构建完成！\n")
    cat("结果已保存到 results/ 目录\n")
    
  } else {
    cat("警告: 在目标蛋白质之间未找到足够的相互作用\n")
  }
  
} else {
  cat("错误: 未能映射任何靶点到STRING数据库\n")
}

cat("网络构建分析完成！\n") 