#!/usr/bin/env Rscript
# Compound-Target Interaction Network Visualization

cat("=== Compound-Target Interaction Network Visualization ===\n")
cat("Start time:", as.character(Sys.time()), "\n")

# Load necessary packages
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(stringr)
})

# Set working directory and create output directory
if (!file.exists("data") && file.exists("../../zwsjk")) {
  setwd("../..")
}
output_dir <- "results/figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# Define input file path
target_data_file <- "data/processed/lithospermum_targets.tsv"
ingredient_data_file <- "data/processed/lithospermum_ingredients_filtered.tsv"

# Check if input files exist
if (!file.exists(target_data_file)) {
  cat("❌ ERROR: Target data file not found at", target_data_file, "\n")
  stop("Missing target data file.")
}
if (!file.exists(ingredient_data_file)) {
  cat("❌ ERROR: Ingredient data file not found at", ingredient_data_file, "\n")
  stop("Missing ingredient data file.")
}

# Read target data
cat("Reading target data from:", target_data_file, "\n")
target_data <- readr::read_tsv(target_data_file, show_col_types = FALSE)
cat("✓ Loaded", nrow(target_data), "target associations.\n")

# Read ingredient data to get compound names
cat("Reading ingredient data from:", ingredient_data_file, "\n")
ingredient_data <- readr::read_tsv(ingredient_data_file, show_col_types = FALSE) %>%
  select(Ingredient_ID = np_id, Compound_Name = pref_name)
cat("✓ Loaded", nrow(ingredient_data), "ingredients with names.\n")

# Filter for experimentally validated interactions
# Assuming Activity_Value and Activity_Unit are present and indicate experimental validation
# Also ensuring Gene_Symbol is not NA
validated_interactions <- target_data %>%
  filter(!is.na(Activity_Value) & !is.na(Activity_Unit) & Activity_Unit != "" & !is.na(Gene_Symbol)) %>%
  select(Ingredient_ID, Gene_Symbol, Activity_Value, Activity_Unit)

cat("Found", nrow(validated_interactions), "interactions with activity values.\n")

# Create unique compound-target pairs
# We use Ingredient_ID and Gene_Symbol to define a unique pair
unique_pairs <- validated_interactions %>%
  distinct(Ingredient_ID, Gene_Symbol)

num_unique_pairs <- nrow(unique_pairs)
cat("✓ Found", num_unique_pairs, "unique experimentally validated compound-target pairs.\n")
# Compare with manuscript claim
if (num_unique_pairs == 77) {
  cat("  This matches the 77 pairs mentioned in the manuscript.\n")
} else {
  cat("  WARNING: This differs from the 77 pairs mentioned in the manuscript. Manuscript needs update to", num_unique_pairs, ".\n")
}

# Prepare data for network graph
edges <- unique_pairs %>%
  rename(from = Ingredient_ID, to = Gene_Symbol)

# Create node list with types (compound or target) and names
compound_nodes <- ingredient_data %>%
  filter(Ingredient_ID %in% edges$from) %>%
  select(id = Ingredient_ID, name = Compound_Name) %>%
  mutate(type = "Compound", label = str_trunc(name, width = 20, side = "right", ellipsis = "...")) %>%
  distinct(id, .keep_all = TRUE)

target_nodes <- tibble(id = unique(edges$to), name = unique(edges$to)) %>%
  mutate(type = "Target", label = name) %>%
  distinct(id, .keep_all = TRUE)

nodes <- bind_rows(compound_nodes, target_nodes) %>%
  distinct(id, .keep_all = TRUE) # Ensure unique nodes if a compound ID somehow matches a gene symbol

# Create igraph object
interaction_graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

# Simplify graph (remove multiple edges between same nodes if any, and loops)
interaction_graph <- simplify(interaction_graph, remove.multiple = TRUE, remove.loops = TRUE)
cat("Network has", vcount(interaction_graph), "nodes and", ecount(interaction_graph), "edges after simplification.\n")

# Basic network plot using ggraph
set.seed(123) # for reproducibility
p_network <- ggraph(interaction_graph, layout = "fr") + # Fruchterman-Reingold layout
  geom_edge_link(aes(alpha = 0.8), color = "grey50", show.legend = FALSE) +
  geom_node_point(aes(color = type, size = degree(interaction_graph, mode = "all"))) +
  geom_node_text(aes(label = label), repel = TRUE, size = 2.5, max.overlaps = Inf) +
  scale_color_manual(values = c("Compound" = "#FF6B6B", "Target" = "#4ECDC4")) + # Distinct colors
  scale_size_continuous(range = c(3, 10)) + # Node size based on degree
  labs(
    title = "Experimentally Validated Compound-Target Interaction Network",
    subtitle = paste0("Visualizing ", ecount(interaction_graph), " unique interactions between ",
                     length(unique(nodes$id[nodes$type=='Compound'])), " compounds and ",
                     length(unique(nodes$id[nodes$type=='Target'])), " targets"),
    caption = "Data source: Lithospermum erythrorhizon study"
  ) +
  theme_graph(base_family = "sans") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, size=12)
  )

# Save the plot
output_png_path <- file.path(output_dir, "Figure7_Compound_Target_Network.png")
output_pdf_path <- file.path(output_dir, "Figure7_Compound_Target_Network.pdf")

ggsave(filename = output_png_path, plot = p_network, width = 12, height = 10, dpi = 300, units = "in")
ggsave(filename = output_pdf_path, plot = p_network, width = 12, height = 10, units = "in")

cat("✓ Network visualization saved to:", output_png_path, "and", output_pdf_path, "\n")

# Output a summary file with the number of pairs
summary_file_path <- file.path(dirname(output_dir), "compound_target_network_summary.txt") # Save in results/
summary_text <- c(
  paste("Compound-Target Interaction Network Summary (Figure 7)"),
  paste("Date Generated:", Sys.Date()),
  paste("Total unique experimentally validated compound-target pairs:", num_unique_pairs),
  paste("Total compounds in network:", length(unique(nodes$id[nodes$type=='Compound']))),
  paste("Total targets in network:", length(unique(nodes$id[nodes$type=='Target']))),
  paste("Total nodes in graph:", vcount(interaction_graph)),
  paste("Total edges in graph:", ecount(interaction_graph))
)
writeLines(summary_text, summary_file_path)
cat("✓ Network summary saved to:", summary_file_path, "\n")

cat("End time:", as.character(Sys.time()), "\n")
cat("=== Compound-Target Interaction Network Visualization Complete ===\n")

