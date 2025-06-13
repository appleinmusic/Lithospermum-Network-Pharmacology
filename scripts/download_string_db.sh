#!/bin/bash
# STRING Database Download Script for Lithospermum Network Pharmacology Project
# This script downloads the required STRING v12.0 database files for human protein interactions

echo "=== STRING Database Download Script ==="
echo "Downloading STRING v12.0 database files for human proteins (Homo sapiens, tax ID: 9606)"
echo "Total download size: ~500MB (compressed)"
echo ""

# Create the directory if it doesn't exist
mkdir -p data/string_db
cd data/string_db

# Base URL for STRING downloads
BASE_URL="https://stringdb-downloads.org/download"

# Download URLs for STRING v12.0
PROTEIN_LINKS_URL="${BASE_URL}/protein.links.v12.0/9606.protein.links.v12.0.txt.gz"
PROTEIN_INFO_URL="${BASE_URL}/protein.info.v12.0/9606.protein.info.v12.0.txt.gz"
PROTEIN_ALIASES_URL="${BASE_URL}/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz"

# Function to download with progress and error checking
download_file() {
    local url=$1
    local filename=$2
    local description=$3
    
    echo "Downloading ${description}..."
    echo "URL: ${url}"
    
    if command -v wget >/dev/null 2>&1; then
        wget --progress=bar:force:noscroll -O "${filename}" "${url}"
    elif command -v curl >/dev/null 2>&1; then
        curl -L --progress-bar -o "${filename}" "${url}"
    else
        echo "Error: Neither wget nor curl is available. Please install one of them."
        exit 1
    fi
    
    if [ $? -eq 0 ]; then
        echo "✓ Successfully downloaded ${filename}"
        echo "File size: $(du -h "${filename}" | cut -f1)"
        echo ""
    else
        echo "✗ Failed to download ${filename}"
        echo "Please check your internet connection and try again."
        exit 1
    fi
}

# Check if files already exist
if [ -f "9606.protein.links.v12.0.txt.gz" ] && \
   [ -f "9606.protein.info.v12.0.txt.gz" ] && \
   [ -f "9606.protein.aliases.v12.0.txt.gz" ]; then
    echo "STRING database files already exist. Checking file integrity..."
    
    # Check if files are not empty
    for file in "9606.protein.links.v12.0.txt.gz" "9606.protein.info.v12.0.txt.gz" "9606.protein.aliases.v12.0.txt.gz"; do
        if [ ! -s "$file" ]; then
            echo "Warning: $file exists but is empty. Re-downloading..."
            rm -f "$file"
            break
        fi
    done
fi

# Download protein interactions file (~400MB)
if [ ! -f "9606.protein.links.v12.0.txt.gz" ]; then
    download_file "${PROTEIN_LINKS_URL}" "9606.protein.links.v12.0.txt.gz" "Protein-Protein Interactions"
fi

# Download protein information file (~50MB)
if [ ! -f "9606.protein.info.v12.0.txt.gz" ]; then
    download_file "${PROTEIN_INFO_URL}" "9606.protein.info.v12.0.txt.gz" "Protein Information"
fi

# Download protein aliases file (~50MB)
if [ ! -f "9606.protein.aliases.v12.0.txt.gz" ]; then
    download_file "${PROTEIN_ALIASES_URL}" "9606.protein.aliases.v12.0.txt.gz" "Protein Aliases"
fi

echo "=== Download Summary ==="
echo "All STRING database files have been downloaded successfully!"
echo ""
echo "Downloaded files:"
ls -lh *.txt.gz

echo ""
echo "Total size:"
du -sh .

echo ""
echo "✓ STRING database setup complete!"
echo "You can now run the network pharmacology analysis scripts."
echo ""
echo "Next steps:"
echo "1. Navigate back to project root: cd ../.."
echo "2. Run the analysis: Rscript scripts/R/01_complete_data_loading.R"
echo ""
echo "Note: These files are compressed and will be automatically"
echo "      decompressed by R when reading the data."
