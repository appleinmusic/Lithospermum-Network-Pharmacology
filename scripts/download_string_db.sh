#!/bin/bash
# STRING Database Download Script
# Downloads required STRING v12.0 database files for network analysis

echo "=== STRING Database Download Script ==="
echo "Downloading STRING v12.0 database files for Homo sapiens (9606)..."

# Create data directory if it doesn't exist
mkdir -p data/string_db
cd data/string_db

# Base URL for STRING database
BASE_URL="https://stringdb-downloads.org/download"

# Download protein information file
echo "1. Downloading protein information file..."
if [ ! -f "9606.protein.info.v12.0.txt.gz" ]; then
    curl -O "${BASE_URL}/protein.info.v12.0/9606.protein.info.v12.0.txt.gz"
    echo "✓ Downloaded: 9606.protein.info.v12.0.txt.gz"
else
    echo "✓ File already exists: 9606.protein.info.v12.0.txt.gz"
fi

# Download protein aliases file
echo "2. Downloading protein aliases file..."
if [ ! -f "9606.protein.aliases.v12.0.txt.gz" ]; then
    curl -O "${BASE_URL}/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz"
    echo "✓ Downloaded: 9606.protein.aliases.v12.0.txt.gz"
else
    echo "✓ File already exists: 9606.protein.aliases.v12.0.txt.gz"
fi

# Download protein links file (main interaction data)
echo "3. Downloading protein interaction links file..."
if [ ! -f "9606.protein.links.v12.0.txt.gz" ]; then
    curl -O "${BASE_URL}/protein.links.v12.0/9606.protein.links.v12.0.txt.gz"
    echo "✓ Downloaded: 9606.protein.links.v12.0.txt.gz"
else
    echo "✓ File already exists: 9606.protein.links.v12.0.txt.gz"
fi

# Optional: Download detailed links file (contains evidence scores)
echo "4. Downloading detailed protein links file (optional)..."
if [ ! -f "9606.protein.links.detailed.v12.0.txt.gz" ]; then
    curl -O "${BASE_URL}/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.txt.gz"
    echo "✓ Downloaded: 9606.protein.links.detailed.v12.0.txt.gz"
else
    echo "✓ File already exists: 9606.protein.links.detailed.v12.0.txt.gz"
fi

echo ""
echo "=== Download Summary ==="
echo "Files downloaded to: $(pwd)"
ls -lh *.gz

echo ""
echo "=== File Descriptions ==="
echo "• 9606.protein.info.v12.0.txt.gz      - Protein names and descriptions"
echo "• 9606.protein.aliases.v12.0.txt.gz   - Alternative protein names and IDs"
echo "• 9606.protein.links.v12.0.txt.gz     - Protein-protein interaction scores"
echo "• 9606.protein.links.detailed.v12.0.txt.gz - Detailed evidence scores (optional)"

echo ""
echo "=== Usage Instructions ==="
echo "These files will be automatically used by the R analysis scripts."
echo "Total approximate size: ~500MB"
echo "No need to extract - R will read .gz files directly"

echo ""
echo "✅ STRING database download completed successfully!"
echo "You can now run the network analysis scripts."