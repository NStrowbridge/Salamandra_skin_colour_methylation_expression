#!/bin/bash

# Activate Kraken2 conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate kraken2

# Define paths and filenames
base_dir="base_dir"
input_file="input_fastq"
kraken_db="Kraken2_db_directory"
threads=32
confidence=0.90

# Extract base name (without extension)
cd "$base_dir" || exit
basename="${input_file%.*}"

# Define output filenames
report_file="${basename}_kraken_conf${confidence}.txt"
unclassified_file="${basename}_kraken_conf${confidence}.fastq"

# Run Kraken2
kraken2 \
  --db "$kraken_db" \
  --threads "$threads" \
  --confidence "$confidence" \
  --report "$report_file" \
  --output "$report_file" \
  --unclassified-out "$unclassified_file" \
  "$input_file"

echo "Kraken2 classification complete for $input_file"
