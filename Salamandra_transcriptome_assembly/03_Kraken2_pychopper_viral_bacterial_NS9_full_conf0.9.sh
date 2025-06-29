#!/bin/bash

# Activate Kraken2 conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate kraken2

# Define paths and filenames
base_dir="/export/home4/2694872s/DirectRNA/NS9/pychopper_NS9_dorado_qscore9_output"
input_file="NS9_dorado_qscore9_pychopper_full.fq"
kraken_db="/export/home4/2694872s/DirectRNA/kraken2_db/kraken2_db"
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
