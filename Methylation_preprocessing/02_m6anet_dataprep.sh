#!/bin/bash

# Define variables
CONDA_ENV="f5c_eventalign"
SAMPLE_DIR="Sample/Directory"   # <-- Replace with actual sample directory
EVENTALIGN_OUTPUT="sample_eventalign.txt"              # <-- Replace with actual eventalign output file
M6A_OUT_DIR="m6A_dataprep"                             # <-- Output directory for m6Anet dataprep

# Activate the appropriate conda environment
conda activate "$CONDA_ENV"

# Navigate to the sample directory
cd "$SAMPLE_DIR" || { echo "Directory not found: $SAMPLE_DIR"; exit 1; }

# Run m6Anet dataprep
m6anet dataprep \
  --eventalign "$EVENTALIGN_OUTPUT" \
  --n_processes 48 \
  --readcount_max 1000000 \
  --out_dir "$M6A_OUT_DIR"

echo "m6Anet dataprep complete for: $EVENTALIGN_OUTPUT"
