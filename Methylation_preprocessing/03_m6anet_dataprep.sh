#!/bin/bash
# Define variables
CONDA_ENV="f5c_eventalign"
SAMPLE_DIR="/export/home4/2694872s/DirectRNA/NS2_dor"
EVENTALIGN_OUTPUT="sample_eventalign.txt"
m6A_OUT_DIR="m6A_dataprep"

# Activate conda environment
conda activate $CONDA_ENV

# Move to directory of sample
cd $SAMPLE_DIR

# Run xpore dataprep
m6anet dataprep \
--eventalign $EVENTALIGN_OUTPUT \
--n_processes 48 \
--readcount_max 1000000 \
--out_dir $m6A_OUT_DIR
