#!/bin/bash

# Define paths
BASE_DIR="base_dir" #Change depending on sample

# activate conda environment with pod5
conda activate pod5_conversion

# Move to directory with yaml files

cd $BASE_DIR

pod5 convert fast5 Sample_fast5/*.fast5 -t 16 --output Sample_converted.pod5
