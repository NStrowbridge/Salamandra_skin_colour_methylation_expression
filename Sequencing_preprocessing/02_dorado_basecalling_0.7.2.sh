#!/bin/bash

# Activate Dorado basecalling environment
conda activate dorado_basecalling

# Define paths and filenames
base_dir="/export/home4/2694872s/DirectRNA/sample"
input_pod5="${base_dir}/sample_converted.pod5"
dorado_model="rna002_70bps_hac@v3"
output_fastq="${base_dir}/dorado_0.7.2_sample.fastq"

# Move to base directory
cd "$base_dir" || exit

# Download the Dorado model
dorado download --model "$dorado_model"

# Run Dorado basecaller
dorado basecaller "$dorado_model" -v --emit-fastq "$input_pod5" > "$output_fastq"
