#!/bin/bash

# Activate Dorado basecalling environment
conda activate dorado_basecalling  # Ensure this environment is documented

# Define paths and filenames
project_directory="/export/home4/2694872s/DirectRNA/sample"
sample_name="sample"
input_pod5="${project_directory}/${sample_name}_converted.pod5"
dorado_model_name="RNA002_or_RNA004_model"
output_fastq="${project_directory}/dorado_0.7.2_${sample_name}.fastq"

# Move to project directory
cd "$project_directory" || { echo "Directory not found: $project_directory"; exit 1; }

# Download the Dorado model (if not already downloaded)
dorado download --model "$dorado_model_name"

# Run Dorado basecaller
dorado basecaller "$dorado_model_name" -v --emit-fastq "$input_pod5" > "$output_fastq"
