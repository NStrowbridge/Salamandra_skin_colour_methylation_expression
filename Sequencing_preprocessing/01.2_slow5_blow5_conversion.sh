#!/bin/bash

# Activate slow5tools Conda environment
conda activate slow5tools  # Ensure this environment is documented

# Define sample-specific variables
base_dir="/Volumes/LaCie/Nic_DirectRNA_raw_data/2024_DirectRNA_data"
sample_folder="Sample_folder"              # <-- Replace with actual sample folder name
sample_name="Sample_name"                  # <-- Replace with actual sample name
subdir="Subdirectory_with_fast5_folder"    # <-- Replace with actual subdirectory name

# Construct full path to the directory containing .fast5 files
full_path="${base_dir}/${sample_folder}/${sample_name}/${subdir}"
output_dir="${full_path}/blow5_output"
output_file="${sample_name}.blow5"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Navigate to the directory and run the conversion
cd "$full_path" || { echo "Directory not found: $full_path"; exit 1; }

# Convert a directory of FAST5 files into BLOW5 format (default compression: zlib+svb-zd)
slow5tools f2s fast5_dir -d "$output_dir"

echo "Conversion complete for sample: $sample_name"
