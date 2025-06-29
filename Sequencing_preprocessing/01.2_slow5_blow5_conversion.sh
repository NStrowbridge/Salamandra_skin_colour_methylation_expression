#!/bin/bash

# Activate slow5tools environment
conda activate slow5tools

# Define sample-specific variables
base_dir="/Volumes/LaCie/Nic_DirectRNA_raw_data/2024_DirectRNA_data"
sample_folder="Sample_folder" #Change this depending on sample
sample_name="Sample_name"
subdir="Subdirectory_with_fast5_folder"

# Construct full path
full_path="${base_dir}/${sample_folder}/${sample_name}/${subdir}"
output_file="${sample_name}.blow5"

# Move to sample directory and run conversion
cd "$full_path" || { echo "Directory not found: $full_path"; exit 1; }

# convert a directory of fast5 files into BLOW5 files (default compression: zlib+svb-zd)
slow5tools f2s fast5_dir -d blow5_dir

echo "Conversion complete for $sample_name"
