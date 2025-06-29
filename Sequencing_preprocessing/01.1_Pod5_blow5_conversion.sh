#!/bin/bash

#Install blue-crab to convert pod5 to blow5
python3 -m venv ./blue-crab-venv
source ./blue-crab-venv/bin/activate
python3 -m pip install --upgrade pip

pip install blue-crab

# Activate blue-crab virtual environment
source ./blue-crab-venv/bin/activate

# Define sample-specific variables
base_dir="/Volumes/LaCie/Nic_DirectRNA_raw_data/2024_DirectRNA_data"
sample_folder="Sample_folder" #Change this depending on sample
sample_name="Sample_name"
subdir="Subdirectory_with_pod5_folder"


# Construct full path
full_path="${base_dir}/${sample_folder}/${sample_name}/${subdir}"
output_file="${sample_name}.blow5"

# Move to sample directory and run conversion
cd "$full_path" || { echo "Directory not found: $full_path"; exit 1; }
blue-crab p2s pod5 --output "$output_file"

echo "Conversion complete for $sample_name"
