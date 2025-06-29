#!/bin/bash

# Set up and activate Blue Crab virtual environment
python3 -m venv ./blue-crab-venv
source ./blue-crab-venv/bin/activate
python3 -m pip install --upgrade pip
pip install blue-crab

# Define sample-specific variables
base_dir="/Volumes/LaCie/Nic_DirectRNA_raw_data/2024_DirectRNA_data"
sample_folder="Sample_folder"         # <-- Replace with actual sample folder name
sample_name="Sample_name"             # <-- Replace with actual sample name
subdir="Subdirectory_with_pod5_folder" # <-- Replace with actual subdirectory name

# Construct full path to the directory containing .pod5 files
full_path="${base_dir}/${sample_folder}/${sample_name}/${subdir}"
output_file="${sample_name}.blow5"

# Navigate to the directory and run the conversion
cd "$full_path" || { echo "Directory not found: $full_path"; exit 1; }

# Run Blue Crab conversion from POD5 to BLOW5
blue-crab p2s pod5 --output "$output_file"

echo "Conversion complete for sample: $sample_name"
