#!/bin/bash

# Activate blue-crab environment
source ./blue-crab-venv/bin/activate

# Define base directory for RNA004 DirectRNA data
BASE_DIR="Base_dir"

# Define sample identifiers
SAMPLES=()

# Run blue-crab conversion
for ana in "${SAMPLES[@]}"; do
  cd "${BASE_DIR}/${ana}"
  blue-crab p2s pod5 --output "${ana}.blow5"
done

# Activate slow5tools environment
conda create -n slow5tools -y
conda activate slow5tools
conda install -y slow5tools ont_vbz_hdf_plugin

# Define base directory for RNA002 DirectRNA data
BASE_DIR="Base_directory"
SAMPLES=()

# Run slow5tools conversion and merge
for ana in "${SAMPLES[@]}"; do
  cd "${BASE_DIR}/${ana}_raw_data"
  slow5tools f2s "${ana}_fast5" -d "${ana}_blow5"
  slow5tools merge "${ana}_blow5" -o "${ana}.blow5"
done
