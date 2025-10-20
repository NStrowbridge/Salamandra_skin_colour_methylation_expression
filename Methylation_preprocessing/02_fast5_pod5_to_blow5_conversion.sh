#!/bin/bash

# Activate blue-crab environment
source ./blue-crab-venv/bin/activate

# Define base directory for RNA004 DirectRNA data
BASE_DIR="Base_directory"

# Define sample identifiers
SAMPLES=("ELT14079_dor" "ELT14081_lat" "ELT14081_dor" "ELT14082_lat" "ELT14082_dor" "ELT14084_lat" "ELT14084_dor" "NS12_lat" "NS12_dor" "NS11_lat" "NS11_dor")

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
SAMPLES=("NS7_lat" "NS7_dor" "NS1_lat" "NS2_lat" "NS2_dor" "NS3_lat" "NS3_dor" "NS5_lat" "NS5_dor" "NS8_lat" "NS8_dor" "NS10_lat" "NS10_dor")

# Run slow5tools conversion and merge
for ana in "${SAMPLES[@]}"; do
  cd "${BASE_DIR}/${ana}_raw_data"
  slow5tools f2s "${ana}_fast5" -d "${ana}_blow5"
  slow5tools merge "${ana}_blow5" -o "${ana}.blow5"
done
