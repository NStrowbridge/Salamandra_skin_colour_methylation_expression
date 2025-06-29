#!/bin/bash

# Define paths
BASE_DIR="/export/home4/2694872s/DirectRNA/NS9" #Directory with fast5 raw data file folder from nanopore

# activate conda environment with pod5
conda activate pod5_conversion

# Move to directory with yaml files

cd $BASE_DIR

pod5 convert fast5 NS9_fast5/*.fast5 -t 16 --output NS9_converted.pod5
