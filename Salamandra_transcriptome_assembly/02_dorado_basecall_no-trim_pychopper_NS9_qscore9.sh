#!/bin/bash

# Activate Dorado basecalling environment
conda activate dorado_basecalling

# Define paths and filenames
base_dir="/export/home4/2694872s/DirectRNA/NS9"
output_dir="${base_dir}/20250305_pychopper_NS9_output"
mkdir -p "$output_dir"

input_pod5="${base_dir}/NS9_converted.pod5"
dorado_model="dna_r9.4.1_e8_hac@v3.3"
dorado_output="${base_dir}/dorado_0.7.2_NS9_qscore9.fastq"

# Move to base directory and download model
cd "$base_dir" || exit
dorado download --model "$dorado_model"

# Run Dorado basecaller
dorado basecaller "$dorado_model" -v --no-trim --min-qscore 9 --emit-fastq "$input_pod5" > "$dorado_output"

# Check if Dorado finished successfully
if [ $? -eq 0 ]; then
    echo "Dorado basecalling completed successfully."

    # Activate Pychopper environment
    conda activate pychopper

    # Define Pychopper output filenames
    cd "$output_dir" || exit
    pychopper_input="$dorado_output"
    pychopper_prefix="NS9_dorado_qscore9"

    pychopper -k PCS109 \
        -r "${pychopper_prefix}_report.pdf" \
        -u "${pychopper_prefix}_unclassified.fq" \
        -S "${pychopper_prefix}_stats_output.tsv" \
        -w "${pychopper_prefix}_rescued.fq" \
        "$pychopper_input" \
        "${pychopper_prefix}_full_length_output.fq"

    echo "Pychopper processing completed."

else
    echo "Dorado basecalling failed. Exiting."
    exit 1
fi
