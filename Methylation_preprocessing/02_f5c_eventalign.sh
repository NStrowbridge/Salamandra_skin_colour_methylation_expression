#!/bin/bash

# Define variables
CONDA_ENV="f5c_eventalign"
SAMPLE_DIR="Sample/Directory"                     # <-- Replace with actual sample directory
SLOW5_FILE="Sample_blow5.blow5"                   # <-- Replace with actual BLOW5 file
FASTQ_FILE="Sample.fastq"                         # <-- Replace with actual FASTQ file
BAM_FILE="Sample_alignment.bam"                   # <-- Replace with aligned BAM file
REFERENCE_TRANSCRIPTOME="transcriptome_path"
KMER_MODEL="rna004.nucleotide.5mer.model"         # <-- Use only if using RNA004 chemistry
EVENTALIGN_OUTPUT="Sample_eventalign.txt"         # <-- Output file for eventalign results

# Activate the appropriate conda environment
conda activate "$CONDA_ENV"

# Navigate to the sample directory
cd "$SAMPLE_DIR" || { echo "Directory not found: $SAMPLE_DIR"; exit 1; }

# Index the FASTQ reads using the BLOW5 file
f5c index --slow5 "$SLOW5_FILE" "$FASTQ_FILE"

# Run f5c eventalign with signal scaling and signal index output
f5c eventalign \
    --threads 48 \
    --min-mapq 0 \
    --rna \
    -b "$BAM_FILE" \
    -g "$REFERENCE_TRANSCRIPTOME" \
    -r "$FASTQ_FILE" \
    -o "$EVENTALIGN_OUTPUT" \
    --slow5 "$SLOW5_FILE" \
    --kmer-model "$KMER_MODEL" \
    --signal-index \
    --scale-events

echo "Eventalign complete for sample: $FASTQ_FILE"
