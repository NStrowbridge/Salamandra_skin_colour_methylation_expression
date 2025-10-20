#!/bin/bash

# Define base directory and reference transcriptome
BASE_DIR="base_directory"
REFERENCE="${BASE_DIR}/reference_file"

# List of sample identifiers
SAMPLES=()

# Activate minimap2 environment and align reads
conda activate minimap2
for SAMPLE in "${SAMPLES[@]}"; do
  SAMPLE_DIR="${BASE_DIR}/${SAMPLE}"
  cd "$SAMPLE_DIR" || { echo "Directory not found: $SAMPLE_DIR"; continue; }

  minimap2 -ax splice -t 48 -uf -k14 "$REFERENCE" "dorado_0.7.2_${SAMPLE}.fastq" > "${SAMPLE}_kmer19_transcriptome_aln.sam"
done

# Activate samtools environment and process alignments
conda activate samtools
for SAMPLE in "${SAMPLES[@]}"; do
  SAMPLE_DIR="${BASE_DIR}/${SAMPLE}"
  cd "$SAMPLE_DIR" || { echo "Directory not found: $SAMPLE_DIR"; continue; }

  samtools flagstat --threads 48 "${SAMPLE}_kmer19_transcriptome_aln.sam" > "${SAMPLE}_flagstat_output_denovo.txt"
  samtools sort "${SAMPLE}_kmer19_transcriptome_aln.sam" -o "${SAMPLE}_kmer19_transcriptome_aln.bam"
  samtools index "${SAMPLE}_kmer19_transcriptome_aln.bam"
done

# Activate NanoPlot environment and generate plots
conda activate Nanoplot
for SAMPLE in "${SAMPLES[@]}"; do
  SAMPLE_DIR="${BASE_DIR}/${SAMPLE}"
  cd "$SAMPLE_DIR" || { echo "Directory not found: $SAMPLE_DIR"; continue; }
