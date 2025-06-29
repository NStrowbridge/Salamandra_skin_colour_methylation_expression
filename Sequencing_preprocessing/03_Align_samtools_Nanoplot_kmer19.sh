#!/bin/bash

# Define base directory and reference transcriptome
BASE_DIR="/export/home4/2694872s/DirectRNA"
REFERENCE="${BASE_DIR}/NS9/pychopper_NS9_dorado_qscore9_output/NS9_dorado_qscore9_trimmed_reads_kraken_conf0.9_rRNA/RNAbloom_assembly_qscore9_pychopper_kraken2_rRNAremoved_k19/rnabloom.transcripts.fa"

# List of sample identifiers
SAMPLES=(
  "NS2_dor" "NS2_lat" "NS3_dor" "NS3_lat" "NS5_dor" "NS5_lat"
  "NS7_dor" "NS7_lat" "NS8_dor" "NS8_lat" "NS10_dor" "NS10_lat"
  "NS11_dor" "NS11_lat" "NS12_dor" "NS12_lat"
  "ELT_14079_dor" "ELT_14079_lat" "ELT_14081_dor" "ELT_14081_lat"
  "ELT_14082_dor" "ELT_14082_lat" "ELT_14084_dor" "ELT_14084_lat"
)

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
