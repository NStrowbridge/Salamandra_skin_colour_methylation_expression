#!/bin/bash

# Define directories and other variables
DIRECTORY="/export/home4/2694872s/DirectRNA"
REFERENCE="rnabloom_transcriptome_kmer19.fasta"
SAMPLES=()

# Activate minimap2 environment
conda activate minimap2

# Align with minimap2 v. 2.24
for ana in "${SAMPLES[@]}"; do
  cd "${DIRECTORY}/${ana}"
  minimap2 -ax map-ont -t 48 "${REFERENCE}" "dorado_0.7.2_${ana}.fastq" > "${ana}_kmer19_transcriptome_aln.sam"
done

# Activate salmon environment

conda activate salmon

# Run salmon quant

for ana in "${array[@]}";do
  cd "${DIRECTORY}"
  mkdir salmon_quant
  cd salmon_quant
  salmon quant -t ${REFERENCE} -l A -a "${DIRECTORY}/${ana}/${ana}_kmer19_transcriptome_aln.sam" -o ${ana}_salmon_quant -p 48 --ont
done

# Activate samtools environment
conda activate samtools

# Run samtools flagstat and output to text file
for ana in "${SAMPLES[@]}"; do
  cd "${DIRECTORY}/${ana}"
  samtools flagstat --threads 48 "${ana}_kmer19_transcriptome_aln.sam" > "${ana}_flagstat_output_denovo.txt"
done

# Sort and index alignments
for ana in "${SAMPLES[@]}"; do
  cd "${DIRECTORY}/${ana}"
  samtools sort "${ana}_kmer19_transcriptome_aln.sam" > "${ana}_kmer19_transcriptome_aln.bam"
  samtools index "${ana}_kmer19_transcriptome_aln.bam"
done

# Activate NanoPlot environment
conda activate Nanoplot

# Produce NanoPlots
for ana in "${SAMPLES[@]}"; do
  cd "${DIRECTORY}/${ana}"
  NanoPlot -t 48 --color yellow --bam "${ana}_kmer19_transcriptome_aln.bam" -o "${ana}_denovo_kmer19_transcriptome_nanoplot"
done
