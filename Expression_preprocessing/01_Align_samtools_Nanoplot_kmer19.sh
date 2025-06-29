#!/bin/bash

# Define directories and other variables
DIRECTORY="/export/home4/2694872s/DirectRNA"
REFERENCE="/export/home4/2694872s/DirectRNA/NS9/pychopper_NS9_dorado_qscore9_output/NS9_dorado_qscore9_trimmed_reads_kraken_conf0.9_rRNA/RNAbloom_assembly_qscore9_pychopper_kraken2_rRNAremoved_k19/rnabloom.transcripts.fa"
SAMPLES=("NS2_dor" "NS2_lat" "NS3_dor" "NS3_lat" "NS5_dor" "NS5_lat" "NS7_dor" "NS7_lat" "NS8_dor" "NS8_lat" "NS10_dor" "NS10_lat" "NS11_dor" "NS11_lat" "NS12_dor" "NS12_lat" "ELT_14079_dor" "ELT_14079_lat" "ELT_14081_dor" "ELT_14081_lat" "ELT_14082_dor" "ELT_14082_lat" "ELT_14084_dor" "ELT_14084_lat")

# Activate minimap2 environment
conda activate minimap2

# Align with minimap2 v. 2.24
for ana in "${SAMPLES[@]}"; do
  cd "${DIRECTORY}/${ana}"
  minimap2 -ax splice -t 48 -uf -k14 "${REFERENCE}" "dorado_0.7.2_${ana}.fastq" > "${ana}_kmer19_transcriptome_aln.sam"
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
