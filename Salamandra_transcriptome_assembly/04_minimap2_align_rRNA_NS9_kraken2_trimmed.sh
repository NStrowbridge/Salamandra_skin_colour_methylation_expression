#!/bin/bash

###Activate minimap2  environment
conda activate minimap2

# Define paths
BASE_DIR="base_directory_with_fastq"
READS="input_reads_trimmed_from_kraken2"
RRNA_DB="rRNA_db_directory/rRNA_fasta_file"
OUTPUT_DIR="output_directory"
NON_RRNA_READS="${OUTPUT_DIR}/non_rrna_reads.fastq"
RRNA_READS="${OUTPUT_DIR}/rrna_reads.sam"

# Create output directory if it doesn't exist
cd $BASE_DIR
mkdir -p $OUTPUT_DIR

# Align reads to rRNA database using Minimap2
minimap2 -ax map-ont $RRNA_DB $READS > $RRNA_READS

# Extract non-rRNA reads
samtools view -b -f 4 $RRNA_READS | samtools fastq - > $NON_RRNA_READS

echo "Non-rRNA reads have been saved to $NON_RRNA_READS"
