#!/bin/bash

###Activate minimap2  environment
conda activate minimap2

# Define paths
BASE_DIR="/export/home4/2694872s/DirectRNA/NS9/pychopper_NS9_dorado_qscore9_output"
READS="NS9_dorado_qscore9_trimmed_reads_kraken_conf0.9.fastq"
RRNA_DB="/export/home4/2694872s/DirectRNA/rRNA_databases_v4/smr_v4.3_default_db.fasta"
OUTPUT_DIR="NS9_dorado_qscore9_trimmed_reads_kraken_conf0.9_rRNA"
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
