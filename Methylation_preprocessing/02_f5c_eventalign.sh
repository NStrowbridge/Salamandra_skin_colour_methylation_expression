#!/bin/bash

# Define variables
CONDA_ENV="f5c_eventalign"
SAMPLE_DIR="Sample/Directory"
SLOW5_FILE="Sample_blow5.blow5"
FASTQ_FILE="Sample.fastq"
BAM_FILE="Sample_alignment.bam"
REFERENCE_TRANSCRIPTOME="/export/home4/2694872s/DirectRNA/NS9/pychopper_NS9_dorado_qscore9_output/NS9_dorado_qscore9_trimmed_reads_kraken_conf0.9_rRNA/RNAbloom_assembly_qscore9_pychopper_kraken2_rRNAremoved_k19/rnabloom_transcriptome_kmer19.fasta"
KMER_MODEL="rna004.nucleotide.5mer.model"
EVENTALIGN_OUTPUT="Sample_eventalign.txt"
XPORE_OUT_DIR="xpore_dataprep"

# Activate conda environment
conda activate $CONDA_ENV

# Move to directory of sample
cd $SAMPLE_DIR

# Index with FASTQ reads
f5c index --slow5 $SLOW5_FILE $FASTQ_FILE

# Eventalign and write to the named pipe
f5c eventalign --threads 48 --min-mapq 0 --rna \
    -b $BAM_FILE \
    -g $REFERENCE_TRANSCRIPTOME \
    -r $FASTQ_FILE \
    -o $EVENTALIGN_OUTPUT \
    --slow5 $SLOW5_FILE \
    --signal-index --scale-events
