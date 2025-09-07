#!/bin/bash

# Define k-mer sizes to test, 25 is default, so 6 either side and then 10 above
kmer_sizes=(19 25 31 41 51 61)

# Define input files and directories
base_dir="/export/home4/2694872s/DirectRNA/NS9/pychopper_NS9_dorado_qscore9_output/NS9_dorado_qscore9_trimmed_reads_kraken_conf0.9_rRNA"
pychopper_trimmed="${base_dir}/NS9_dorado_qscore9_kraken+rRNAfiltered.fastq" #File was renamed from NON_RRNA_READS of previous script
output_base_dir="RNAbloom_assembly_qscore9_pychopper_kraken2_rRNAremoved"
busco_lineage="tetrapoda_odb10" # Change depending on lineage of interest

# Change to the specified directory
cd $base_dir

# Loop through each k-mer size
for k in "${kmer_sizes[@]}"; do
    # Create output directory for each k-mer size
    output_dir="${output_base_dir}_k${k}"
    mkdir -p $output_dir

    # Activate RNA-Bloom conda environment and run RNA-Bloom
    conda activate rnabloom
    rnabloom -long $pychopper_trimmed -k $k -t 48 -outdir $output_dir

    # Activate BUSCO conda environment and run BUSCO
    conda activate BUSCO.1
    busco -m transcriptome -c 48 -i ${output_dir}/rnabloom.transcripts.fa -o ${output_dir}_busco -l $busco_lineage -f

    # Activate TransRate conda environment and run TransRate
    conda activate /export/home4/2694872s/anaconda3/envs/transcriptome_QC
    transrate --assembly ${output_dir}/rnabloom.transcripts.fa --left $pychopper_trimmed --output ${output_dir}_transrate
done

echo "All k-mer size tests completed."
