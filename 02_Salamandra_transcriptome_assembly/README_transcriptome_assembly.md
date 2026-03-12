The following provides a concise overview of each analysis script for transcriptome assembly, including its purpose and software versions.

**01_NS9_fast5_to_pod5.sh**

Description: A script converts raw output files from Nanopore (fast5) to pod5 files for quick computation.
Packages: Pod5 0.3.23

**02_dorado_basecall_no-trim_pychopper_NS9_qscore9.sh**

Description: A script that basecalls the raw Nanopore files using no-trim, and then trims the sequences using pychopper.
Packages: Dorado 0.7.2, Pychopper 2

**03_Kraken2_pychopper_viral_bacterial_NS9_full_conf0.9.sh**

Description: A script that uses Kraken2 to remove contamination from reads using a database of commonly found bacteria and viruses.
Packages: Kraken 2.1.3

**04_minimap2_align_rRNA_NS9_kraken2_trimmed.sh**

Description: A script that removes rRNAs via alignment to the SortMeRNA standard database.
Packages: Minimap2 2.28

**05_Transcriptome_assembly+QC_iterate_kmer.sh**

Description: A script that performs an iterative transcriptome assembly using RNAbloom2, compares completeness of transcriptomes using BUSCO and finally compares assembly statistics using Transrate.
Packages: RNAbloom 2.0.1, BUSCO 6.0.0, transrate 1.0.3
