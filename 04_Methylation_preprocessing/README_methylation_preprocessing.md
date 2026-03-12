The following provides a concise overview of each analysis script for pre-processing the methylation data, including its purpose and software versions.

**01_fast5_pod5_to_blow5_conversion.sh**

Description: A script converts raw output files from Nanopore (fast5) to blow5 for creation of eventalign files.
Packages: blue-crab 0.2.0, slow5tools 1.3.0

**02_f5c_eventalign.sh**

Description: A script that creates eventalign files to detect methylation marks in Nanopore directRNA sequencing.
Packages: F5C 1.5

**03_m6anet_dataprep.sh**

Description: A script that uses m6anet to identify m⁶A RNA methylation modifications on each nucleotide.
Packages: m6anet 2.1.0

**04_m6anet_methylated_reads_per_transcript.sh**

Description: A script that uses awk to get number of methylated reads per transcript.

**04.1_m6anet_combine_m6a_reads.R**

Description: A script that imports all m6Anet output files, merges them by transcript, fills missing values, and produces a single matrix of per‑sample m⁶A‑modified read counts.
R version: 4.5.1
Packages: tibble_3.3.1, purrr_1.2.1, dplyr_1.2.0, readr_2.2.0
