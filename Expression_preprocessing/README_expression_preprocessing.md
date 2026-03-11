The following provides a concise overview of each analysis script for pre-processing the expression data, including its purpose and software versions.

**01_Align_count_samtools_Nanoplot_kmer19.sh**

Description: A script that aligns raw reads to de-novo transcriptome using minimap2, performs transcript counting per transcript using Salmon, using samtools to index and sort alignments, and NanoPlot to perform sequencing QC.
Packages: minimap2 2.28, salmon 1.10.1, samtools 1.18, NanoPlot 1.19.0
