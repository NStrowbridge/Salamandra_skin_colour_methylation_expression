The following provides a concise overview of each analysis script for differential methylation, including its purpose, dependencies, and software versions.

**01_get_functional_annotation_combined_RNAsamba_TransDecoder.R**

Description: A script that merges TransDecoder and RNAsamba eggNOG annotations, resolves conflicts by score, and outputs a unified functional transcriptome annotation set.
R version: 4.5.1
Packages: dplyr_1.2.0, readxl_1.4.5

**02.1_edgeR_methylation_transcript_level.R**

Description: A script which performs a full edgeR differential methylation analysis pipeline, from loading filtered transcript‑level methylation counts, normalising and batch‑correcting them, to fitting GLMs and generating differential methylation results between skin colours. For the manuscript the Yellow vs black skin and black vs brown skin comparisons were used.
R version: 4.5.1
Packages: sva_3.56.0, BiocParallel_1.42.2, genefilter_1.90.0, mgcv_1.9-4, nlme_3.1-168, locfit_1.5-9.12, svglite_2.2.2, edgeR_4.6.3, limma_3.64.3, tximport_1.36.1, GenomicFeatures_1.60.0, AnnotationDbi_1.70.0, Biobase_2.68.0, GenomicRanges_1.60.0, GenomeInfoDb_1.44.3, IRanges_2.42.0, S4Vectors_0.46.0, BiocGenerics_0.54.1, generics_0.1.4, BiocManager_1.30.27

**02.2_edgeR_methylation_transcript_level_morph.R**

Description: A script which performs a full edgeR differential methylation analysis pipeline, from loading filtered transcript‑level methylation counts, normalising and batch‑correcting them, to fitting GLMs and generating differential methylation results between colour morphs. For the manuscript the Yellow vs brown colour morph comparison was used.
R version: 4.5.1
Packages: sva_3.56.0, BiocParallel_1.42.2, genefilter_1.90.0, mgcv_1.9-4, nlme_3.1-168, locfit_1.5-9.12, svglite_2.2.2, edgeR_4.6.3, limma_3.64.3, tximport_1.36.1, GenomicFeatures_1.60.0, AnnotationDbi_1.70.0, Biobase_2.68.0, GenomicRanges_1.60.0, GenomeInfoDb_1.44.3, IRanges_2.42.0, S4Vectors_0.46.0, BiocGenerics_0.54.1, generics_0.1.4, BiocManager_1.30.27

**02.4_edgeR_methylation_transcript_level_stridorvstrilat.R**

Description: A script which performs a full edgeR differential methylation analysis pipeline, from loading filtered transcript‑level methylation counts, normalising and batch‑correcting them, to fitting GLMs and generating differential methylation results between striped black and striped yellow samples.
R version: 4.5.1
Packages: sva_3.56.0, BiocParallel_1.42.2, genefilter_1.90.0, mgcv_1.9-4, nlme_3.1-168, locfit_1.5-9.12, svglite_2.2.2, edgeR_4.6.3, limma_3.64.3, tximport_1.36.1, GenomicFeatures_1.60.0, AnnotationDbi_1.70.0, Biobase_2.68.0, GenomicRanges_1.60.0, GenomeInfoDb_1.44.3, IRanges_2.42.0, S4Vectors_0.46.0, BiocGenerics_0.54.1, generics_0.1.4, BiocManager_1.30.27

**03.1_dm_analysis_skin_colour_transcript_level_colour.R**

Description: A script which takes the edgeR differential‑methylation output and CPM matrix, merges them, identifies significant genes, and generates organised result tables and gene lists for each pairwise colour comparison. Again, for the manuscript the Yellow vs black skin and black vs brown skin comparisons were used.
R version: 4.5.1
Packages: reshape2_1.4.5, amap_0.8-20, ggrepel_0.9.7, ggplot2_4.0.2, rstudioapi_0.18.0


**03.2_dm_analysis_skin_colour_transcript_level_morph.R**

Description: A script which takes the edgeR differential‑methylation output and CPM matrix, merges them, identifies significant genes, and generates organised result tables and gene lists for each pairwise colour comparison. Again, for the manuscript the Yellow vs brown colour morph comparison was used.
R version: 4.5.1
Packages: reshape2_1.4.5, amap_0.8-20, ggrepel_0.9.7, ggplot2_4.0.2, rstudioapi_0.18.0


**04.1_dm_map_preferrednamextranscript_skin_colour.R**

Description: A script which walks through each differential‑methylation results folder and maps transcript IDs to their preferred gene names, producing matched gene‑name files for all gene lists. Again, for the manuscript the Yellow vs black skin and black vs brown skin comparisons were used.
R version: 4.5.1
Packages: dplyr_1.2.0


**04.2_dm_map_preferrednamextranscript_skin_morph.R**

Description: A script which walks through each differential‑methylation results folder and maps transcript IDs to their preferred gene names, producing matched gene‑name files for all gene lists. Again, for the manuscript the Yellow vs brown colour morph comparison was used.
R version: 4.5.1
Packages: dplyr_1.2.0


**05.1_BrovBla_methylation_get_transcriptID_compare_mapped_vs_nonmapped.py**

Description: A script that extracts FDR‑significant transcripts lacking a preferred‑name match in Brown vs Black colour skin comparison.
Python version: 3.11.8

**05.2_YelvBla_methylation_get_transcriptID_compare_mapped_vs_nonmapped.py**

Description: A script that extracts FDR‑significant transcripts lacking a preferred‑name match in Yellow vs Black colour skin comparison.
Python version: 3.11.8

**05.3_YelvBroM_methylation_get_transcriptID_compare_mapped_vs_nonmapped.py**

Description: A script that extracts FDR‑significant transcripts lacking a preferred‑name match in Yellow vs Brown morph comparison.
Python version: 3.11.8

**06_unique_transcripts_for_all_comparisons_DE+DM.py**

Description: A script which pulls in all FDR-significant transcripts that lack a preferred-name match from all methylation and expression comparisons, merges them, removes duplicates, and outputs one unified list of unique missing Transcript_ID–Preferred_name pairs across all analyses.
Python version: 3.11.8

**07_unique_transcripts_for_blast_split.py**

Description: A scrip that takes the full list of missing Transcript_ID–Preferred_name pairs and splits it into 13 evenly sized text files so they can be processed in smaller batches for blast.
Python version: 3.11.8

**08_get_transcript_sequence_loop.sh**

Description: A shell script which loops through each split missing Transcript_ID–Preferred_name pair lists, searches the transcriptome FASTA, and writes out the matching full sequences into separate .fa files.

**09_dm_graphs_analysis_colour_loci_gene_level.R**

Description: A script which creates the finalized plots for the published paper. Including volcano plots with pigment genes labelled, boxplots of pigment gene methylation, and combined Figure 2.
R version: 4.5.1
Packages: cowplot_1.2.0, reshape2_1.4.5, amap_0.8-20, ggrepel_0.9.7, gplot2_4.0.2, studioapi_0.18.0

**10_dm_+_dge_scatterplots_+_lm_+_fishertest.R**

Description: A script examines the relationship between differential methylation and differential expression using scatterplots, linear regression, and Fisher Exact tests.
R version: 4.5.1
Packages: purrr_1.2.1, patchwork_1.3.2, ggrepel_0.9.7, ggplot2_4.0.2, dplyr_1.2.0  
