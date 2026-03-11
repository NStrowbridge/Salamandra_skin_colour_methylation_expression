The following provides a concise overview of each analysis script for differential expression, including its purpose, dependencies, and software versions.

**01_get_functional_annotation_combined_RNAsamba_TransDecoder.R**

Description: A script that merges TransDecoder and RNAsamba eggNOG annotations, resolves conflicts by score, and outputs a unified functional transcriptome annotation set.
R version: 4.5.1
Packages: dplyr_1.2.0, readxl_1.4.5

**02.1_edgeR_methylation_transcript_level.R**

Description: A script which performs a full edgeR differential expression analysis pipeline, from loading filtered transcript‑level expression counts, normalising and batch‑correcting them, to fitting GLMs and generating differential expression results between skin colours. For the manuscript the Yellow vs black skin and black vs brown skin comparisons were used.
R version: 4.5.1
Packages: sva_3.56.0, BiocParallel_1.42.2, genefilter_1.90.0, mgcv_1.9-4, nlme_3.1-168, locfit_1.5-9.12, svglite_2.2.2, edgeR_4.6.3, limma_3.64.3, tximport_1.36.1, GenomicFeatures_1.60.0, AnnotationDbi_1.70.0, Biobase_2.68.0, GenomicRanges_1.60.0, GenomeInfoDb_1.44.3, IRanges_2.42.0, S4Vectors_0.46.0, BiocGenerics_0.54.1, generics_0.1.4, BiocManager_1.30.27

**02.2_edgeR_methylation_transcript_level_morph.R**

Description: A script which performs a full edgeR differential expression analysis pipeline, from loading filtered transcript‑level expression counts, normalising and batch‑correcting them, to fitting GLMs and generating differential expression results between colour morphs. For the manuscript the Yellow vs brown colour morph comparison was used.
R version: 4.5.1
Packages: sva_3.56.0, BiocParallel_1.42.2, genefilter_1.90.0, mgcv_1.9-4, nlme_3.1-168, locfit_1.5-9.12, svglite_2.2.2, edgeR_4.6.3, limma_3.64.3, tximport_1.36.1, GenomicFeatures_1.60.0, AnnotationDbi_1.70.0, Biobase_2.68.0, GenomicRanges_1.60.0, GenomeInfoDb_1.44.3, IRanges_2.42.0, S4Vectors_0.46.0, BiocGenerics_0.54.1, generics_0.1.4, BiocManager_1.30.27

**02.3_edgeR_methylation_transcript_level_stridorvstrilat.R**

Description: A script which performs a full edgeR differential methylation analysis pipeline, from loading filtered transcript‑level methylation counts, normalising and batch‑correcting them, to fitting GLMs and generating differential methylation results between striped black and striped yellow samples.
R version: 4.5.1
Packages: sva_3.56.0, BiocParallel_1.42.2, genefilter_1.90.0, mgcv_1.9-4, nlme_3.1-168, locfit_1.5-9.12, svglite_2.2.2, edgeR_4.6.3, limma_3.64.3, tximport_1.36.1, GenomicFeatures_1.60.0, AnnotationDbi_1.70.0, Biobase_2.68.0, GenomicRanges_1.60.0, GenomeInfoDb_1.44.3, IRanges_2.42.0, S4Vectors_0.46.0, BiocGenerics_0.54.1, generics_0.1.4, BiocManager_1.30.27

**03.1_dge_analysis_skin_colour.R**

Description: A script which takes the edgeR differential‑methylation output and CPM matrix, merges them, identifies significant genes, and generates organised result tables and gene lists for each pairwise colour comparison. Again, for the manuscript the Yellow vs black skin and black vs brown skin comparisons were used.
R version: 4.5.1
Packages: reshape2_1.4.5, amap_0.8-20, ggrepel_0.9.7, ggplot2_4.0.2, rstudioapi_0.18.0

**03.2_dge_analysis_skin_colour_morph.R**

Description: A script which takes the edgeR differential‑methylation output and CPM matrix, merges them, identifies significant genes, and generates organised result tables and gene lists for each pairwise colour comparison. Again, for the manuscript the Yellow vs brown colour morph comparison was used.
R version: 4.5.1
Packages: reshape2_1.4.5, amap_0.8-20, ggrepel_0.9.7, ggplot2_4.0.2, rstudioapi_0.18.0

**04.1_dge_map_preferrednamextranscript_skin_colour.R**

Description: A script which walks through each differential‑methylation results folder and maps transcript IDs to their preferred gene names, producing matched gene‑name files for all gene lists. Again, for the manuscript the Yellow vs black skin and black vs brown skin comparisons were used.
R version: 4.5.1
Packages: dplyr_1.2.0

**04.2_dge_map_preferrednamextranscript_skin_morph.R**

Description: A script which walks through each differential‑methylation results folder and maps transcript IDs to their preferred gene names, producing matched gene‑name files for all gene lists. Again, for the manuscript the Yellow vs brown colour morph comparison was used.
R version: 4.5.1
Packages: dplyr_1.2.0

**05.1_BrovBla_expression_get_transcriptID_compare_mapped_vs_nonmapped.py**

Description: A script that extracts FDR‑significant transcripts lacking a preferred‑name match in Brown vs Black colour skin comparison.
Python version: 3.11.8

**05.2_YelvBla_expression_get_transcriptID_compare_mapped_vs_nonmapped.py**

Description: A script that extracts FDR‑significant transcripts lacking a preferred‑name match in Yellow vs Black colour skin comparison.
Python version: 3.11.8

**05.3_YelvBroM_expression_get_transcriptID_compare_mapped_vs_nonmapped.py**

Description: A script that extracts FDR‑significant transcripts lacking a preferred‑name match in Yellow vs Brown morph comparison.
Python version: 3.11.8

**06_dge_graphs_analysis_colour_loci_gene_level.R**

Description: A script which creates the finalized plots for the published paper. Including volcano plots with pigment genes labelled, boxplots of pigment gene methylation, and combined Figure 2.
R version: 4.5.1
Packages: cowplot_1.2.0, reshape2_1.4.5, amap_0.8-20, ggrepel_0.9.7, gplot2_4.0.2, studioapi_0.18.0

**07_UpSet_plot_graph.R**

Description: A script that compiles methylation and expression gene lists across all comparisons and visualises their overlap using a customised multi‑set Upset plot.
R version: 4.5.1
Packages: stringr_1.6.0, ggtext_0.1.2, ggplot2_4.0.2, ComplexUpset_1.3.3, dplyr_1.2.0    

**08_StridorvStrilat_graphs+annotation**

Description: A script that identifies extreme log‑fold‑change transcripts in striped yellow vs striped black, plots their DE/DM P‑value patterns, annotates them (including manual BLAST fixes), and outputs summarised gene‑level results.
R version: 4.5.1
Packages: dplyr_1.2.0, ggrepel_0.9.7, ggplot2_4.0.2
