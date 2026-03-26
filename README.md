The project in this repository is from a paper published in Molecular Ecology (Strowbridge et al. 2026), which investigates how gene expression and m⁶A RNA methylation shape colour variation in the fire salamander, Salamandra salamandra bernardezi. Using long‑read PCR cDNA and direct RNA sequencing, the authors assembled a de novo transcriptome and compared black, yellow, and brown skin. They identified 129 differentially expressed genes and 281 differentially methylated genes, including many core pigmentation genes such as TYR, TYRP1, PMEL, and MLANA. Expression and methylation changes were positively correlated, and several pigmentation genes showed both effects, suggesting that epitranscriptomic regulation contributes to amphibian colour diversity. The work highlights RNA modifications as an important, previously overlooked mechanism for colour variation in non‑model vertebrates.

The following provides a concise overview of each folder, with a small description of what each contains.

**01_Sequencing_preprocessing**

Description: Scripts to convert and then basecall raw Nanopore signal.

**02_Salamandra_transcriptome_assembly**

Description: Scripts for basecalling Nanopore PCR cDNA sequencing, contamination removal, rRNA removal, and iterative transcriptome assembly and QC.

**03_Expression_preprocessing**

Description: Script for aligning direct RNA sequencing to de-novo transcriptome, producing transcript counts, and QC.

**04_Methylation_preprocessing**

Description: Scripts for creating eventalign files, generating m⁶A methylation site data, generating per transcript m⁶A counts, and combining those counts in a dataframe for differential methylation analysis.

**05_Differential_expression**

Description: Scripts for merging annotations across two different programs, differential expression analysis, processing transcripts which were not annotated in the transcriptome annotation, differential expression plot creation, and UpSet plot creation.

**06_Differential_methylation**

Description: Scripts for merging annotations across two different programs, differential methylation analysis, processing transcripts which were not annotated in the transcriptome annotation, extraction of sequence from non-annotated transcripts, differential methylation plot creation, differential methylation and differential expression scatterplot creation, linear models, and Fisher Exact Tests. 
