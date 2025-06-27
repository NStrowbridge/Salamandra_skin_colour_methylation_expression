####~load libraries~~~~~~~~~~####
library(dplyr)
library(ComplexUpset)
library(ggplot2)
library(dplyr)
library(ggtext)
library(stringr)

####~housekeeping~~~~~~~~~~~~####

rm(list=ls()) #clear the environment
setwd("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/01_scripts") #set wd

###~~output dir~~~~~~~~~~~~~~####
output = "../05.2_dge_graphs_analysis_colour_loci_gene_level/"
dir.create(output) #make folder for output
setwd(output)

####~~~load data~~~####

# Load methylation data
meth_yelvbla <- read.csv("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_methylation/02_reference_data/YelvBla_DM_anno.csv")
colnames(meth_yelvbla) = meth_yelvbla[1, ]
meth_yelvbla <- meth_yelvbla[-c(1),]
meth_brovbla <- read.csv("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_methylation/02_reference_data/BrovBla_DM_anno.csv")
colnames(meth_brovbla) = meth_brovbla[1, ]
meth_brovbla <- meth_brovbla[-c(1),]
meth_yelvbro <- read.csv("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_methylation/02_reference_data/YelvBroM_DM_anno.csv")
colnames(meth_yelvbro) = meth_yelvbro[1, ]
meth_yelvbro <- meth_yelvbro[-c(1),]

# Load expression data
expr_yelvbla <- read.csv("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/02_reference_data/YelvBla_DE_anno.csv")
colnames(expr_yelvbla) = expr_yelvbla[1, ]
expr_yelvbla <- expr_yelvbla[-c(1),]
expr_brovbla <- read.csv("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/02_reference_data/BrovBla_DE_anno.csv")
colnames(expr_brovbla) = expr_brovbla[1, ]
expr_brovbla <- expr_brovbla[-c(1),]
expr_yelvbro <- read.csv("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/02_reference_data/YelvBroM_DE_anno.csv")
colnames(expr_yelvbro) = expr_yelvbro[1, ]
expr_yelvbro <- expr_yelvbro[-c(1),]

# Extract Putative Gene columns
genes_yelvbla_meth <- meth_yelvbla$`Putative Gene`
genes_yelvbla_expr <- expr_yelvbla$`Putative Gene`
genes_brovbla_meth <- meth_brovbla$`Putative Gene`
genes_brovbla_expr <- expr_brovbla$`Putative Gene`
genes_yelvbro_meth <- meth_yelvbro$`Putative Gene`
genes_yelvbro_expr <- expr_yelvbro$`Putative Gene`

# Combine all unique genes
all_genes <- unique(c(
  genes_yelvbla_meth, genes_yelvbla_expr,
  genes_brovbla_meth, genes_brovbla_expr,
  genes_yelvbro_meth, genes_yelvbro_expr
))

# Create presence/absence matrix
gene_matrix <- data.frame(Gene = all_genes) %>%
  mutate(
    "Yellow vs Black Methylation" = Gene %in% genes_yelvbla_meth,
    "Yellow vs Black Expression" = Gene %in% genes_yelvbla_expr,
    "Brown vs Black Methylation" = Gene %in% genes_brovbla_meth,
    "Brown vs Black Expression" = Gene %in% genes_brovbla_expr,
    "Yellow vs Brown Methylation" = Gene %in% genes_yelvbro_meth,
    "Yellow vs Brown Expression" = Gene %in% genes_yelvbro_expr
  )

# Create the plot
upset_plot_meth_expr <- upset(
  gene_matrix,
  intersect = c(
    "Yellow vs Black Methylation", "Yellow vs Black Expression",
    "Brown vs Black Methylation", "Brown vs Black Expression",
    "Yellow vs Brown Methylation", "Yellow vs Brown Expression"
  ),
  name = "Putative Gene Overlap",
  width_ratio = 0.25,
  base_annotations=list(
    'Intersection size'=intersection_size(counts=FALSE)
  )
    )

ggsave("upset_plot_meth_expr.svg", upset_plot_meth_expr, width = 14, height = 6)
