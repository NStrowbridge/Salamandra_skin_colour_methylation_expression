# Load required libraries
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

####~housekeeping~~~~~~~~~~~~####

rm(list=ls()) #clear the environment
setwd("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_methylation/01_scripts") #set wd

###~~output dir~~~~~~~~~~~~~~####
output = "../05.2_dm_graphs_analysis_colour_loci_gene_level/"
dir.create(output) #make folder for output
setwd(output)

# Load methylation data
meth_yelvbla <- read.csv("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_methylation/06_dm_analysis_skin_transcript_level/03_YelvBla/master.csv", row.names = 1)
meth_brovbla <- read.csv("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_methylation/06_dm_analysis_skin_transcript_level/01_BrovBla/master.csv", row.names = 1)
meth_yelvbro <- read.csv("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_methylation/06_dm_analysis_skin_transcript_level_morph/02_YelvBroM/master.csv", row.names = 1)

# Load expression data
expr_yelvbla <- read.csv("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/06_dge_analysis_skin_colour/03_YelvBla/master.csv", row.names = 1)
expr_brovbla <- read.csv("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/06_dge_analysis_skin_colour/01_BrovBla/master.csv", row.names = 1)
expr_yelvbro <- read.csv("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/06_dge_analysis_skin_morph/02_YelvBroM/master.csv", row.names = 1)


genes_yelvbla <- c("rb_77721", "rb_79752", "rb_62964", "rb_92274", "rb_78195", "rb_15163")
map_yelvbla <- c(
  "rb_77721" = "NOTCH1",
  "rb_79752" = "TYRP1",
  "rb_62964" = "TYR.I",
  "rb_92274" = "MLANA",
  "rb_78195" = "PMEL",
  "rb_15163" = "DGAT2"
)

genes_brovbla <- c("rb_1531", "rb_78661", "rb_79752", "rb_78195", "rb_91883", "rb_7890", "rb_77721", "rb_81552", "rb_15163")
map_brovbla <- c(
  "rb_1531" = "TYR.II",
  "rb_78661" = "GPNMB.I",
  "rb_79752" = "TYRP1",
  "rb_78195" = "PMEL",
  "rb_91883" = "ATOX1",
  "rb_7890" = "ATP7A",
  "rb_77721" = "NOTCH1",
  "rb_15163" = "DGAT2"
)

genes_yelvbro <- c("rb_78661", "rb_27044", "rb_91883", "rb_81618", "rb_76127", "rb_5385", "rb_79752", "rb_92274", "rb_78195")
map_yelvbro <- c(
  "rb_78661" = "GPNMB.I",
  "rb_27044" = "GPNMB.II",
  "rb_91883" = "ATOX1",
  "rb_81618" = "QDPR",
  "rb_76127" = "PAH",
  "rb_5385" = "TXN",
  "rb_79752" = "TYRP1",
  "rb_92274" = "MLANA",
  "rb_78195" = "PMEL"
)

# Define color palette
category_colors <- c(
  "Hyper Up Significant" = "#D73027",
  "Hyper Up Not Significant" = "#F46D43",
  "Hypo Down Significant" = "#4575B4",
  "Hypo Down Not Significant" = "#74ADD1",
  "Hyper Down Significant" = "#1A9850",
  "Hyper Down Not Significant" = "#A6D96A",
  "Hypo Up Significant" = "#762A83",
  "Hypo Up Not Significant" = "#C2A5CF",
  "Other Significant" = "#000000",
  "Other Not Significant" = "#BDBDBD"
)

# Helper function to prepare data
prepare_combined_data <- function(expr_df, meth_df, genes_to_label, gene_name_map, comparison_label) {
  merged <- inner_join(
    data.frame(gene = rownames(expr_df), expr_logFC = expr_df$logFC, expr_sig = expr_df$sigFDR),
    data.frame(gene = rownames(meth_df), meth_logFC = meth_df$logFC, meth_sig = meth_df$sigFDR),
    by = "gene"
  ) %>%
    mutate(
      category = case_when(
        expr_logFC > 1 & meth_logFC > 1 ~ "Hyper Up",
        expr_logFC < -1 & meth_logFC < -1 ~ "Hypo Down",
        expr_logFC < -1 & meth_logFC > 1 ~ "Hyper Down",
        expr_logFC > 1 & meth_logFC < -1 ~ "Hypo Up",
        TRUE ~ "Other"
      ),
      sig_status = ifelse(expr_sig & meth_sig, "Significant", "Not Significant"),
      category_sig = factor(paste(category, sig_status), levels = names(category_colors)),
      label = ifelse(gene %in% genes_to_label, gene_name_map[gene], NA),
      comparison = comparison_label
    )
  return(merged)
}

# Prepare all datasets
df1 <- prepare_combined_data(expr_yelvbla, meth_yelvbla, genes_yelvbla, map_yelvbla, "Yellow vs Black skin")
df2 <- prepare_combined_data(expr_brovbla, meth_brovbla, genes_brovbla, map_brovbla, "Brown vs Black skin")
df3 <- prepare_combined_data(expr_yelvbro, meth_yelvbro, genes_yelvbro, map_yelvbro, "Yellow vs Brown skin")

# Combine all into one data frame
combined_df <- bind_rows(df1, df2, df3)

# Set desired facet order
combined_df$comparison <- factor(combined_df$comparison, levels = c("Yellow vs Black skin", "Brown vs Black skin", "Yellow vs Brown skin"))

# Plot with facets
p <- ggplot(combined_df, aes(x = expr_logFC, y = meth_logFC, color = category_sig)) +
  geom_point(alpha = 0.7) +
  geom_label_repel(
    data = subset(combined_df, !is.na(label)),
    aes(label = label),
    fill = "white",
    size = 3,
    box.padding = 0.6,
    point.padding = 0.4,
    segment.color = "grey50",
    segment.size = 0.4,
    max.overlaps = Inf,
    force = 2,
    force_pull = 0.5,
    nudge_y = 0.2,
    show.legend = FALSE
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted", color = "grey40") +
  scale_color_manual(values = category_colors, drop = TRUE) +
  facet_wrap(~comparison) +
  labs(x = "Expression log2(Fold change)",
    y = "m6A methylation log2(Fold change)",
    color = "Classification"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  xlim(-10, 10) +
  ylim(-10, 10)

# Save the plot
ggsave("combined_facet_scatterplot.svg", p, width = 18, height = 6)