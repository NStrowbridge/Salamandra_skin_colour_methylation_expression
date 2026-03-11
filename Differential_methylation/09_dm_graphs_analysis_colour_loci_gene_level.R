# Volcano plots and boxplots for differential expression analysis
# by Nic
# for analysis of DGE and DGM in yellow, brown and black skin

####~~~~load libraries~~~~####
library(rstudioapi)
library(ggplot2)
library(ggrepel)
#BiocManager::install("amap")
library(amap)
library(reshape2)
library(cowplot)

####~~~~housekeeping~~~~####

rm(list=ls()) #clear the environment
setwd("script_folder") #set wd

####~~~~output dir~~~~####
output = "../05.2_dm_graphs_analysis_colour_loci_gene_level/"
dir.create(output) #make folder for output
setwd(output)

####~~~~make combined differential methylation volcano plots~~~~####

# Function to create volcano plot without legend
create_volcano_plot <- function(data, genes_to_label, gene_name_mapping, title) {
  data$gene_label <- gene_name_mapping[rownames(data)]
  
  ggplot(data, aes(x = logFC, y = mlog10p, colour = direction)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50", size = 0.5) +
    scale_colour_manual(
      values = c("down" = "blue3", "ns" = "black", "up" = "red3"),
      labels = c("Hypomethylated", "Not Significant", "Hypermethylated"),
      name = NULL
    ) +
    geom_label_repel(
      data = subset(data, rownames(data) %in% genes_to_label),
      aes(label = gene_label),
      size = 5,
      box.padding = 0.5,
      point.padding = 0.5,
      segment.color = "gray50",
      segment.size = 0.4,
      arrow = arrow(length = unit(0.01, "npc"), type = "closed", ends = "last"),
      force = 2,
      max.overlaps = Inf,
      nudge_y = 0.5,
      show.legend = FALSE
    ) +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 P-value") +
    coord_cartesian(xlim = c(-10, 10), ylim = c(0, 7.5)) +
    theme_minimal(base_size = 20) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
}


# Load your data
yelvbla <- read.csv("/DirectRNA_Colour_bernardezi/Differential_methylation/06_dm_analysis_skin_transcript_level/03_YelvBla/master.csv", row.names = 1)
brovbla <- read.csv("/DirectRNA_Colour_bernardezi/Differential_methylation/06_dm_analysis_skin_transcript_level/01_BrovBla/master.csv", row.names = 1)
yelvbro <- read.csv("/DirectRNA_Colour_bernardezi/Differential_methylation/06_dm_analysis_skin_transcript_level_morph/02_YelvBroM/master.csv", row.names = 1)

# Define genes and mappings

genes_yelvbla <- c("rb_77721", "rb_79752")
map_yelvbla <- c("rb_77721" = "NOTCH1", "rb_79752" = "TYRP1")

genes_brovbla <- c("rb_1531", "rb_78661", "rb_79752", "rb_78195", "rb_91883", "rb_7890", "rb_7890", "rb_77721")
map_brovbla <- c("rb_1531" = "TYR.II", "rb_78661" = "GPNMB.I", "rb_79752" = "TYRP1", "rb_78195" = "PMEL", "rb_91883" = "ATOX1", "rb_7890" = "ATP7A", "rb_77721" = "NOTCH1")

genes_yelvbro <- c("rb_78661", "rb_27044", "rb_91883", "rb_81618")
map_yelvbro <- c("rb_78661" = "GPNMB.I", "rb_27044" = "GPNMB.II", "rb_91883" = "ATOX1", "rb_81618" = "QDPR" )

# Create plots
plot1 <- create_volcano_plot(yelvbla, genes_yelvbla, map_yelvbla, "Yellow vs Black skin")
plot2 <- create_volcano_plot(brovbla, genes_brovbla, map_brovbla, "Brown vs Black skin")
plot3 <- create_volcano_plot(yelvbro, genes_yelvbro, map_yelvbro, "Yellow vs Brown skin")

# Create a dummy plot to extract the legend
legend_plot <- ggplot(yelvbla, aes(x = logFC, y = mlog10FDR, colour = direction)) +
  geom_point() +
  scale_colour_manual(
    values = c("down" = "blue3", "ns" = "black", "up" = "red3"),
    labels = c("Hypomethylated", "Not Significant", "Hypermethylated"),
    name = NULL
  ) +
  theme_void() +
  theme(legend.position = "right", legend.text = element_text(size = 18))

# Extract the legend
legend <- get_legend(legend_plot)

combined <- plot_grid( plot1, plot2, plot3, ncol = 3, align = "hv" ) 

final_plot <- plot_grid( combined, legend, ncol = 2, rel_widths = c(3, 0.5) )

# Save
ggsave("combined_volcano_with_shared_legend.svg", final_plot, width = 20, height = 6)


####~~~~make combined differential methylation boxplots~~~~####
# Define color mapping
color_mapping <- c("Black" = "grey3", "Yellow" = "yellow3", "Brown" = "brown4")

# Function to create a boxplot for a given dataset
create_boxplot <- function(methylation_data, sample_sheet, gene_ids, gene_name_map, groups, title) {
  # Remove non-numeric columns
  methylation_data <- methylation_data[, sapply(methylation_data, is.numeric)]
  
  # Ensure sample names match between methylation data and sample sheet
  common_samples <- intersect(colnames(methylation_data), sample_sheet$SampleID)
  methylation_data <- methylation_data[, common_samples]
  sample_sheet <- sample_sheet[sample_sheet$SampleID %in% common_samples, ]
  
  # Reorder sample sheet to match methylation data
  sample_sheet <- sample_sheet[match(common_samples, sample_sheet$SampleID), ]
  
  # Log2 transform
  master_log2 <- log2(methylation_data + 1)
  
  # Extract and transpose candidate genes
  candidate_genes <- master_log2[gene_ids, , drop = FALSE]
  gene_data <- data.frame(t(candidate_genes))
  gene_data$sample_group <- sample_sheet$Condition
  
  # Filter for specified groups
  gene_data_filtered <- subset(gene_data, sample_group %in% groups)
  gene_data_filtered.m <- melt(gene_data_filtered, id.vars = "sample_group")
  gene_data_filtered.m$variable <- gene_name_map[gene_data_filtered.m$variable]
  
  # Create boxplot
  p <- ggplot(gene_data_filtered.m, aes(x = variable, y = value, fill = sample_group)) +
    geom_boxplot(outlier.size = 0) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 20)) +
    labs(title = title, y = "log2(Methylation + 1)", x = "Gene") +
    scale_fill_manual(values = color_mapping, drop = FALSE) +
    ylim(0,10)
  
  return(p)
}

# Load sample sheet

ss <- read.csv("/DirectRNA_Colour_bernardezi/Differential_methylation/02_reference_data/sample_sheet.csv")

# Load DM data

yelvbla <- read.csv("/DirectRNA_Colour_bernardezi/Differential_methylation/06_dm_analysis_skin_transcript_level/03_YelvBla/master.csv", row.names = 1)
brovbla <- read.csv("/DirectRNA_Colour_bernardezi/Differential_methylation/06_dm_analysis_skin_transcript_level/01_BrovBla/master.csv", row.names = 1)
yelvbro <- read.csv("/DirectRNA_Colour_bernardezi/Differential_methylation/06_dm_analysis_skin_transcript_level_morph/02_YelvBroM/master.csv", row.names = 1)

# Define genes and mappings

genes_yelvbla <- c("rb_77721", "rb_79752")
map_yelvbla <- c("rb_77721" = "NOTCH1", "rb_79752" = "TYRP1")

genes_brovbla <- c("rb_1531", "rb_78661", "rb_79752", "rb_78195", "rb_91883", "rb_7890", "rb_77721")
map_brovbla <- c("rb_1531" = "TYR.II", "rb_78661" = "GPNMB.I", "rb_79752" = "TYRP1", "rb_78195" = "PMEL", "rb_91883" = "ATOX1", "rb_7890" = "ATP7A", "rb_77721" = "NOTCH1")

genes_yelvbro <- c("rb_78661", "rb_27044", "rb_91883", "rb_81618")
map_yelvbro <- c("rb_78661" = "GPNMB.I", "rb_27044" = "GPNMB.II", "rb_91883" = "ATOX1", "rb_81618" = "QDPR" )



# Generate plots
plot1 <- create_boxplot(yelvbla, ss, genes_yelvbla, map_yelvbla, c("Yellow", "Black"),"")
plot2 <- create_boxplot(brovbla, ss, genes_brovbla, map_brovbla, c("Brown", "Black"),"")
plot3 <- create_boxplot(yelvbro, ss, genes_yelvbro, map_yelvbro, c("Yellow", "Brown"),"")

# Ensure consistent factor levels across all plots
all_levels <- c("Black", "Brown", "Yellow")
plot1 <- plot1 + scale_fill_manual(values = color_mapping, drop = FALSE, limits = all_levels) + theme(axis.title.x = element_blank())
plot2 <- plot2 + scale_fill_manual(values = color_mapping, drop = FALSE, limits = all_levels) + theme(axis.title.y = element_blank())
plot3 <- plot3 + scale_fill_manual(values = color_mapping, drop = FALSE, limits = all_levels) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

# Create a dummy plot to extract the legend

dummy_plot <- ggplot(
  data.frame(skin_colour = factor(c("Black", "Brown", "Yellow"), levels = c("Black", "Brown", "Yellow"))),
  aes(x = skin_colour, fill = skin_colour)
) +
  geom_bar() +
  scale_fill_manual(
    values = c("Black" = "grey3", "Brown" = "brown4", "Yellow" = "yellow3"),
    drop = FALSE,
    name = "Skin colour"
  ) +
  theme_void() +
  theme(legend.position = "right", legend.text = element_text(size = 16), legend.title = element_text(size = 18))


# Extract the legend
legend <- get_legend(dummy_plot)

combined <- plot_grid( plot1, plot2, plot3, ncol = 3, align = "hv" ) 

final_plot <- plot_grid( combined, legend, ncol = 2, rel_widths = c(3, 0.5) )

# Save the final plot
ggsave("combined_boxplot_with_single_legend.svg", final_plot, width = 20, height = 6)

####~~~~make combined differential methylation volcano and boxplot together~~~~####
# Rename volcano plots
volcano1 <- create_volcano_plot(yelvbla, genes_yelvbla, map_yelvbla, "Yellow vs Black skin")
volcano2 <- create_volcano_plot(brovbla, genes_brovbla, map_brovbla, "Brown vs Black skin")
volcano3 <- create_volcano_plot(yelvbro, genes_yelvbro, map_yelvbro, "Yellow vs Brown skin")
volcano1 <- volcano1 + ylab("-log10(pvalue)") + xlab("")
volcano2 <- volcano2 + ylab("") + xlab("log2(Fold Change)")
volcano3 <- volcano3 + ylab("") + xlab("")

# Rename boxplots
box1 <- create_boxplot(yelvbla, ss, genes_yelvbla, map_yelvbla, c("Yellow", "Black"), "")
box2 <- create_boxplot(brovbla, ss, genes_brovbla, map_brovbla, c("Brown", "Black"), "")
box3 <- create_boxplot(yelvbro, ss, genes_yelvbro, map_yelvbro, c("Yellow", "Brown"), "")
box1 <- box1 + theme(legend.position = "none")
box2 <- box2 + theme(legend.position = "none")
box3 <- box3 + theme(legend.position = "none")
box1 <- box1 + ylab("log2(m6A methylation + 1)") + xlab("")
box2 <- box2 + ylab("") + xlab("Gene")
box3 <- box3 + ylab("") + xlab("")

# Arrange volcano plots in a row
volcano_row <- plot_grid(volcano1, volcano2, volcano3, ncol = 3, align = "v")

# Arrange boxplots in a row
boxplot_row <- plot_grid(box1, box2, box3, ncol = 3, align = "v")

# Stack volcano and boxplot rows
right_column <- plot_grid(volcano_row, boxplot_row, ncol = 1, rel_heights = c(1, 1))

# Extract legends
volcano_legend <- get_legend(legend_plot)
boxplot_legend <- get_legend(dummy_plot)

# Combine legends
combined_legend <- plot_grid(volcano_legend, boxplot_legend, ncol = 1)

# Final layout with wider left column for Venn
final_layout <- plot_grid(
  plot_grid(right_column, combined_legend, ncol = 2, rel_widths = c(2.5, 0.4)),
  ncol = 1,
  rel_widths = c(1.4, 3) 
)

# Add labels directly to the ggdraw canvas
final_labeled <- ggdraw(final_layout) +
  draw_plot_label(
    label = c("a)", "b)"),
    x = c(0.02, 0.02),
    y = c(0.98, 0.5),
    fontface = "bold",
    size = 18
  )

# Save the final figure
ggsave("final_combined_figure_with_legends_labeled.svg", final_labeled, width = 17, height = 10)

####~~~~Supplemental boxplots for all pigment genes~~~~####

# Define gene set and mapping
all_genes <- unique(c(genes_yelvbla, genes_brovbla, genes_yelvbro))
all_gene_map <- c(map_yelvbla, map_brovbla, map_yelvbro)

color_mapping <- c("Black" = "grey3", "Brown" = "brown4", "Yellow" = "yellow3")

#Filter Yellow samples from yelvbro using Morph column
yellow_samples <- ss[ss$Morph == "Yellow", "SampleID"]
yellow_methylation <- yelvbro[, yellow_samples]
yellow_methylation <- yellow_methylation[, sapply(yellow_methylation, is.numeric)]

# Filter Black and Brown samples from brovbla using Condition column
bb_samples <- ss[ss$Condition %in% c("Black", "Brown"), "SampleID"]
bb_samples <- intersect(bb_samples, colnames(brovbla))
bb_methylation <- brovbla[, bb_samples]
bb_methylation <- bb_methylation[, sapply(bb_methylation, is.numeric)]

# Align gene sets and combine
bb_methylation_subset <- bb_methylation[all_genes, , drop = FALSE]
yellow_methylation_subset <- yellow_methylation[all_genes, , drop = FALSE]

# Order rows for consistency
bb_methylation_subset <- bb_methylation_subset[order(rownames(bb_methylation_subset)), ]
yellow_methylation_subset <- yellow_methylation_subset[order(rownames(yellow_methylation_subset)), ]

# Ensure both datasets have the same row names
common_rows <- intersect(rownames(bb_methylation_subset), rownames(yellow_methylation_subset))

# Subset both datasets to only include matching rows
bb_matched <- bb_methylation_subset[common_rows, ]
yellow_matched <- yellow_methylation_subset[common_rows, ]

# Combine the matched datasets column-wise
combined_methylation <- cbind(bb_matched, yellow_matched)

# Save the final list of filtered genes
filtered_genes <- common_rows

# Match and reorder sample sheet
combined_sample_sheet <- ss[ss$SampleID %in% colnames(combined_methylation), ]
combined_sample_sheet <- combined_sample_sheet[match(colnames(combined_methylation), combined_sample_sheet$SampleID), ]

# Log2 transform and reformat
combined_log2 <- log2(combined_methylation + 1)
candidate_genes <- combined_log2[filtered_genes, , drop = FALSE]
gene_data <- data.frame(t(candidate_genes))
gene_data$sample_group <- combined_sample_sheet$Condition

# Create filtered gene map using the filtered gene list
filtered_gene_map <- all_gene_map[filtered_genes]

gene_data_m <- melt(gene_data, id.vars = "sample_group")
gene_data_m$variable <- filtered_gene_map[gene_data_m$variable]


# Plot boxplots
boxplot_all <- ggplot(gene_data_m, aes(x = variable, y = value, fill = sample_group)) +
  geom_boxplot(outlier.size = 0) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 16)) +
  labs(y = "log2(methylation + 1)", x = "Gene") +
  scale_fill_manual(values = color_mapping, drop = FALSE)

# Save plot

# Add legend to the combined_all_colours_boxplot
legend_dummy <- ggplot(
  data.frame(skin_colour = factor(c("Black", "Brown", "Yellow"), levels = c("Black", "Brown", "Yellow"))),
  aes(x = skin_colour, fill = skin_colour)
) +
  geom_bar() +
  scale_fill_manual(
    values = c("Black" = "grey3", "Brown" = "brown4", "Yellow" = "yellow3"),
    drop = FALSE,
    name = "Skin colour"
  ) +
  theme_void() +
  theme(legend.position = "right")

legend_all <- get_legend(legend_dummy)
boxplot_all_with_legend <- plot_grid(boxplot_all + theme(legend.position = "none"), legend_all, ncol = 2, rel_widths = c(3.5, 0.5))
ggsave("combined_all_colours_boxplot.svg", boxplot_all_with_legend, width = 13, height = 6)
