####~load libraries~~~~~~~~~~####
library(rstudioapi)
library(ggplot2)
library(ggrepel)
library(grid)
library(patchwork)
library(cowplot)
#BiocManager::install("amap")
library(amap)
library(reshape2)

####~housekeeping~~~~~~~~~~~~####

rm(list=ls()) #clear the environment
setwd("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/01_scripts") #set wd

###~~output dir~~~~~~~~~~~~~~####
output = "../05.2_dge_graphs_analysis_colour_loci_gene_level/"
dir.create(output) #make folder for output
setwd(output)

#~~~~make combined volcano plots~~~~~~~~~~~~###
# Function to create volcano plot without legend
create_volcano_plot <- function(data, genes_to_label, gene_name_mapping, title) {
  data$gene_label <- gene_name_mapping[rownames(data)]
  
  ggplot(data, aes(x = logFC, y = mlog10p, colour = direction)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50", size = 0.5) +
    scale_colour_manual(
      values = c("down" = "blue3", "ns" = "black", "up" = "red3"),
      labels = c("Downregulated", "Not Significant", "Upregulated"),
      name = NULL
    ) +
    geom_label_repel(
      data = subset(data, rownames(data) %in% genes_to_label),
      aes(label = gene_label),
      size = 3,
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
    coord_cartesian(xlim = c(-6, 6), ylim = c(0, 10)) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
  
}


# Load data from DGE analysis
yelvbla <- read.csv("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/06_dge_analysis_skin_colour/03_YelvBla/master.csv", row.names = 1)
brovbla <- read.csv("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/06_dge_analysis_skin_colour/01_BrovBla/master.csv", row.names = 1)
yelvbro <- read.csv("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/06_dge_analysis_skin_morph/02_YelvBroM/master.csv", row.names = 1)

# Define genes and mappings
genes_yelvbla <- c("rb_62964", "rb_79752", "rb_92274", "rb_78195", "rb_15163")
map_yelvbla <- c("rb_62964" = "TYR.I", "rb_79752" = "TYRP1", "rb_92274" = "MLANA", "rb_78195" = "PMEL", "rb_15163" = "DGAT2")

genes_brovbla <- c("rb_1531", "rb_78661", "rb_81552", "rb_15163")
map_brovbla <- c("rb_1531" = "TYR.II", "rb_78661" = "GPNMB.I", "rb_15163" = "DGAT2")

genes_yelvbro <- c("rb_76127", "rb_5385", "rb_78661", "rb_79752", "rb_92274", "rb_78195")
map_yelvbro <- c("rb_76127" = "PAH", "rb_5385" = "TXN", "rb_78661" = "GPNMB.I", "rb_79752" = "TYRP1", "rb_92274" = "MLANA", "rb_78195" = "PMEL")

# Create plots
plot1 <- create_volcano_plot(yelvbla, genes_yelvbla, map_yelvbla, "Yellow vs Black skin")
plot2 <- create_volcano_plot(brovbla, genes_brovbla, map_brovbla, "Brown vs Black skin")
plot3 <- create_volcano_plot(yelvbro, genes_yelvbro, map_yelvbro, "Yellow vs Brown skin")

# Create a dummy plot to extract the legend
legend_plot <- ggplot(yelvbla, aes(x = logFC, y = mlog10FDR, colour = direction)) +
  geom_point() +
  scale_colour_manual(
    values = c("down" = "blue3", "ns" = "black", "up" = "red3"),
    labels = c("Downregulated", "Not Significant", "Upregulated"),
    name = NULL
  ) +
  theme_void() +
  theme(legend.position = "right")

# Extract legend
legend <- get_legend(legend_plot)

# Combine plots with shared legend
combined_plot <- (plot1 + plot2 + plot3) + plot_layout(guides = "collect") & theme(legend.position = "right")

# Save
ggsave("combined_volcano_with_shared_legend.svg", combined_plot, width = 20, height = 6)

#~~~~ Combine boxplots into a single figure with unified legend and axis labels ~~~~#
# Define color mapping
color_mapping <- c("Black" = "grey3", "Yellow" = "yellow3", "Brown" = "brown4")

# Function to create a boxplot for a given dataset
create_boxplot <- function(expression_data, sample_sheet, gene_ids, gene_name_map, groups, title) {
  # Remove non-numeric columns
  expression_data <- expression_data[, sapply(expression_data, is.numeric)]
  
  # Ensure sample names match between expression data and sample sheet
  common_samples <- intersect(colnames(expression_data), sample_sheet$SampleID)
  expression_data <- expression_data[, common_samples]
  sample_sheet <- sample_sheet[sample_sheet$SampleID %in% common_samples, ]
  
  # Reorder sample sheet to match expression data
  sample_sheet <- sample_sheet[match(common_samples, sample_sheet$SampleID), ]
  
  # Log2 transform
  master_log2 <- log2(expression_data + 1)
  
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
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 16)) +
    labs(title = title, y = "log2(Expression + 1)", x = "Gene") +
    scale_fill_manual(values = color_mapping, drop = FALSE) +
    ylim(0,10)
  
  return(p)
}

# Load sample sheet
ss <- read.csv("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/02_reference_data/sample_sheet.csv")

# Load DGE data
yelvbla <- read.csv("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/06_dge_analysis_skin_colour/03_YelvBla/master.csv", row.names = 1)
brovbla <- read.csv("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/06_dge_analysis_skin_colour/01_BrovBla/master.csv", row.names = 1)
yelvbro <- read.csv("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/06_dge_analysis_skin_morph/02_YelvBroM/master.csv", row.names = 1)

# Define gene sets and mappings
genes_yelvbla <- c("rb_62964", "rb_79752", "rb_92274", "rb_78195", "rb_15163")
map_yelvbla <- c("rb_62964" = "TYR.I", "rb_79752" = "TYRP1", "rb_92274" = "MLANA", "rb_78195" = "PMEL", "rb_15163" = "DGAT2")

genes_brovbla <- c("rb_1531", "rb_78661", "rb_81552", "rb_15163")
map_brovbla <- c("rb_1531" = "TYR.II", "rb_78661" = "GPNMB.I", "rb_15163" = "DGAT2")

genes_yelvbro <- c("rb_76127", "rb_5385", "rb_78661", "rb_79752", "rb_92274", "rb_78195")
map_yelvbro <- c("rb_76127" = "PAH", "rb_5385" = "TXN", "rb_78661" = "GPNMB.I", "rb_79752" = "TYRP1", "rb_92274" = "MLANA", "rb_78195" = "PMEL")

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
  theme(legend.position = "right")


# Extract the legend
legend <- get_legend(dummy_plot)

# Combine plots with shared legend
combined_boxplot <- (plot1 | plot2 | plot3) + plot_layout(guides = "collect") & theme(legend.position = "none")

# Add the legend to the combined plot
final_plot <- plot_grid(combined_boxplot, legend, ncol = 2, rel_widths = c(3, 0.5))

# Save the final plot
ggsave("combined_boxplot_with_single_legend.svg", final_plot, width = 20, height = 6)

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
box1 <- box1 + ylab("log2(Expression + 1)") + xlab("")
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
  plot_grid(right_column, combined_legend, ncol = 2, rel_widths = c(3, 0.4)),
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
ggsave("final_combined_figure_with_legends_labeled.svg", final_labeled, width = 16, height = 10)
