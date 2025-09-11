#### Script to analyze striped yellow vs striped black DE & DM
#### By Nic

####~load libraries~~~~~~~~~~####
library(ggplot2)
library(ggrepel)
library(dplyr)

####~housekeeping~~~~~~~~~~~~####
rm(list = ls()) # clear the environment
setwd("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/01_scripts") # set wd

####~~output dir~~~~~~~~~~~~~~####
output <- "../05.3_dge+dm_stridorvstrilat_pvalue_plot_annotation/"
dir.create(output, showWarnings = FALSE) # make folder for output
setwd(output)

####~import edgeR results~~~~~~~~####
dge_path <- "../04_edge_R_transcript_level_skin_stridorxstrilat/results_Stri_dorvStri_lat.csv"
dge_data <- read.csv(dge_path, header = TRUE, stringsAsFactors = FALSE)

DM_path <- "../../Differential_methylation/04_edgeR_methylation_transcript_level_stridorvstrilat/results_Stri_dorvStri_lat.csv"
DM_data <- read.csv(DM_path, header = TRUE, stringsAsFactors = FALSE)

####~rename transcript column if needed~~~~~~~~####
colnames(DM_data)[1] <- "Transcript_ID"
colnames(dge_data)[1] <- "Transcript_ID"

####~flag logFC > 6 or < -6~~~~~~~~####
dge_data$logFC_gt6_or_lt_neg6 <- abs(dge_data$logFC) > 6
DM_data$logFC_gt6_or_lt_neg6 <- abs(DM_data$logFC) > 6

####~plot PValue vs logFC (DE)~~~~~~~~####
pval_fc_plot_DE <- ggplot(dge_data, aes(x = logFC, y = PValue, color = logFC_gt6_or_lt_neg6)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "firebrick")) +
  theme_minimal(base_size = 14) +
  labs(
    x = "log2(Fold Change)",
    y = "Expression (PValue)",
    color = "logFC > 6 or < -6"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
  )

####~plot PValue vs logFC (DM)~~~~~~~~####
pval_fc_plot_DM <- ggplot(DM_data, aes(x = logFC, y = PValue, color = logFC_gt6_or_lt_neg6)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "firebrick")) +
  theme_minimal(base_size = 14) +
  labs(
    x = "log2(Fold Change)",
    y = "Methylation (PValue)",
    color = "logFC > 6 or < -6"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
  )

####~save plots~~~~~~~~~~~~~####
ggsave("PValue_vs_logFC_plot_DE_colored.png", pval_fc_plot_DE, width = 8, height = 6, dpi = 300)
ggsave("PValue_vs_logFC_plot_DM_colored.png", pval_fc_plot_DM, width = 8, height = 6, dpi = 300)

#### Script to extract and annotate transcripts with logFC > 6 or < -6
####~load transcript-to-gene mapping~~~~~~~~####
txname_preferredname <- read.csv("../02_reference_data/txname_preferredname.csv", stringsAsFactors = FALSE)
txname_preferredname <- txname_preferredname[, -1] # drop first column if it's an index

####~filter transcripts with logFC > 6 or < -6~~~~~~~~####
dge_filtered <- dge_data %>% filter(abs(logFC) > 6)
DM_filtered  <- DM_data %>% filter(abs(logFC) > 6)

####~retain only first column~~~~~~~~####
dge_filtered <- dge_filtered[, 1, drop = FALSE]
DM_filtered  <- DM_filtered[, 1, drop = FALSE]
                            
####~annotate filtered transcripts~~~~~~~~####
dge_annotated <- dge_filtered %>%
  left_join(txname_preferredname, by = "Transcript_ID")

DM_annotated <- DM_filtered %>%
  left_join(txname_preferredname, by = "Transcript_ID")

####~save annotated tables~~~~~~~~####
write.csv(dge_annotated, "DGE_transcripts_logFC_gt6_or_lt_neg6_annotated.csv", row.names = FALSE)
write.csv(DM_annotated, "DM_transcripts_logFC_gt6_or_lt_neg6_annotated.csv", row.names = FALSE)

####~extract rows with NA Preferred_name~~~~~~~~####
dge_to_annotate <- dge_annotated %>% filter(is.na(Preferred_name))
DM_to_annotate  <- DM_annotated %>% filter(is.na(Preferred_name))

####~save transcripts needing annotation~~~~~~~~####
write.csv(dge_to_annotate, "DGE_to_annotate.csv", row.names = FALSE)
write.csv(DM_to_annotate, "DM_to_annotate.csv", row.names = FALSE)

#Manually annotated the outputs without annotations using blast search, now reimport to join to others
####~DGE: re-import and join~~~~~~~~####
DGE_auto <- read.csv("DGE_transcripts_logFC_gt6_or_lt_neg6_annotated.csv", stringsAsFactors = FALSE)
DGE_manual <- read.csv("DGE_manually_annotate_stridorvstrilat.csv", stringsAsFactors = FALSE)

DGE_combined <- DGE_auto %>%
  left_join(DGE_manual %>% select(Transcript_ID, Preferred_name_manual = Preferred_name), by = "Transcript_ID") %>%
  mutate(Preferred_name = ifelse(!is.na(Preferred_name_manual), Preferred_name_manual, Preferred_name)) %>%
  select(-Preferred_name_manual)

write.csv(DGE_combined, "DGE_transcripts_logFC_gt6_or_lt_neg6_combined_annotation.csv", row.names = FALSE)

####~DM: re-import and join~~~~~~~~####
DM_auto <- read.csv("DM_transcripts_logFC_gt6_or_lt_neg6_annotated.csv", stringsAsFactors = FALSE)
DM_manual <- read.csv("DM_manually_annotate_stridorvstrilat.csv", stringsAsFactors = FALSE)

DM_combined <- DM_auto %>%
  left_join(DM_manual %>% select(Transcript_ID, Preferred_name_manual = Preferred_name), by = "Transcript_ID") %>%
  mutate(Preferred_name = ifelse(!is.na(Preferred_name_manual), Preferred_name_manual, Preferred_name)) %>%
  select(-Preferred_name_manual)

write.csv(DM_combined, "DM_transcripts_logFC_gt6_or_lt_neg6_combined_annotation.csv", row.names = FALSE)

####~rejoin to original data and label expression~~~~~~~~####
DM_data_annotated <- DM_data %>%
  inner_join(DM_combined, by = "Transcript_ID") %>%
  mutate(expression_label = ifelse(logFC < 0, "Hypo",
                                   ifelse(logFC > 0, "Hyper", NA)))

dge_data_annotated <- dge_data %>%
  inner_join(DGE_combined, by = "Transcript_ID") %>%
  mutate(expression_label = ifelse(logFC < 0, "Down",
                                   ifelse(logFC > 0, "Up", NA)))

####~extract summary columns~~~~~~~~####
DM_summary <- DM_data_annotated %>%
  select(Transcript_ID, Preferred_name, expression_label)

write.csv(DM_summary, "DM_transcripts_logFC_gt6_or_lt_neg6_summary.csv", row.names = FALSE)

DGE_summary <- dge_data_annotated %>%
  select(Transcript_ID, Preferred_name, expression_label)

write.csv(DGE_summary, "DGE_transcripts_logFC_gt6_or_lt_neg6_summary.csv", row.names = FALSE)



