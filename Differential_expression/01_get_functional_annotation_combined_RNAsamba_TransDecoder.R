#Get functional gene annotation using eggnogg protein mappings 
# overlapping annotation created with Transdecoder then blastn and RNAsamba and blastp (eggnogmapper: diamond)
# R script
# by Nic

####~load libraries~~~~~~~~~~####
library(readxl)
library(dplyr)

####~housekeeping~~~~~~~~~~~~####
rm(list=ls()) #clear the environment
setwd("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/02_reference_data") #set wd to reference data folder

####~get eggnogg annotation~~####
#transdecoder
eggnog_transdecoder = read_excel("eggnogmapper_diamond_output_transdecoder_kmer19_transcriptID.emapper.annotations.xlsx", skip = 4)
eggnog_transdecoder = eggnog_transdecoder[-c(19103,19104,19105),] #delete empty rows and header row
colnames(eggnog_transdecoder)[1] = "query"
eggnog_transdecoder <- eggnog_transdecoder %>%
  group_by(Transcript_ID) %>%
  slice_max(score, n = 1) %>%
  ungroup()

#rnasamba
eggnog_rnasamba <- read.delim("eggnogmapper_diamond_output_rnasamba_kmer19.emapper.annotations.tsv", skip = 4)
eggnog_rnasamba = eggnog_rnasamba[-c(18948,18949,18950),] #delete empty rows and header row
colnames(eggnog_rnasamba)[1] = "Transcript_ID" 

# Function to overlap two dataframes by Transcript_ID and check overlaps in other columns
overlap_dataframes <- function(df1, df2) {
  # Merge dataframes on Transcript_ID column with all rows from both dataframes
  merged_df <- merge(df1, df2, by = "Transcript_ID", suffixes = c("_transdecoder", "_rnasamba"), all = TRUE)
  
  # Check overlaps in other columns
  overlap_columns <- intersect(names(df1), names(df2))
  overlap_columns <- setdiff(overlap_columns, "Transcript_ID")
  
  match_counts <- data.frame(Column = character(), TRUE_Count = integer(), FALSE_Count = integer(), NA_Count = integer(), stringsAsFactors = FALSE)
  for (col in overlap_columns) {
    match_col <- paste0(col, "_match")
    merged_df[[match_col]] <- merged_df[[paste0(col, "_transdecoder")]] == merged_df[[paste0(col, "_rnasamba")]]
    counts <- table(factor(merged_df[[match_col]], levels = c(TRUE, FALSE, NA)))
    match_counts <- rbind(match_counts, data.frame(Column = match_col, TRUE_Count = counts["TRUE"], FALSE_Count = counts["FALSE"], NA_Count = sum(is.na(merged_df[[match_col]]))))
  }
  
  # Reorder columns in merged dataframe
  reordered_columns <- c("Transcript_ID")
  for (col in overlap_columns) {
    reordered_columns <- c(reordered_columns, paste0(col, "_transdecoder"), paste0(col, "_rnasamba"))
  }
  reordered_columns <- c(reordered_columns, grep("_match$", names(merged_df), value = TRUE))
  merged_df <- merged_df[, reordered_columns]
  
  list(merged_df = merged_df, match_counts = match_counts)
}

# Example usage with the same dataframe for demonstration
result <- overlap_dataframes(eggnog_transdecoder, eggnog_rnasamba)

# Now making dataframe with final annotated transcripts 

final_annotated_transcripts <- function(merged_df) {
  # Define the columns to retain
  columns_to_retain <- c("Transcript_ID", "seed_ortholog", "Preferred_name", "GOs", "Description", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "PFAMs",  "score_transdecoder", "score_rnasamba")
  
  # Create a new dataframe to store the final list of annotated transcripts
  final_df <- merged_df %>%
    mutate(
      seed_ortholog = ifelse(
        !is.na(seed_ortholog_transdecoder) & !is.na(seed_ortholog_rnasamba),
        ifelse(score_transdecoder >= score_rnasamba,
               seed_ortholog_transdecoder, seed_ortholog_rnasamba),
        ifelse(!is.na(seed_ortholog_transdecoder),
               seed_ortholog_transdecoder, seed_ortholog_rnasamba)
      ),
      
      Preferred_name = ifelse(
        !is.na(Preferred_name_transdecoder) & !is.na(Preferred_name_rnasamba),
        ifelse(score_transdecoder >= score_rnasamba,
               Preferred_name_transdecoder, Preferred_name_rnasamba),
        ifelse(!is.na(Preferred_name_transdecoder),
               Preferred_name_transdecoder, Preferred_name_rnasamba)
      ),
      
      GOs = ifelse(
        !is.na(GOs_transdecoder) & !is.na(GOs_rnasamba),
        ifelse(score_transdecoder >= score_rnasamba,
               GOs_transdecoder, GOs_rnasamba),
        ifelse(!is.na(GOs_transdecoder),
               GOs_transdecoder, GOs_rnasamba)
      ),
      
      Description = ifelse(
        !is.na(Description_transdecoder) & !is.na(Description_rnasamba),
        ifelse(score_transdecoder >= score_rnasamba,
               Description_transdecoder, Description_rnasamba),
        ifelse(!is.na(Description_transdecoder),
               Description_transdecoder, Description_rnasamba)
      ),
      
      KEGG_ko = ifelse(
        !is.na(KEGG_ko_transdecoder) & !is.na(KEGG_ko_rnasamba),
        ifelse(score_transdecoder >= score_rnasamba,
               KEGG_ko_transdecoder, KEGG_ko_rnasamba),
        ifelse(!is.na(KEGG_ko_transdecoder),
               KEGG_ko_transdecoder, KEGG_ko_rnasamba)
      ),
      
      KEGG_Pathway = ifelse(
        !is.na(KEGG_Pathway_transdecoder) & !is.na(KEGG_Pathway_rnasamba),
        ifelse(score_transdecoder >= score_rnasamba,
               KEGG_Pathway_transdecoder, KEGG_Pathway_rnasamba),
        ifelse(!is.na(KEGG_Pathway_transdecoder),
               KEGG_Pathway_transdecoder, KEGG_Pathway_rnasamba)
      ),
      
      KEGG_Module = ifelse(
        !is.na(KEGG_Module_transdecoder) & !is.na(KEGG_Module_rnasamba),
        ifelse(score_transdecoder >= score_rnasamba,
               KEGG_Module_transdecoder, KEGG_Module_rnasamba),
        ifelse(!is.na(KEGG_Module_transdecoder),
               KEGG_Module_transdecoder, KEGG_Module_rnasamba)
      ),
      
      KEGG_Reaction = ifelse(
        !is.na(KEGG_Reaction_transdecoder) & !is.na(KEGG_Reaction_rnasamba),
        ifelse(score_transdecoder >= score_rnasamba,
               KEGG_Reaction_transdecoder, KEGG_Reaction_rnasamba),
        ifelse(!is.na(KEGG_Reaction_transdecoder),
               KEGG_Reaction_transdecoder, KEGG_Reaction_rnasamba)
      ),
      
      KEGG_rclass = ifelse(
        !is.na(KEGG_rclass_transdecoder) & !is.na(KEGG_rclass_rnasamba),
        ifelse(score_transdecoder >= score_rnasamba,
               KEGG_rclass_transdecoder, KEGG_rclass_rnasamba),
        ifelse(!is.na(KEGG_rclass_transdecoder),
               KEGG_rclass_transdecoder, KEGG_rclass_rnasamba)
      ),
      
      PFAMs = ifelse(
        !is.na(PFAMs_transdecoder) & !is.na(PFAMs_rnasamba),
        ifelse(score_transdecoder >= score_rnasamba,
               PFAMs_transdecoder, PFAMs_rnasamba),
        ifelse(!is.na(PFAMs_transdecoder),
               PFAMs_transdecoder, PFAMs_rnasamba)
      )
    ) %>%
    dplyr::select(
      Transcript_ID, seed_ortholog, Preferred_name, GOs,
      Description, KEGG_ko, KEGG_Pathway, KEGG_Module,
      KEGG_Reaction, KEGG_rclass, PFAMs
    ) %>%
    filter(
      !is.na(seed_ortholog),
      !is.na(Preferred_name),
      !is.na(GOs)
    )
  
  return(final_df)
}


# Example usage with the merged dataframe from the previous function
result <- overlap_dataframes(eggnog_transdecoder, eggnog_rnasamba)
final_result <- final_annotated_transcripts(result$merged_df)

duplicates <- final_result %>% filter(duplicated(Transcript_ID))
print("\nDuplicates in Transcript_ID column:")
print(duplicates)

transcriptome_annotaton_combined_RNAsamba_Transdecoder <- final_result %>% distinct(Transcript_ID, .keep_all = TRUE)


transcriptome_annotaton_combined_RNAsamba_Transdecoder_dups <- transcriptome_annotaton_combined_RNAsamba_Transdecoder %>% filter(duplicated(Transcript_ID))
print("\nDuplicates in Transcript_ID column:")
print(transcriptome_annotaton_combined_RNAsamba_Transdecoder)

###~~make txname*geneid~~~~~~~~~~####
# Rename seed_ortholog to gene name
transcriptome_annotaton_combined_RNAsamba_Transdecoder <- transcriptome_annotaton_combined_RNAsamba_Transdecoder %>% rename("GENEID" = "seed_ortholog")
print(head(transcriptome_annotaton_combined_RNAsamba_Transdecoder))
txname_geneid = transcriptome_annotaton_combined_RNAsamba_Transdecoder[,c("Transcript_ID","GENEID")]

###~~make annotation files~~~~~~~~~~####
txname_Preferredname_GOcat = transcriptome_annotaton_combined_RNAsamba_Transdecoder[,c("Transcript_ID","Preferred_name", "GOs")]
txname_Preferredname = transcriptome_annotaton_combined_RNAsamba_Transdecoder[,c("Transcript_ID","Preferred_name")]
txname_all_annotations = transcriptome_annotaton_combined_RNAsamba_Transdecoder[,c("Transcript_ID", "GENEID", "Preferred_name", "GOs", "Description", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "PFAMs")]

####~map file for topGO~~~~~~####
gene2GO = transcriptome_annotaton_combined_RNAsamba_Transdecoder[,c("GENEID","GOs")] #make gene2GO annotation list

####~write out annotation~~~~####
write.csv(transcriptome_annotaton_combined_RNAsamba_Transdecoder, file = "../02_reference_data/transcriptome_annotation_combined_RNAsamba_Transdecoder.csv", row.names = FALSE)
write.csv(txname_geneid, file = "../02_reference_data/txname_geneid.csv")
write.csv(txname_Preferredname_GOcat, file = "../02_reference_data/txname_Preferredname_GOcat.csv")
write.csv(txname_Preferredname, file = "../02_reference_data/txname_Preferredname.csv")
write.csv(txname_all_annotations, file = "../02_reference_data/txname_all_annotations.csv")
write_tsv(gene2GO, file = "../02_reference_data/gene2GO.map", col_names = FALSE)
write.csv(result$merged_df, file = "../02_reference_data/merged_df.csv")
####~fin~~~~~~~~~~~~~~~~~~~~~####
closeAllConnections()
