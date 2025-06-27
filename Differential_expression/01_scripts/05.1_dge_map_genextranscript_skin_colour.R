#DGE analyses script
# by NS
#intended for use with output from edgeR2.R script

# Load necessary libraries
library(dplyr)

####~housekeeping~~~~~~~~~~~~####

rm(list=ls()) #clear the environment
setwd("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/01_scripts") #set wd
count = 00

###~~output dir~~~~~~~~~~~~~~####
output = "../06_dge_analysis_skin_colour/"
dir.create(output) #make folder for output
setwd(output)

# Load the txname_geneid.csv file
txname_geneid <- read.csv("../02_reference_data/txname_geneid.csv")
txname_geneid <- txname_geneid[ , -1]

# Define the analyses
analyses <- c("BrovBla", "YelvBro", "YelvBla")

# Loop through each analysis folder
for (count in seq_along(analyses)) {
  analysis <- analyses[count]
  
  # Determine the directory name
  if (count < 10) {
    anal_dir <- paste0("0", count, "_", analysis)
  } else {
    anal_dir <- paste0(count, "_", analysis)
  }
  
  # Change to the analysis directory
  setwd(anal_dir)
  
  # Define the gene list files
  gene_list_files <- c(
    "gene_universe.txt",
    "genes_sig.txt",
    "genes_non_sig.txt",
    "genes_sig_up.txt",
    "genes_sig_down.txt",
    "genes_sig_fdr.txt",
    "genes_fdr_up.txt",
    "genes_fdr_down.txt"
  )
  
  # Loop through each gene list file
  for (gene_list_file in gene_list_files) {
    # Check if the file exists
    if (file.exists(gene_list_file)) {
      # Load the gene list file
      gene_list <- readLines(gene_list_file)
      
      # Map transcript names to gene IDs
      gene_ids <- txname_geneid %>%
        filter(Transcript_ID %in% gene_list) %>%
        pull(GENEID)
      
      # Save the mapped gene IDs to a new file
      mapped_file <- paste0("mapped_", gene_list_file)
      writeLines(gene_ids, mapped_file)
    }
  }
  
  # Change back to the parent directory
  setwd("..")
}

