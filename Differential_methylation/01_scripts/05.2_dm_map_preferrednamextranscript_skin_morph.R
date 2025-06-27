#Mapping preferred name to transcript name analyses script
# by NS

# Load necessary libraries
library(dplyr)

####~housekeeping~~~~~~~~~~~~####

rm(list=ls()) #clear the environment
setwd("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_methylation/01_scripts") #set wd
count = 00

###~~output dir~~~~~~~~~~~~~~####
output = "../06_dm_analysis_skin_transcript_level_morph/"
dir.create(output) #make folder for output
setwd(output)

# Load the txname_preferredname.csv file
txname_preferredname <- read.csv("../02_reference_data/txname_preferredname.csv")
txname_preferredname <- txname_preferredname[ , -1]

# Define the analyses
analyses <- c("BrovStri", "YelvBroM", "YelvStri")

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
  print(paste("Changing to directory:", anal_dir))
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
      print(paste("Processing file:", gene_list_file))
      # Load the gene list file
      gene_list <- readLines(gene_list_file)
      
      # Map transcript names to gene IDs
      mapped_genes <- txname_preferredname %>%
        filter(Transcript_ID %in% gene_list) %>%
        dplyr::select(Transcript_ID, Preferred_name)
      
      # Save the mapped gene IDs to a new file with .csv extension
      mapped_file <- paste0("preferredname_", sub(".txt", ".csv", gene_list_file))
      write.csv(mapped_genes, mapped_file, row.names = FALSE)
    } else {
      print(paste("File not found:", gene_list_file))
    }
  }
  
  # Change back to the parent directory
  setwd("..")
}
