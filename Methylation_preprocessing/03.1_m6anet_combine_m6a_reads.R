#### m6Anet combine data and adjust to get modified reads ####

####~housekeeping~~~~~~~~~~~~####
rm(list=ls()) #clear the environment
setwd("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_methylation/01_scripts") #set wd to Scripts folder

###~~output directory~~~~~~~~####
output = "../04_n_methylated_reads" #specify where the output should go
dir.create(output) #create directory for output
setwd(output) #set the new output directory as the working directory

# Load necessary library
library(readr)
library(dplyr)
library(purrr)
library(tibble)


# Import dataframes with space as separator and rename columns
NS3_lat_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS3_lat/NS3_lat_m6anet_results/NS3_lat_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "NS3_lat"))
NS8_lat_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS8_lat/NS8_lat_m6anet_results/NS8_lat_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "NS8_lat"))
ELT_14082_lat_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/ELT_14082_lat/ELT_14082_lat_m6anet_results/ELT_14082_lat_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "ELT_14082_lat"))
ELT_14084_lat_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/ELT_14084_lat/ELT_14084_lat_m6anet_results/ELT_14084_lat_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "ELT_14084_lat"))

NS3_dor_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS3_dor/NS3_dor_m6anet_results/NS3_dor_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "NS3_dor"))
NS8_dor_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS8_dor/NS8_dor_m6anet_results/NS8_dor_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "NS8_dor"))
ELT_14082_dor_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/ELT_14082_dor/ELT_14082_dor_m6anet_results/ELT_14082_dor_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "ELT_14082_dor"))
ELT_14084_dor_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/ELT_14084_dor/ELT_14084_dor_m6anet_results/ELT_14084_dor_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "ELT_14084_dor"))

NS7_lat_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS7_lat/NS7_lat_m6anet_results/NS7_lat_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "NS7_lat"))
NS10_lat_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS10_lat/NS10_lat_m6anet_results/NS10_lat_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "NS10_lat"))
ELT_14079_lat_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/ELT_14079_lat/ELT_14079_lat_m6anet_results/ELT_14079_lat_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "ELT_14079_lat"))
ELT_14081_lat_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/ELT_14081_lat/ELT_14081_lat_m6anet_results/ELT_14081_lat_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "ELT_14081_lat"))

NS7_dor_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS7_dor/NS7_dor_m6anet_results/NS7_dor_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "NS7_dor"))
NS10_dor_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS10_dor/NS10_dor_m6anet_results/NS10_dor_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "NS10_dor"))
ELT_14079_dor_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/ELT_14079_dor/ELT_14079_dor_m6anet_results/ELT_14079_dor_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "ELT_14079_dor"))
ELT_14081_dor_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/ELT_14081_dor/ELT_14081_dor_m6anet_results/ELT_14081_dor_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "ELT_14081_dor"))

NS2_lat_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS2_lat/NS2_lat_m6anet_results/NS2_lat_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "NS2_lat"))
NS5_lat_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS5_lat/NS5_lat_m6anet_results/NS5_lat_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "NS5_lat"))
NS11_lat_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS11_lat/NS11_lat_m6anet_results/NS11_lat_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "NS11_lat"))
NS12_lat_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS12_lat/NS12_lat_m6anet_results/NS12_lat_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "NS12_lat"))

NS2_dor_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS2_dor/NS2_dor_m6anet_results/NS2_dor_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "NS2_dor"))
NS5_dor_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS5_dor/NS5_dor_m6anet_results/NS5_dor_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "NS5_dor"))
NS11_dor_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS11_dor/NS11_dor_m6anet_results/NS11_dor_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "NS11_dor"))
NS12_dor_m6a_reads <- read_delim("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS12_dor/NS12_dor_m6anet_results/NS12_dor_m6acounts.csv", delim = " ", col_names = c("Transcript_ID", "NS12_dor"))


# Combine dataframes by Transcript_ID
combined_df <- list(NS3_lat_m6a_reads, NS8_lat_m6a_reads, ELT_14082_lat_m6a_reads, ELT_14084_lat_m6a_reads,
                    NS3_dor_m6a_reads, NS8_dor_m6a_reads, ELT_14082_dor_m6a_reads, ELT_14084_dor_m6a_reads,
                    NS7_lat_m6a_reads, NS10_lat_m6a_reads, ELT_14079_lat_m6a_reads, ELT_14081_lat_m6a_reads,
                    NS7_dor_m6a_reads, NS10_dor_m6a_reads, ELT_14079_dor_m6a_reads, ELT_14081_dor_m6a_reads,
                    NS2_lat_m6a_reads, NS5_lat_m6a_reads, NS11_lat_m6a_reads, NS12_lat_m6a_reads,
                    NS2_dor_m6a_reads, NS5_dor_m6a_reads, NS11_dor_m6a_reads, NS12_dor_m6a_reads) %>%
  reduce(full_join, by = "Transcript_ID")

# Replace NAs with 0s
combined_df[is.na(combined_df)] <- 0

# Change column 1 to row names
combined_df <- column_to_rownames(combined_df, var = "Transcript_ID")

# Save combined_df as n_methylated_reads.csv
write.csv(combined_df, "n_methylated_reads.csv", row.names = TRUE)


