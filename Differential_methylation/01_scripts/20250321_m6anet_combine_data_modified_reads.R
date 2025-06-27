#### m6Anet combine data and adjust to get modified reads ####
# Import dataframes, name them by their site, colour, and number
# Load necessary library
library(readr)
#Striped lateral
# Import the NS3_lat_CSV file
Striped_lat_m6a_NS3 <- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS3_lat/NS3_lat_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Striped_lat_m6a_NS3)
# Import the NS8_lat_CSV file
Striped_lat_m6a_NS8 <- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS8_lat/NS8_lat_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Striped_lat_m6a_NS8)
# Import the ELT_14082_lat_CSV file
Striped_lat_m6a_14082 <- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/ELT_14082_lat/ELT_14082_lat_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Striped_lat_m6a_14082)
# Import the ELT_14084_lat_CSV file
Striped_lat_m6a_14084 <- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/ELT_14084_lat/ELT_14084_lat_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Striped_lat_m6a_14084)

#Striped dorsal
# Import the NS3_dor_CSV file
Striped_dor_m6a_NS3 <- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS3_dor/NS3_dor_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Striped_dor_m6a_NS3)
# Import the NS8_dor_CSV file
Striped_dor_m6a_NS8 <- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS8_dor/NS8_dor_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Striped_dor_m6a_NS8)
# Import the ELT_14082_dor_CSV file
Striped_dor_m6a_14082 <- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/ELT_14082_dor/ELT_14082_dor_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Striped_dor_m6a_14082)
# Import the ELT_14084_dor_CSV file
Striped_dor_m6a_14084 <- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/ELT_14084_dor/ELT_14084_dor_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Striped_dor_m6a_14084)

#Yellow lateral
# Import the NS7_lat_CSV file
Yellow_lat_m6a_NS7 <- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS7_lat/NS7_lat_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Yellow_lat_m6a_NS7)
# Import the NS10_lat_CSV file
Yellow_lat_m6a_NS10 <- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS10_lat/NS10_lat_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Yellow_lat_m6a_NS10)
# Import the ELT_14079_lat_CSV file
Yellow_lat_m6a_14079 <- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/ELT_14079_lat/ELT_14079_lat_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Yellow_lat_m6a_14079)
# Import the ELT_14081_lat_CSV file
Yellow_lat_m6a_14081<- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/ELT_14081_lat/ELT_14081_lat_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Yellow_lat_m6a_14081)

#Yellow dorsal
# Import the NS7_dor_CSV file
Yellow_dor_m6a_NS7 <- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS7_dor/NS7_dor_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Yellow_dor_m6a_NS7)
# Import the NS10_dor_CSV file
Yellow_dor_m6a_NS10 <- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS10_dor/NS10_dor_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Yellow_dor_m6a_NS10)
# Import the ELT_14079_dor_CSV file
Yellow_dor_m6a_14079 <- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/ELT_14079_dor/ELT_14079_dor_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Yellow_dor_m6a_14079)
# Import the ELT_14081_dor_CSV file
Yellow_dor_m6a_14081 <- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/ELT_14081_dor/ELT_14081_dor_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Yellow_dor_m6a_14081)

#Brown lateral
# Import the NS2_lat_CSV file
Brown_lat_m6a_NS2 <- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS2_lat/NS2_lat_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Brown_lat_m6a_NS2)
# Import the NS5_lat_CSV file
Brown_lat_m6a_NS5 <- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS5_lat/NS5_lat_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Brown_lat_m6a_NS5)
# Import the NS11_lat_CSV file
Brown_lat_m6a_NS11 <- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS11_lat/NS11_lat_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Brown_lat_m6a_NS11)
# Import the NS12_lat_CSV file
Brown_lat_m6a_NS12 <- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS12_lat/NS12_lat_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Brown_lat_m6a_NS12)

#Brown dorsal
# Import the NS2_dor_CSV file
Brown_dor_m6a_NS2 <- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS2_dor/NS2_dor_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Brown_dor_m6a_NS2)
# Import the NS5_dor_CSV file
Brown_dor_m6a_NS5 <- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS5_dor/NS5_dor_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Brown_dor_m6a_NS5)
# Import the NS11_dor_CSV file
Brown_dor_m6a_NS11 <- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS11_dor/NS11_dor_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Brown_dor_m6a_NS11)
# Import the NS12_dor_CSV file
Brown_dor_m6a_NS12 <- read_csv("/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data/NS12_dor/NS12_dor_m6anet_results/data.site_proba.csv")
# View the first few rows of the dataframe to confirm import
head(Brown_dor_m6a_NS12)

####Combine dataframes into one ####
library(ggplot2)
# List of dataframes
file_paths <- list(
  Striped_lat_m6a_NS3 = Striped_lat_m6a_NS3,
  Striped_dor_m6a_NS3 = Striped_dor_m6a_NS3,
  Striped_lat_m6a_NS8 = Striped_lat_m6a_NS8,
  Striped_dor_m6a_NS8 = Striped_dor_m6a_NS8,
  Striped_lat_m6a_14082 = Striped_lat_m6a_14082,
  Striped_dor_m6a_14082 = Striped_dor_m6a_14082,
  Striped_lat_m6a_14084 = Striped_lat_m6a_14084,
  Striped_dor_m6a_14084 = Striped_dor_m6a_14084,
  Yellow_lat_m6a_NS7 = Yellow_lat_m6a_NS7, 
  Yellow_dor_m6a_NS7 = Yellow_dor_m6a_NS7,
  Yellow_lat_m6a_NS10 = Yellow_lat_m6a_NS10,
  Yellow_dor_m6a_NS10 = Yellow_dor_m6a_NS10,
  Yellow_lat_m6a_14079 = Yellow_lat_m6a_14079,
  Yellow_dor_m6a_14079 = Yellow_dor_m6a_14079,
  Yellow_lat_m6a_14081 = Yellow_lat_m6a_14081,
  Yellow_dor_m6a_14081 = Yellow_dor_m6a_14081,
  Brown_lat_m6a_NS2 = Brown_lat_m6a_NS2,
  Brown_dor_m6a_NS2 = Brown_dor_m6a_NS2, 
  Brown_lat_m6a_NS5 = Brown_lat_m6a_NS5,
  Brown_dor_m6a_NS5 = Brown_dor_m6a_NS5,
  Brown_lat_m6a_NS11 = Brown_lat_m6a_NS11,
  Brown_dor_m6a_NS11 = Brown_dor_m6a_NS11,
  Brown_lat_m6a_NS12 = Brown_lat_m6a_NS12,
  Brown_dor_m6a_NS12 = Brown_dor_m6a_NS12
)

# Function to label data
label_data <- function(data, sample_name) {
  data$sample <- sample_name
  return(data)
}

# Label and combine all samples into one dataframe
combined_data <- bind_rows(lapply(names(file_paths), function(name) label_data(file_paths[[name]], name)))

####Merge datasets only keeping common sites across all samples 
# Split combined data back into individual dataframes
library(dplyr)

# Split the combined data by sample
dfs <- split(combined_data, combined_data$sample)

library(dplyr)

rename_columns <- function(df, sample_name) {
  df <- as_tibble(df)  # Convert to tibble
  df <- df %>% dplyr::select(transcript_id, transcript_position, kmer, probability_modified, mod_ratio, n_reads)
  df <- df %>% mutate(
    mod_ratio = ifelse(probability_modified < 0.9, 0, mod_ratio),
    n_reads = ifelse(probability_modified < 0.9, 0, n_reads)
  )
  df <- df %>% rename_with(~paste(sample_name, ., sep = "_"), c(probability_modified, mod_ratio, n_reads))
  return(df)
}

# Rename columns for each dataframe
dfs <- lapply(names(dfs), function(name) rename_columns(dfs[[name]], name))

# Combine the dataframes back into one
merged_df_full <- Reduce(function(x, y) full_join(x, y, by = c("transcript_id", "transcript_position", "kmer")), dfs)

# Replace all NA values with zeros in the merged dataframe
merged_df_full[is.na(merged_df_full)] <- 0

# View the first few rows of the updated dataframe to confirm
head(merged_df_full)

# Create new column for each sample named "sample_name_n_mod_reads" by multiplying n_reads * mod_ratio
create_n_mod_reads_column <- function(df, sample_name) {
  df[[paste(sample_name, "n_mod_reads", sep = "_")]] <- df[[paste(sample_name, "n_reads", sep = "_")]] * df[[paste(sample_name, "mod_ratio", sep = "_")]]
  return(df)
}

# Apply the function to the merged dataframe
for (sample_name in names(file_paths)) {
  merged_df_full <- create_n_mod_reads_column(merged_df_full, sample_name)
}

# View the first few rows of the updated dataframe to confirm
head(merged_df_full)

# Save the merged dataframe to a CSV file
write.csv(merged_df_full, "merged_df_full.csv", row.names = FALSE)

# Create a new dataframe with only n_mod_reads columns
n_mod_reads_columns <- grep("_n_mod_reads$", colnames(merged_df_full), value = TRUE)
n_mod_reads_df <- merged_df_full[, c("transcript_id", "transcript_position", "kmer", n_mod_reads_columns)]

# Convert tibble to data frame
n_mod_reads_df <- as.data.frame(n_mod_reads_df)

# Sum n_mod_reads for each sample at each transcript_id
sum_n_mod_reads <- function(df, sample_name) {
  df %>%
    group_by(transcript_id) %>%
    summarise(!!paste(sample_name, "total_n_mod_reads", sep = "_") := sum(!!sym(paste(sample_name, "n_mod_reads", sep = "_")), na.rm = TRUE))
}

# Apply the function to each sample
summed_dfs <- lapply(names(file_paths), function(name) sum_n_mod_reads(n_mod_reads_df, name))

# Combine the summed dataframes back into one
merged_summed_df <- Reduce(function(x, y) full_join(x, y, by = "transcript_id"), summed_dfs)

# View the first few rows of the updated dataframe to confirm
head(merged_summed_df)

# Convert tibble to data frame
filtered_merged_summed_df <- as.data.frame(merged_summed_df)

# Create row names by combining transcript_id
filtered_merged_summed_df$row_name <- paste(filtered_merged_summed_df$transcript_id)

# Set row names
rownames(filtered_merged_summed_df) <- filtered_merged_summed_df$row_name

filtered_merged_summed_df <- filtered_merged_summed_df[, -c(1, ncol(filtered_merged_summed_df))]  # Remove the original columns and the row_name column

# View the first few rows of the updated dataframe to confirm
head(filtered_merged_summed_df)

# Save the n_mod_reads dataframe to a CSV file
write.csv(filtered_merged_summed_df, "n_mod_reads_summedxtranscript_df.csv", row.names = TRUE)

# Make row names be a combination of transcript_id, transcript_position, and kmer
n_mod_reads_df$row_name <- paste(n_mod_reads_df$transcript_id, n_mod_reads_df$transcript_position, n_mod_reads_df$kmer, sep = "_")
rownames(n_mod_reads_df) <- n_mod_reads_df$row_name
n_mod_reads_df <- n_mod_reads_df[, -c(1:3, ncol(n_mod_reads_df))]  # Remove the original columns and the row_name column

# View the first few rows of the updated dataframe to confirm
head(n_mod_reads_df)

# Save the n_mod_reads dataframe to a CSV file
write.csv(n_mod_reads_df, "n_mod_reads_df.csv", row.names = TRUE)

# Filter rows that contain at least >= 5 in at least 60% of samples (14 samples)
filtered_n_mod_reads_df <- n_mod_reads_df[rowSums(n_mod_reads_df >= 5) >= 13, ]

# View the first few rows of the updated dataframe to confirm
head(filtered_n_mod_reads_df)

# Save the filtered n_mod_reads dataframe to a CSV file
write.csv(filtered_n_mod_reads_df, "filtered_n_mod_reads_df.csv", row.names = TRUE)



