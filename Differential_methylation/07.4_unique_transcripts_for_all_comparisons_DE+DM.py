import pandas as pd

# Function to read 'Transcript_ID' and 'Preferred_name' columns and drop duplicates
def read_transcripts(file_path):
    return pd.read_csv(file_path, usecols=['Transcript_ID', 'Preferred_name']).drop_duplicates()

# Read and process each file
transcripts1 = read_transcripts('/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_methylation/06_dm_analysis_skin_transcript_level/01_BrovBla/preferredname_genes_sig_fdr.csv')
transcripts2 = read_transcripts('/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_methylation/06_dm_analysis_skin_transcript_level/03_YelvBla/preferredname_genes_sig_fdr.csv')
transcripts3 = read_transcripts('/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_methylation/06_dm_analysis_skin_transcript_level_morph/02_YelvBroM/preferredname_genes_sig_fdr.csv')
transcripts4 = read_transcripts('/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/06_dge_analysis_skin_colour/01_BrovBla/preferredname_genes_sig_fdr.csv')
transcripts5 = read_transcripts('/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/06_dge_analysis_skin_colour/03_YelvBla/preferredname_genes_sig_fdr.csv')
transcripts6 = read_transcripts('/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/06_dge_analysis_skin_morph/02_YelvBroM/preferredname_genes_sig_fdr.csv')

# Combine all and drop duplicates across the full dataset
combined_transcripts = pd.concat([transcripts1, transcripts2, transcripts3, transcripts4, transcripts5, transcripts6]).drop_duplicates()

# Save to file
combined_transcripts.to_csv('unique_transcripts_with_preferred_names_allcomparisons.txt', index=False)

print("Unique Transcript_ID and Preferred_name pairs saved to 'unique_transcripts_with_preferred_names_allcomparisons.txt'")
