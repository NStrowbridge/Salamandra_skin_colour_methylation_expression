import csv

# Define input and output files
transcript_names_file = "/DirectRNA_Colour_bernardezi/Differential_expression/06_dge_analysis_skin_colour/03_YelvBla/genes_sig_fdr.txt"
preferred_name_file = "/DirectRNA_Colour_bernardezi/Differential_expression/06_dge_analysis_skin_colour/03_YelvBla/preferredname_genes_sig_fdr.csv"
output_file = "/DirectRNA_Colour_bernardezi/Differential_expression/06_dge_analysis_skin_colour/YelvBla_missing_genes_sigFDR.fa"

def compare_transcripts(transcript_file, csv_file, output_file):
    # Read the transcript names from the text file
    with open(transcript_file, 'r') as file:
        transcripts = file.read().splitlines()

    # Read the Transcript_IDs from the CSV file
    transcript_ids = set()
    with open(csv_file, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            transcript_ids.add(row['Transcript_ID'])

    # Find transcripts not in the CSV file
    missing_transcripts = [transcript for transcript in transcripts if transcript not in transcript_ids]

    # Write the missing transcripts to the output file
    with open(output_file, 'w') as file:
        for transcript in missing_transcripts:
            file.write(transcript + '\n')

# Example usage
compare_transcripts(transcript_names_file, preferred_name_file, output_file)
