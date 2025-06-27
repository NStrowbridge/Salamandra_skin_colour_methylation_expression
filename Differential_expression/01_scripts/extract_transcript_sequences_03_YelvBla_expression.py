def extract_transcripts(transcript_names_file, transcriptome_file, output_file):
    # Read the transcript names into a set
    with open(transcript_names_file, 'r') as f:
        transcript_names = set(line.strip() for line in f)

    # Open the transcriptome file and the output file
    with open(transcriptome_file, 'r') as transcriptome, open(output_file, 'w') as output:
        write_sequence = False
        for line in transcriptome:
            if line.startswith('>'):
                # Check if the header contains any of the transcript names
                transcript_id = line.split()[0][1:]  # Extract the transcript ID from the header
                if transcript_id in transcript_names:
                    write_sequence = True
                    output.write(line)  # Write the header to the output file
                else:
                    write_sequence = False
            elif write_sequence:
                output.write(line)  # Write the sequence to the output file

# Define the file paths
transcript_names_file="/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/06_dge_analysis_skin_colour/03_YelvBla/genes_sig_fdr.txt"
transcriptome_file="/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/NS9_transcriptome/rnabloom_transcriptome_kmer19.fasta"
output_file="/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/06_dge_analysis_skin_colour/YelvBla_extracted_sequences_sigFDR.fa"

# Call the function with the defined file paths
extract_transcripts(transcript_names_file, transcriptome_file, output_file)
