#!/bin/bash

# Change to the directory containing the transcript name files
cd /Code || exit

# Define input and output files
transcript_names_files=(
  unique_transcripts_split_part_1.txt unique_transcripts_split_part_2.txt unique_transcripts_split_part_3.txt
  unique_transcripts_split_part_4.txt unique_transcripts_split_part_5.txt unique_transcripts_split_part_6.txt
  unique_transcripts_split_part_7.txt unique_transcripts_split_part_8.txt unique_transcripts_split_part_9.txt
  unique_transcripts_split_part_10.txt unique_transcripts_split_part_11.txt unique_transcripts_split_part_12.txt unique_transcripts_split_part_13.txt
)

transcriptome_file="/DirectRNA_Colour_bernardezi/NS9_transcriptome/rnabloom_transcriptome_kmer19.fasta"

# Process each transcript name file
for ana in "${transcript_names_files[@]}"; do
  output_file="/DirectRNA_Colour_bernardezi/${ana}.fa"
  > "$output_file"  # Clear or create the output file

  while IFS= read -r transcript_name; do
    grep -n "^>${transcript_name}[[:space:]]" "$transcriptome_file" | cut -d: -f1 | while read -r line_number; do
      awk -v start="$line_number" 'NR>=start {if (/^>/ && NR!=start) exit; print}' "$transcriptome_file" >> "$output_file"
    done
  done < "$ana"

  echo "Sequences have been extracted to $output_file"
done
