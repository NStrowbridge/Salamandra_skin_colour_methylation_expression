import pandas as pd
import numpy as np

# Load the transcript IDs from the file
file_path = 'unique_transcripts_for_annotation_allcomparisons_expression+methylation.txt'
transcripts = pd.read_csv(file_path, header=None, names=['TranscriptID'])

# Determine how many rows each split should have
num_splits = 13
rows_per_split = int(np.ceil(len(transcripts) / num_splits))

# Split the DataFrame and save each part
for i in range(num_splits):
    split_df = transcripts.iloc[i * rows_per_split : (i + 1) * rows_per_split]
    split_df.to_csv(f'unique_transcripts_split_part_{i+1}.txt', index=False, header=False)

print("The transcript IDs have been split into 13 parts and saved.")
