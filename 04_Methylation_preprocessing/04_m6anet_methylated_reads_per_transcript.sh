#!/bin/bash

# Define directory and samples
BASE_DIR="Base_dir"
SAMPLES=()

# Loop through to get m6a counts
for ana in "${SAMPLES[@]}"; do
  RESULT_DIR="${BASE_DIR}/${ana}/${ana}_m6anet_results"
  INPUT_FILE="${RESULT_DIR}/data.indiv_proba.csv"
  OUTPUT_FILE="${RESULT_DIR}/${ana}_m6acounts.csv"

  if [[ -f "$INPUT_FILE" ]]; then
    awk -F, '{if ($4 >= 0.033379376) { sites_per_transcript_read[$3"^"$1]++ } } END { for (i in sites_per_transcript_read ) { print i  }  }' "$INPUT_FILE" \
    | awk -F^ '{ count_per_transcript[$2]++ } END { for (i in count_per_transcript) { print i, count_per_transcript[i] }}' \
    > "$OUTPUT_FILE"
    echo "Processed: $ana"
  else
    echo "Missing file: $INPUT_FILE"
  fi
done
