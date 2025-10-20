#!/bin/bash

# Define directory and samples
BASE_DIR="/Volumes/LaCie/Nic_DirectRNA_raw_data/m6anet_data"
SAMPLES=("NS2_dor" "NS2_lat" "NS3_dor" "NS3_lat" "NS5_dor" "NS5_lat" "NS7_dor" "NS7_lat" "NS8_dor" "NS8_lat" "NS10_dor" "NS10_lat" "NS11_dor" "NS11_lat" "NS12_dor" "NS12_lat" "ELT_14079_dor" "ELT_14079_lat" "ELT_14081_dor" "ELT_14081_lat" "ELT_14082_dor" "ELT_14082_lat" "ELT_14084_dor" "ELT_14084_lat")

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
