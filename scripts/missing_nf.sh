#!/bin/bash

# Usage: ./check_samples_report.sh <results_dir> <sample_list.csv>
RESULTS_DIR=$1
CSV_FILE=$2
REPORT_FILE="${PWD}/missing_samples_report.txt"

# Check for required arguments
if [[ -z "$RESULTS_DIR" || -z "$CSV_FILE" ]]; then
    echo "Usage: $0 <results_dir> <sample_list.csv>"
    exit 1
fi

# Initialize/Clear the report file
echo "NEXTFLOW AUDIT - $(date)" > "$REPORT_FILE"
echo "------------------------------------------------" >> "$REPORT_FILE"
printf "%-20s | %-40s\n" "SAMPLE_ID" "EXPECTED_LOCATION" >> "$REPORT_FILE"
echo "------------------------------------------------" >> "$REPORT_FILE"

# Use mapfile for safer array creation from CSV
mapfile -t SAMPLES < <(tail -n +2 "$CSV_FILE" | cut -d',' -f1 | tr -d '\r')
MISSING_COUNT=0
TOTAL_COUNT=0

# Loop through every directory inside the RESULTS_DIR
for FOLDER_PATH in "$RESULTS_DIR"/*/; do
    # Remove trailing slash for cleaner reporting
    FOLDER_NAME=$(basename "$FOLDER_PATH")

    for SAMPLE in "${SAMPLES[@]}"; do
        # Use a glob check that works
        # This checks if any file matching the pattern exists and has size > 0
        MATCHES=($FOLDER_PATH/${SAMPLE}*)
        
        if [[ -e "${MATCHES[0]}" && -s "${MATCHES[0]}" ]]; then
            # Success, file exists and is not empty
            : 
        else
            ((MISSING_COUNT++))
            echo " - Missing or Empty: ${FOLDER_NAME}/${SAMPLE}" >> "$REPORT_FILE"
        fi
        ((TOTAL_COUNT++))
    done
done

# Append summary footer to report
echo "------------------------------------------------" >> "$REPORT_FILE"
echo "Summary: $MISSING_COUNT out of $TOTAL_COUNT checks failed (missing/empty)." >> "$REPORT_FILE"

echo "Done! Report generated: $REPORT_FILE"
