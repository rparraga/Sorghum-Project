#!/bin/bash

# Define your reference files
METADATA="/software/projects/pawsey1157/${USER}/setonix/GitHub/Sorghum-Project/search/sorghum_runs_projects_merged_filtered.csv"  # The table with PRJNA numbers
LOG_DIR="/scratch/pawsey1157/modanilevicz/sorghum/nxf_work/data_stages/hisat2_align/"         # Directory containing your *.log files

echo "sample_id,overall_align_%,project" > "/software/projects/pawsey1157/${USER}/setonix/GitHub/Sorghum-Project/search/hisat2_alignment_summary.csv"

# Loop through every HISAT2 log file
for logfile in ${LOG_DIR}/*_align.log; do
    
    # 1. Extract the Sample ID from the filename (e.g., SRR1976113)
    sample_id=$(basename "$logfile" | cut -d'_' -f1)

    # 2. Extract the overall alignment rate using awk
    # Looks for the line containing "overall alignment rate" and takes the first word
    align_rate=$(awk '/overall alignment rate/ {sub(/%/, ""); print $1}' "$logfile")
    
    # 3. Match the Sample ID to the BioProject in the metadata table
    # We use awk to find the row where column 1 matches sample_id, then print column 2
    project=$(awk -F',' -v id="$sample_id" '$1 == id {print $2}' "$METADATA")

    # If the project wasn't found in the metadata, label it 'Unknown'
    if [ -z "$project" ]; then project="Unknown"; fi

    # 4. Output the CSV row
    echo "${sample_id},${align_rate},${project}" >> "/software/projects/pawsey1157/${USER}/setonix/GitHub/Sorghum-Project/search/hisat2_alignment_summary.csv"

done
