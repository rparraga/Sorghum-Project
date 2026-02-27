#!/bin/bash --login

#SBATCH --account=pawsey1157
#SBATCH --nodes=1               # each subtask uses 1 node
#SBATCH --ntasks=1              # 1 subtask per file in the array-subtask
#SBATCH --cpus-per-task=8      # cpus per subtask
#SBATCH --mem=10GB              # Needed memory per array-subtask 
#SBATCH --time=24:00:00         # time per subtask
#SBATCH --partition=copy

set -euo pipefail

# --- Configuration & Argument Checks ---
module load rclone/1.68.1

# 1. Check for Target Directory
if [ -z "$1" ]; then
    echo "Usage: $0 <source_directory> <rclone_remote_name> <rclone_bucket_path>"
    echo "Example: $0 /home/user/data my_s3_remote my-backups/data_zips"
    exit 1
fi

TARGET_DIR="$1"
RCLONE_REMOTE="$2"
RCLONE_BUCKET="$3"

echo "1/2 Zipping folder: ${TARGET_DIR} to ${TARGET_DIR}.zip..."
zip -r -q "${TARGET_DIR}.zip" "${TARGET_DIR}"

# Backup the compressed file using rclone 
rclone copy "${TARGET_DIR}.zip" "${RCLONE_REMOTE}:${RCLONE_BUCKET}" --progress --stats 1s

echo "--> Backup Success: '$ZIP_FILE_PATH' successfully uploaded."
            
