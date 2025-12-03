#!/bin/bash --login

#SBATCH --account=pawsey1157
#SBATCH --nodes=1               # each subtask uses 1 node
#SBATCH --ntasks=1              # 1 subtask per file in the array-subtask
#SBATCH --cpus-per-task=8      # cpus per subtask
#SBATCH --mem=10GB              # Needed memory per array-subtask 
#SBATCH --time=12:00:00         # time per subtask

set -euo pipefail

# This script finds all subdirectories within a specified target directory,
# zips each one, and then backs up the resulting zip files to an rclone remote.

# --- Configuration & Argument Checks ---
module load rclone/1.68.1

# 1. Check for Target Directory
if [ -z "$1" ]; then
    echo "Usage: $0 <target_directory> <rclone_remote_name> <rclone_bucket_path>"
    echo "Example: $0 /home/user/data my_s3_remote my-backups/data_zips"
    exit 1
fi

TARGET_DIR="$1"
RCLONE_REMOTE="$2"
RCLONE_BUCKET="$3"

# 2. Check for Rclone Remote Name
if [ -z "$RCLONE_REMOTE" ]; then
    echo "Error: Please provide the rclone remote configuration name (e.g., 'my_s3_remote')."
    exit 1
fi

# 3. Check for Rclone Bucket/Path
if [ -z "$RCLONE_BUCKET" ]; then
    echo "Error: Please provide the destination bucket or path (e.g., 'my-backups/data_zips')."
    exit 1
fi

# 4. Validate Target Directory
if [ ! -d "$TARGET_DIR" ]; then
    echo "Error: The path '$TARGET_DIR' is not a valid directory."
    exit 1
fi

# 5. Check if rclone is installed
if ! command -v rclone &> /dev/null
then
    echo "Error: 'rclone' is not installed or not in the PATH."
    echo "Please install rclone to enable the cloud backup functionality."
    exit 1
fi

echo "--- Starting Backup Process ---"
echo "Target Directory: $TARGET_DIR"
echo "Rclone Destination: $RCLONE_REMOTE:$RCLONE_BUCKET"
echo "-------------------------------"

# 6. Use 'find' to locate directories and pipe them to a loop
find "$TARGET_DIR" -mindepth 1 -type d -print0 | while IFS= read -r -d $'\0' folder; do
    
    FOLDER_BASENAME=$(basename "$folder")
    ZIP_FILE_PATH="${FOLDER_BASENAME}.zip"

    echo "1/2 Zipping folder: '$folder' to '$ZIP_FILE_PATH'..."

    # Create the zip archive
    # -r: recursive, -q: quiet (avoids excessive output)
    zip -r -q "$ZIP_FILE_PATH" "$folder"

    if [ $? -eq 0 ]; then
        echo "--> Zip Success: '$ZIP_FILE_PATH' created."

        # 7. Backup the compressed file using rclone
        # rclone copy automatically creates the destination bucket/path if it doesn't exist.
        echo "2/2 Backing up '$ZIP_FILE_PATH' to cloud storage..."
        
        # Use --progress and --stats 1s for visual feedback during large transfers
        rclone copy "$ZIP_FILE_PATH" "$RCLONE_REMOTE:$RCLONE_BUCKET" --progress --stats 1s

        if [ $? -eq 0 ]; then
            echo "--> Backup Success: '$ZIP_FILE_PATH' successfully uploaded."
            
            # 8. Clean up the local zip file after a successful backup
            echo "Cleaning up local file: '$ZIP_FILE_PATH'"
            rm -f "$ZIP_FILE_PATH"
        else
            echo "--> ERROR: rclone backup failed for '$ZIP_FILE_PATH'."
            echo "The local zip file has been left for manual inspection."
        fi
        
    else
        echo "--> ZIP ERROR: Failed to create zip for '$folder'. Skipping backup."
    fi
    
    echo "-------------------------------" # Separator for next iteration

done

echo "--- Script execution complete. ---"
