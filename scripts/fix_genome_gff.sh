#!/bin/bash

# ---
# A script to fix GFF files where the strand column (column 7)
# contains a '?' instead of '+', '-', or '.'.
# It replaces the '?' with a '.'
# ---

# Check if the user provided exactly two arguments
if [ "$#" -ne 2 ]; then
    echo "Error: Invalid number of arguments."
    echo "Usage: $0 <input_file.gff> <output_file.gff>"
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="$2"

# Check if the input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found."
    exit 1
fi

# Inform the user what's happening
echo "Scanning '$INPUT_FILE' to fix strand column..."

# Run the awk command
#  -v OFS='\t' : Sets the Output Field Separator to a tab, same as input.
#  '$7 == "?"' : This is a condition. It checks if the 7th column is exactly "?".
#  '{ $7 = "." }': This is the action. If the condition is true, it sets the 7th column to ".".
#  '1'         : This is an awk shorthand that means "print the current line".
#                It ensures all lines (both modified and unmodified) are printed.
awk 'BEGIN{FS=OFS="\t"} $7 == "?" { $7 = "." } 1' "$INPUT_FILE" > "$OUTPUT_FILE"

# Check if the awk command succeeded
if [ "$?" -eq 0 ]; then
    echo "Success! Fixed file saved as '$OUTPUT_FILE'."
else
    echo "An error occurred during processing."
    exit 1
fi
