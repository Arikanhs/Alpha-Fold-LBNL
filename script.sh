#!/bin/bash

# List of protein codes
protein_codes=("1L2Y" "1UBQ" "2GB1" "1CTF" "4HHB" "1HRC" "3KUY" "1YCR" "2W5I" "2LZM" "2BG9")

# Loop through each protein code and call the Python functions
for code in "${protein_codes[@]}"
do
    echo "Running for PDB code: $code"
    python manipulation_function.py "$code"
    echo " "
done