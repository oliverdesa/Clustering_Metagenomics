#!/bin/bash

## Download PDB files for each cluster

# Directory containing the input files
input_dir="/media/oliver/PGH_Backup/clustering"

# Loop over each file in the input directory
for file in "$input_dir"/*; do
    # Extract the base name of the file
    base_name=$(basename "$file" .txt)

    # Create a directory with the base name
    mkdir -p "/media/oliver/PGH_Backup/clustering/$base_name"

    # Loop over lines in the file
    while IFS= read -r line; do
        # Construct the download URL
        url="https://alphafold.ebi.ac.uk/files/AF-${line}-F1-model_v4.pdb"

        # Download the file and save it in the directory named after the base name
        curl -o "${base_name}/AF-${line}-F1-model_v4.pdb" "$url"
    done < "$file"
done


