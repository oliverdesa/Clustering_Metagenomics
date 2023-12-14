#!/bin/bash

while IFS= read -r line
do
    # Construct the download URL
    url="https://alphafold.ebi.ac.uk/files/AF-${line}-F1-model_v4.pdb"

    # Download the file
    curl -o "AF-${line}-F1-model_v4.pdb" "$url"
done < "uniprot_ids.txt"

