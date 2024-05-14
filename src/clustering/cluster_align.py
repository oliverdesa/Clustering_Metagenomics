#!/usr/bin/env python3

import pandas as pd
import argparse

def get_cluster_sequences(cluster_ids):
    """When given a cluster ID, map it to the mapping file,
       return all UniRef100 IDs in the cluster, and output their sequences as a FASTA file."""
    
    # Read in the clustering data
    mapping_info = pd.read_csv('/Users/odesa/Documents/CRC-Final/src/clustering/plotting_script/data/merged_info.tsv', sep='\t')

    # Initialize an empty dictionary
    uniref_dict = {}

    # Initialize variables to keep track of the current UniRef ID and sequence
    current_uniref_id = None
    current_sequence = []

    # Read the FASTA file and create a dictionary with UniRef ID as the key and sequence as the value
    with open('/Users/odesa/Documents/CRC-Final/data/clustering/dl_endos/dl_endopeptidases_fastas.faa', 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):  # If it's a header line
                # Save the previous sequence under its UniRef ID
                if current_uniref_id and current_sequence:
                    uniref_dict[current_uniref_id] = ''.join(current_sequence)
                    current_sequence = []
                # Extract and clean up the UniRef ID
                current_uniref_id = line.split('_')[2].split('|')[0]
            else:  # If it's a sequence line
                current_sequence.append(line)

        # Save the last sequence under its UniRef ID
        if current_uniref_id and current_sequence:
            uniref_dict[current_uniref_id] = ''.join(current_sequence)

    # Process each cluster ID
    for cluster_id in cluster_ids:
        # Get the UniRef100 IDs for the cluster
        clusters = mapping_info[mapping_info['dl_endopeptidase-foldseek_cluster'] == cluster_id]
        uniref_ids = clusters['Uniref'].values

        if cluster_id not in uniref_ids:
            uniref_ids = list(uniref_ids)
            uniref_ids.append(cluster_id)

        print(f"found cluster IDs for {cluster_id}: {uniref_ids}")

        # Create a filtered dictionary containing only the UniRef IDs present in the list
        filtered_dict = {key: value for key, value in uniref_dict.items() if key in uniref_ids}

        # Write the filtered UniRef IDs and sequences to a new FASTA file
        output_fasta_path = f'/Users/odesa/Desktop/transient_files/down_fastas/{cluster_id}_dl_endopeptidases_fastas.faa'
        with open(output_fasta_path, 'w') as f_out:
            for uniref_id, sequence in filtered_dict.items():
                # Write in FASTA format: header starting with '>', followed by the sequence
                f_out.write(f">{uniref_id}\n")
                # Wrap long sequences for better readability
                for i in range(0, len(sequence), 60):
                    f_out.write(f"{sequence[i:i+60]}\n")

        print(f"Filtered FASTA written to {output_fasta_path}")


# Main entry point to parse command-line arguments
if __name__ == '__main__':
    # Initialize argparse
    parser = argparse.ArgumentParser(description='Retrieve sequences for specified cluster IDs and output them in FASTA format.')
    parser.add_argument('cluster_ids', nargs='+', help='A list of cluster IDs to process, separated by spaces.')

    # Parse the arguments from the command line
    args = parser.parse_args()

    # Call the function with the provided list of cluster IDs
    get_cluster_sequences(args.cluster_ids)


