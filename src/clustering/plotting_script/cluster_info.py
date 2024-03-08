#!/usr/bin/env python3
"""
Author : oliv123 <oliv123@localhost>
Date   : 2024-02-16
Purpose: a script to return cluster information from a given cluster ID
"""

import os
import sys
import glob
import argparse
from typing import NamedTuple, TextIO
import pandas as pd
from pathlib import Path
import warnings
from datetime import datetime
from collections import Counter
import ast
import matplotlib.pyplot as plt
import seaborn as sns



class Args(NamedTuple):
    """ Command-line arguments """
    positional: str

# --------------------------------------------------
    
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Script to return cluster information for a given cluster ID or IDs from a .txt file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('cluster_input',
                        help='Cluster ID or a path to a .txt file with one ID on each line',
                        type=str)

    parser.add_argument('-t','--table',
                        help='return a summary table with information',
                        action='store_true')

    return parser.parse_args()

def process_cluster_input(cluster_input):
    """Process the cluster_input argument to determine if it's a file or a single ID"""
    cluster_file_path = Path(cluster_input)
    if cluster_file_path.is_file():
        with open(cluster_file_path, 'r') as file:
            cluster_ids = [line.strip() for line in file]
    else:
        cluster_ids = [cluster_input]
    
    return cluster_ids

# --------------------------------------------------

def pie_chart(taxa_info: list, output_dir: Path):
    
    # Get taxa info for pie charts
    phyla = []
    family = []
    genus = []

    # format taxa strings to select taxanomic level info
    for tax in taxa_info:
        phyla.append(tax.split(';')[1].replace('p__', ''))
        family.append((tax.split(';')[4]).replace('f__', ''))
        genus.append(tax.split(';')[-2].replace('g__', ''))

    # get the counts of individual taxa
    phyla_counts = Counter(phyla)
    family_counts = Counter(family)
    genus_counts = Counter(genus)

    # Create pie charts
    sns.set_style('whitegrid')

    # Create a figure to hold the subplots
    plt.figure(figsize=(24, 12))

    # Phyla Distribution
    plt.subplot(1, 3, 1)  # 1 row, 3 columns, 1st subplot
    plt.pie(phyla_counts.values(), labels=phyla_counts.keys(), autopct='%1.1f%%')
    plt.title('Phyla Distribution')

    # Family Distribution
    plt.subplot(1, 3, 2)  # 1 row, 3 columns, 2nd subplot
    plt.pie(family_counts.values(), labels=family_counts.keys(), autopct='%1.1f%%')
    plt.title('Family Distribution')

    # Genus Distribution
    plt.subplot(1, 3, 3)  # 1 row, 3 columns, 3rd subplot
    plt.pie(genus_counts.values(), labels=genus_counts.keys(), autopct='%1.1f%%')
    plt.title('Genus Distribution')

    # Show the figure with the subplots
    plt.tight_layout()
    plt.savefig(output_dir / 'pie_chart.png')
    plt.close()


# --------------------------------------------------

def bar_chart(domain_info: list, secretion_info: list, number_proteins: int, output_dir: Path):

    # Get counts of domains and secretion tags
    domain_counter = Counter(domain_info)
    secretion_counter = Counter(secretion_info)

    # Calculate average domain/protein
    domain_average = {key: (value / number_proteins) for key, value in domain_counter.items()}

    # Calculate percent secretion

    secretion_percent = {key: (value / number_proteins) * 100 for key, value in secretion_counter.items()}

    # Preparing data for plotting
    domain_labels = list(domain_average.keys())
    domain_percentages = list(domain_average.values())

    secretion_labels = list(secretion_percent.keys())
    secretion_percentages = list(secretion_percent.values())


    with warnings.catch_warnings():
        # plotting warning regarding palette and hue can be ignored for now
        warnings.simplefilter("ignore", FutureWarning)

        # Set a larger font size for all plot elements
        sns.set_context("talk")

        # Set style
        sns.set_style('whitegrid')

        # Create a larger figure to accommodate the subplots
        plt.figure(figsize=(24, 12))

        # Domain Distribution
        plt.subplot(1, 2, 1)  # 1 row, 2 columns, 1st subplot
        domain_plot = sns.barplot(x=domain_labels, y=domain_percentages, palette="Blues_d", edgecolor='black', linewidth=2)
        plt.title('Average # Domains per Protein', fontsize=20)
        plt.xlabel('Domains', fontsize=18)
        plt.ylabel('Domain per Protein', fontsize=18)  # Adding label for clarity
        plt.xticks(rotation=45, fontsize=14)
        plt.yticks(fontsize=14)

        # Secretion Distribution
        plt.subplot(1, 2, 2)  # 1 row, 2 columns, 2nd subplot
        secretion_plot = sns.barplot(x=secretion_labels, y=secretion_percentages, palette="Greens_d", edgecolor='black', linewidth=2)
        plt.title('Secretion Distribution', fontsize=20)
        plt.xlabel('Secretion Systems', fontsize=18)
        plt.ylabel('Percentage of Cluster Members')  # Intentionally left blank as per your setup
        plt.xticks(rotation=45, fontsize=14)
        plt.yticks(fontsize=14)

        # Adjust the layout
        plt.tight_layout()

        # Show the figure with the subplots
        plt.tight_layout()
        plt.savefig(output_dir / 'bar_chart.png')
        plt.close()

# --------------------------------------------------

def describe_cluster(cluster_id: str, protein_info_table: Path, cluster_info_table: Path, base_output_dir: Path):
    """
    This function takes a cluster ID, a table of information about the clusters and the cluster maps,
    returning pie charts of the taxa information along with bar plots of domains per protein and secretion
    percentage of all proteins in the cluster. 
    """

    # Create output dirs for the plotting functions
    output_dir = base_output_dir / cluster_id
    output_dir.mkdir(parents=True, exist_ok=True)

    # Read in tables as dfs
    cluster_info_table = pd.read_csv(cluster_info_table, sep='\t')
    protein_info_table = pd.read_csv(protein_info_table, sep='\t')

    # given the cluster_id, find the initial uniref IDs and make them into a list
    clustered_dl_list = list(cluster_info_table.loc[cluster_info_table['dl_endopeptidase-foldseek_cluster'] == cluster_id, 'Uniref'])

    number_proteins = len(list(protein_info_table.loc[protein_info_table['Uniref'].isin(clustered_dl_list), '# ID']))

    # search the protein_info_table for the reuturned Uniref ids. return unique species
    species_list = list(set(protein_info_table.loc[protein_info_table['Uniref'].isin(clustered_dl_list), 'Lineage']))

    # Find percentage of proteins that are secreted
    secretion_list = list(protein_info_table.loc[protein_info_table['Uniref'].isin(clustered_dl_list), 'Prediction'])

    # Get domain information and normalize to per protein

    # select the domain info column and make it into a list
    domain_info = list(protein_info_table.loc[protein_info_table['Uniref'].isin(clustered_dl_list), 'Interpro'])

    # flatten the list
    domain_info = [item for sublist in domain_info for item in (sublist if isinstance(sublist, list) else [sublist])]

    # remove any NA values
    domain_info = [item for item in domain_info if not pd.isna(item)]

    # convert strings to lists
    domain_info = [ast.literal_eval(item) for item in domain_info]

    # flatten the list
    domain_info = [item for sublist in domain_info for item in (sublist if isinstance(sublist, list) else [sublist])]

    print(f"Number of proteins in cluster {cluster_id}: {number_proteins}")
    print(f"number unique species in cluster {cluster_id}: {len(species_list)}")



    pie_chart(species_list, output_dir)
    bar_chart(domain_info, secretion_list, number_proteins, output_dir)


# --------------------------------------------------
def main() -> None:
    """ Main function """
    args = get_args()
    if args.table:
        print('table')
    
    cluster_ids = process_cluster_input(args.cluster_input)

    # Create a base output directory for all clusters, e.g., based on the current datetime
    base_output_dir = Path(f'output_{datetime.now().strftime("%Y%m%d_%H%M%S")}')
    base_output_dir.mkdir(parents=True, exist_ok=True)

    protein_info_table = Path('data/protein_info.tsv')
    cluster_info_table = Path('data/cluster_info.tsv')

    for cluster_id in cluster_ids:
        describe_cluster(cluster_id, protein_info_table, cluster_info_table, base_output_dir)


# --------------------------------------------------
if __name__ == '__main__':
    main()