import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

## Functions used for plotting and analysis ##

def scale_values(dataframe):
    # Scale the values in the dataframe
    # Assuming that the columns 'feature' and 'value' exist in the dataframe
    
    dataframe['scaled_value'] = -np.log(dataframe['qval']) * np.sign(dataframe['coef'])
    
    dataframe = dataframe.sort_values(by='scaled_value')

    return dataframe


def create_heatmap(dataframe, title, filename, target_column):
    # Pivot the DataFrame for the heatmap
    # Assuming that the columns 'feature', 'value', and 'scaled_value' exist in the dataframe
    heatmap_data = dataframe.pivot(index='feature', columns='value', values='scaled_value')

    # Sort the pivot data by one of the columns for visual hierarchy in the heatmap
    sorted_heatmap_data = heatmap_data.sort_values(by=target_column, ascending=False)

    # Create the heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(sorted_heatmap_data, annot=True, cmap='coolwarm_r', fmt='.2f', center = 0)

    # Customize the plot
    plt.title(title)
    plt.ylabel('DLE Clusters')
    plt.xlabel('')

    # Save and show the plot
    plt.tight_layout() # Ensure everything fits without overlapping
    plt.savefig(filename, dpi=600, bbox_inches='tight')
    plt.show()


## Functions used for clustering ##
    
# Create a function, when given a foldseek ID, returns all proteins within that cluster

def get_proteins(cluster_id, foldseek_dict, mmseqs_dict):
    ''' Input a cluster id and 2 dictionaries describing cluster
        patterns, return a list of all proteins in that cluster'''
    
    all_proteins = []

    mmseqs_rep = foldseek_dict.get(cluster_id, [])

    # For each mmseqs rep, add the proteins in its cluster to the list
    for id in mmseqs_rep:
        proteins = mmseqs_dict.get(id, [])
        all_proteins.extend(proteins)
    
    return all_proteins


def get_cluster(raw_list, mmseqs_dict, foldseek_dict):
    ''' Input a list of ids and 2 dictionaries describing cluster
        patterns, return a list of foldseek cluster ids'''
    
    mmseqs_list = []
    hits = []

    for value in raw_list:
        for key, content in mmseqs_dict.items():
            if value in content or value == key:
                mmseqs_list.append(key)
                hits.append(value)

    # Need to figure out what to do here if mmseqs dont map
    no_match = [x for x in raw_list if x not in hits]

    # print(f'length of no hits: {len(no_match)}')

    foldseek_list = []
    hits = []

    for value in mmseqs_list:

        for key, content in foldseek_dict.items():
            if value in content or value == key:
                foldseek_list.append(key)
                hits.append(value)

    no_match = [x for x in mmseqs_list if x not in hits]

    foldseek_list.extend(no_match)

    return foldseek_list

def get_cluster_info(ID_list, mmseqs_dict, foldseek_dict, sec_dict):
    '''Input a list of foldseek IDs and create a dataframe containing
       the number of proteins included in each cluster as well as the 
       number of secreted proteins within each cluster'''
    
    new_df = pd.DataFrame(columns=['cluster_id', 'proteins', 'secreted'])

    for idx, id in enumerate(ID_list):
        proteins = get_proteins(id, foldseek_dict, mmseqs_dict)

        sec_count = 0
        for protein in proteins:
            if sec_dict[protein] == 'SP':
                sec_count += 1
        
        new_df.loc[idx] = [id, len(proteins), sec_count]
    
    return new_df

def make_outputs(df, output_name):
    """Make a feather file and tsv from a dataframe"""

    df.to_csv(output_name + '.tsv', sep='\t')

    df.to_feather(output_name + '.feather')

    return output_name
