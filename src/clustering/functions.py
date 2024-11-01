import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import networkx as nx
from collections import defaultdict
import itertools
from sklearn.metrics.pairwise import cosine_similarity
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

## Functions used for plotting and analysis ##

def scale_values(dataframe):
    # Scale the values in the dataframe
    
    dataframe['scaled_value'] = -np.log(dataframe['qval']) * np.sign(dataframe['coef'])
    
    dataframe = dataframe.sort_values(by='scaled_value')

    return dataframe


def create_heatmap(dataframe, title, filename):
    # Pivot the DataFrame for the heatmap
    # Assuming that the columns 'feature', 'value', and 'scaled_value' exist in the dataframe
    heatmap_data = dataframe.pivot(index='feature', columns='value', values='scaled_value')

    # Sort the pivot data by one of the columns for visual hierarchy in the heatmap
    # sorted_heatmap_data = heatmap_data.sort_values(by=target_column, ascending=False)

    # Create the heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(heatmap_data, annot=True, cmap='coolwarm_r', fmt='.2f', center = 0)

    # Customize the plot
    plt.title(title)
    plt.ylabel('DLE Clusters')
    plt.xlabel('')

    # Save and show the plot
    plt.tight_layout() 
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

## Data Processing ##

def make_outputs(df, output_name):
    """Make a feather file and tsv from a dataframe"""

    df.to_csv(output_name + '.tsv', sep='\t')

    df.to_feather(output_name + '.feather')

    return output_name

def clean_table(tsv):
    """Clean the humann table for downstream analysis"""

    #read to tsv
    df = pd.read_csv(tsv, sep='\t')

    # drops dupes including taxa info
    df = df[~df['# Gene Family'].apply(lambda x: '|' in x)]
    
    # Remove everything from the headers following the id, then transpose
    df.columns = [col.split('_')[0] for col in df.columns]
    transposed_df = df.T
    
  
    # make the first row the header
    transposed_df.columns = transposed_df.iloc[0]
    transposed_df = transposed_df.drop(transposed_df.index[0])

    # make the index a column
    transposed_df = transposed_df.reset_index()
    transposed_df = transposed_df.rename(columns={'index': 'sample_id'})

    return transposed_df


def create_maps(mmseqs_tsv, foldseek_tsv, output_name):
    """Create a dictionary mapping uniref100s to mmseqs
       and foldseek clusters. Input clustering file, and
       return a dictionary mapping uniref100s to foldseek"""

    # read in the tsvs
    mmseqs_df = pd.read_csv(mmseqs_tsv, sep='\t')
    foldseek_df = pd.read_csv(foldseek_tsv, sep='\t')

    # rename the columns
    mm_column_names = [f'{output_name}-mmseqs_cluster', f'{output_name}-unclustered']
    fold_column_names = [f'{output_name}-foldseek_cluster', f'{output_name}-mmseqs_unclustered']
    mmseqs_df.columns = mm_column_names
    foldseek_df.columns = fold_column_names

    # parse the DFs keeping only Uniref IDs
    mmseqs_df[f'{output_name}-mmseqs_cluster'] = mmseqs_df[f'{output_name}-mmseqs_cluster'].apply(lambda x: x.split('_')[1])
    mmseqs_df[f'{output_name}-unclustered'] = mmseqs_df[f'{output_name}-unclustered'].apply(lambda x: x.split('_')[1])

    foldseek_df[f'{output_name}-foldseek_cluster'] = foldseek_df[f'{output_name}-foldseek_cluster'].apply(lambda x: x.split('-')[1])
    foldseek_df[f'{output_name}-mmseqs_unclustered'] = foldseek_df[f'{output_name}-mmseqs_unclustered'].apply(lambda x: x.split('-')[1])

    # Merging on 'mmseqs_cluster' from mmseqs_df and 'unclustered' from foldseek_df
    combined_df = mmseqs_df.merge(foldseek_df, left_on=f'{output_name}-mmseqs_cluster', right_on=f'{output_name}-mmseqs_unclustered', how='left')

    # subset to only important columns
    combined_df = combined_df[[f'{output_name}-unclustered', f'{output_name}-mmseqs_cluster', f'{output_name}-foldseek_cluster']]

    combined_df.to_csv(output_name + '.tsv', sep='\t')


def cluster_humann_table(humann_feather, cluster_tsv):
    """Cluster the humann table for each of the PGH enzymes
       input the raw humann df, clustering dataframes, and
       output the clustered humann df"""
    
    # read in the humann table
    humann_df = pd.read_feather(humann_feather)

    # read in the clustering dataframes
    cluster_df = pd.read_csv(cluster_tsv, sep='\t', low_memory=False)

    # list of enzymes
    enzymes = ['DL-endopeptidase', 'LD-carboxypeptidase', 
               'LD-endopeptidase', 'Glucosaminidase',
               'DD-carboxypeptidase', 'DD-endopeptidase',
               'Amidase', 'Muramidase']

    clustered_df = pd.DataFrame()
    for enzyme in enzymes:
        df = humann_df.loc[:, humann_df.columns.str.startswith(enzyme)]
        column_names = df.columns.tolist()

        print(f'{len(column_names)} {enzyme} found')

        # Extract the UniRef IDs from the column names
        column_ids = [x.split('_')[2] for x in column_names]

        # Get the foldseek cluster for each UniRef ID
        results = []
        unclustered = []
        for id in column_ids:
            result = cluster_df.loc[cluster_df[f"{enzyme.replace('-', '_').lower()}-unclustered"] == id, f"{enzyme.replace('-', '_').lower()}-foldseek_cluster"]
            if pd.isna(result).all():
                unclustered.append(f"{id}")
                results.append("unclustered")
            if  not pd.isna(result).all():
                results.append(enzyme + '-' + str(result.iloc[0]))
        
        print(f"{len(results)} {enzyme} found, {len(unclustered)} {enzyme} unclustered")

        # Replace the column names with the foldseek cluster
        df.columns = results

        # Aggregate the columns by foldseek cluster
        agg_df = df.T.groupby(df.columns).sum().T

        # Add the aggregated df to the clustered df
        clustered_df = pd.concat([clustered_df, agg_df], axis=1)
    
    # Add the sample id column back to the dataframe
    clustered_df['sample_id'] = humann_df['sample_id']

    return clustered_df

def cluster_humann_table_improved(humann_feather, cluster_tsv):
    """Cluster the humann table for each of the PGH enzymes and store cluster information."""
    
    # read in the humann table
    humann_df = pd.read_csv(humann_feather, sep='\t')

    # read in the clustering dataframes
    cluster_df = pd.read_csv(cluster_tsv, sep='\t', low_memory=False)

    # list of enzymes
    enzymes = ['DL-endopeptidase', 'LD-carboxypeptidase', 
               'LD-endopeptidase', 'Glucosaminidase',
               'DD-carboxypeptidase', 'DD-endopeptidase',
               'Amidase', 'Muramidase']
    
    extra_classes = ['Saga', 'UC118']

    clustered_df = pd.DataFrame()
    
    # This will store information about each cluster
    cluster_info_list = []

    # Create a mapping for each enzyme's clusters beforehand
    cluster_map = {}
    for enzyme in enzymes:
        enzyme_col = f"{enzyme.replace('-', '_').lower()}-unclustered"
        cluster_col = f"{enzyme.replace('-', '_').lower()}-foldseek_cluster"
        enzyme_cluster_map = cluster_df.set_index(enzyme_col)[cluster_col].to_dict()
        cluster_map[enzyme] = enzyme_cluster_map

    # Process each enzyme
    for enzyme in enzymes:
        df = humann_df.loc[:, humann_df.columns.str.startswith(enzyme)]
        column_names = df.columns.tolist()

        print(f'{len(column_names)} {enzyme} found')

        # Extract the UniRef IDs from the column names
        column_ids = [x.split('_')[2] for x in column_names]

        # Get the foldseek cluster for each UniRef ID
        results = []
        clusters_info = {}
        for id in column_ids:
            result = cluster_map[enzyme].get(id, "unclustered")
            if result != "unclustered":
                cluster_id = f"{enzyme}-{result}"
                clusters_info.setdefault(cluster_id, []).append(id)
            results.append(f"{enzyme}-{result}" if result != "unclustered" else "unclustered")
        
        # Replace the column names with the foldseek cluster
        df.columns = results

        # Aggregate the columns by foldseek cluster
        agg_df = df.T.groupby(df.columns).sum().T

        # Add the aggregated df to the clustered df
        clustered_df = pd.concat([clustered_df, agg_df], axis=1)

        # Collect the cluster information for analysis
        for cluster_id, ids in clusters_info.items():
            # Sum the final abundance for this cluster
            final_abundance = agg_df[cluster_id].sum()

            # Add the cluster info
            cluster_info_list.append({
                'cluster_id': cluster_id,
                'enzyme': enzyme,
                'num_uniref_ids': len(ids),
                'final_abundance': final_abundance
            })
    
    # Aggregate the extra classes (Saga and uc118) into single columns each
    for extra_class in extra_classes:
        df_extra = humann_df.loc[:, humann_df.columns.str.startswith(extra_class)]
        
        if not df_extra.empty:
            print(f'{len(df_extra.columns)} {extra_class} found')
            # Sum all columns for the extra class into one column
            extra_class_agg = df_extra.sum(axis=1)
            clustered_df[f'{extra_class}_aggregated'] = extra_class_agg

            # Collect the info for the extra classes
            cluster_info_list.append({
                'cluster_id': f'{extra_class}_aggregated',
                'enzyme': extra_class,
                'num_uniref_ids': df_extra.shape[1],
                'final_abundance': extra_class_agg.sum()
            })
        else:
            print(f'No {extra_class} found')

    # Add the sample id column back to the dataframe
    clustered_df['sample_id'] = humann_df['sample_id']
    
    # Convert cluster info list to DataFrame
    cluster_info_df = pd.DataFrame(cluster_info_list)
    
    return clustered_df, cluster_info_df


def group_humann_table(humann_table):
    """Group the humann table by enzymes and group all enzyme together"""

    if isinstance(humann_table, str):
        if humann_table.endswith('.tsv') or humann_table.endswith('.csv'):
            humann_df = pd.read_csv(humann_table, sep='\t')
        elif humann_table.endswith('.feather'):
            humann_df = pd.read_feather(humann_table)
        else:
            sys.exit('Invalid file type, must be tsv, csv, or feather')
    elif isinstance(humann_table, pd.DataFrame):
        humann_df = humann_table
    else:
        sys.exit('Invalid input type, must be a file path or a pandas DataFrame')

    headers = humann_df.columns.tolist()

    original_width = len(headers)

    headers_enzymes = [x.split('_')[0] for x in headers]  # Extract enzyme names from column headers

    humann_df.columns = headers_enzymes  # Rename columns with enzyme names

    # Sum all columns with the same enzyme
    grouped_df = humann_df.groupby(humann_df.columns, axis=1).sum()

    new_width = len(grouped_df.columns.tolist())

    print(f'Original width: {original_width}, Grouped width: {new_width}')

    return grouped_df


def plot_network(cluster_ids, significant_edges, association_table, merged_df, title):

    """ Plot the network graph with adjusted node colors and sizes based on association
        values and cluster sizes, respectively. Scale the node colour based on the
        association with the given feature. """
    
    # Make association dict
    association_dict = dict(zip(association_table['feature'], association_table['scaled_value']))

    # Create the network graph
    G_adjusted_similarity = nx.Graph()

    # Add nodes (clusters)
    G_adjusted_similarity.add_nodes_from(cluster_ids)

    # Add edges with weights based on cosine similarity
    G_adjusted_similarity.add_weighted_edges_from(significant_edges)

    # Define color scaling function
    def get_node_color(node):
        if node in association_dict:
            value = association_dict[node]
            return value
        else:
            return None

    # Assign colors to nodes based on their association values
    node_values = [get_node_color(node) for node in G_adjusted_similarity.nodes()]
    abs_max_val = max(abs(val) for val in node_values if val is not None)
    min_val = -abs_max_val
    max_val = abs_max_val
    reversed_bwr = plt.cm.bwr.reversed()  # Reverse the colormap

    node_colors = [
        'grey' if value is None or -3 < value < 3 else reversed_bwr((value - min_val) / (max_val - min_val))
        for value in node_values   
    ]

    # Calculate the size of each cluster
    cluster_sizes = merged_df['dl_endopeptidase-foldseek_cluster'].value_counts().to_dict()

    # Normalize cluster sizes to a suitable range for node sizes in the graph
    min_size = 20
    max_size = 1000
    min_cluster_size = min(cluster_sizes.values())
    max_cluster_size = max(cluster_sizes.values())
    node_sizes = [
        ((cluster_sizes[node] - min_cluster_size) / (max_cluster_size - min_cluster_size) * (max_size - min_size) + min_size)
        if node in cluster_sizes else min_size
        for node in G_adjusted_similarity.nodes()
    ]

    # Visualize the adjusted network
    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(G_adjusted_similarity, seed=42)  # For consistent layout

    # Draw the network
    nodes = nx.draw_networkx_nodes(G_adjusted_similarity, pos, node_color=node_colors, node_size=node_sizes, alpha=0.8, cmap=plt.cm.bwr)
    nx.draw_networkx_edges(G_adjusted_similarity, pos, alpha=0.5)
    nx.draw_networkx_labels(G_adjusted_similarity, pos, font_size=5, alpha=0.7)

    if title is not None:
        plt.title(title)
    else:
        plt.title("Adjusted Network Graph Based on Percentage Domain Inclusion")
        plt.axis('off')

    # Create a ScalarMappable for the colorbar
    sm = ScalarMappable(cmap=reversed_bwr, norm=Normalize(vmin=min_val, vmax=max_val))
    sm.set_array([])  # Dummy array for ScalarMappable
    plt.colorbar(sm, label='Association Value')

    plt.savefig('../../figures/CD_adjusted_network.png', dpi=600, bbox_inches='tight')
            




    

    

       