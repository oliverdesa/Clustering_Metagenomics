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
               'DD-carboxypeptidase', 'Diadenylate-cyclase',
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



            




    

    

       