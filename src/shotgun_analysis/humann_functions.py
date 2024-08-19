import pandas as pd
import xml.etree.ElementTree as ET

def clean_table(tsv):
    """Clean the humann table for downstream analysis"""

    #read to tsv
    df = pd.read_csv(tsv, sep='\t')

    # drops dupes including taxa info
    df = df[~df['# Gene Family'].apply(lambda x: '|' in x)]
    
    # Remove everything from the headers following the id, then transpose
    # df.columns = [col.split('_')[0] for col in df.columns]
    # transposed_df = df.T

    # Edit specifically for bariatric data
    
    df.columns = ['_'.join(col.split('_')[:2]) for col in df.columns]
    transposed_df = df.T


    # make the first row the header
    transposed_df.columns = transposed_df.iloc[0]
    transposed_df = transposed_df.drop(transposed_df.index[0])

    # make the index a column
    transposed_df = transposed_df.reset_index()
    transposed_df = transposed_df.rename(columns={'index': 'sample_id'})

    return transposed_df


def make_outputs(df, output_name):
    """Make a feather file and tsv from a dataframe"""


    df.to_csv(output_name + '.tsv', sep='\t')

    df.to_feather(output_name + '.feather')

    return output_name

def xml_to_spreadsheet_sample(xml_file_path, output_file_path):
    # Parse the XML file
    tree = ET.parse(xml_file_path)
    root = tree.getroot()

    # Create a list to hold the data
    data = []

    # Iterate over each SAMPLE in the XML
    for sample in root.findall('SAMPLE'):
        # Initialize a dictionary for this sample
        sample_info = {
            'Alias': sample.get('alias'),
            'Center Name': sample.get('center_name'),
            'Accession': sample.get('accession'),
            'BioSample ID': sample.find('.//PRIMARY_ID').text if sample.find('.//PRIMARY_ID') is not None else None,
            'Taxon ID': sample.find('.//TAXON_ID').text if sample.find('.//TAXON_ID') is not None else None,
            'Scientific Name': sample.find('.//SCIENTIFIC_NAME').text if sample.find('.//SCIENTIFIC_NAME') is not None else None
        }

        # Extract SAMPLE_ATTRIBUTES
        for attr in sample.findall('.//SAMPLE_ATTRIBUTE'):
            tag = attr.find('TAG').text if attr.find('TAG') is not None else None
            value = attr.find('VALUE').text if attr.find('VALUE') is not None else None
            if tag and value:
                sample_info[tag] = value

        # Append to data list
        data.append(sample_info)

    # Create a DataFrame from the data
    df = pd.DataFrame(data)

    # Export DataFrame to a CSV file
    df.to_csv(output_file_path, index=False)

def xml_to_spreadsheet_run(xml_file_path, output_file_path):
    # Parse the XML file
    tree = ET.parse(xml_file_path)
    root = tree.getroot()

    # Create a list to hold the data
    data = []

    # Iterate over each STUDY in the XML
    for run in root.findall('RUN'):
        # Initialize a dictionary for this study
        run_info = {
            'Alias': run.get('alias'),
            'Center Name': run.get('center_name'),
            'Accession': run.get('accession'),
            'Title': run.find('.//TITLE').text if run.find('.//TITLE') is not None else None,
            
        }

        # Append to data list
        data.append(run_info)

    # Create a DataFrame from the data
    df = pd.DataFrame(data)

    # Export DataFrame to a CSV file
    df.to_csv(output_file_path, index=False)


def group_humann_table(humann_feather):
    """Group the humann table by enzymes and group all enzyme together"""

    # read in the humann table
    humann_df = pd.read_feather(humann_feather)

    headers = humann_df.columns.tolist()

    original_width = len(headers)

    headers_enzymes = [x.split('_')[0] for x in headers]  # Extract enzyme names from column headers

    humann_df.columns = headers_enzymes  # Rename columns with enzyme names

    # Sum all columns with the same enzyme
    grouped_df = humann_df.groupby(humann_df.columns, axis=1).sum()

    new_width = len(grouped_df.columns.tolist())

    print(f'Original width: {original_width}, Grouped width: {new_width}')

    return grouped_df

def cluster_humann_table_improved(humann_feather, cluster_tsv):
    """Cluster the humann table for each of the PGH enzymes and store cluster information."""
    
    # read in the humann table
    humann_df = pd.read_feather(humann_feather)

    # read in the clustering dataframes
    cluster_df = pd.read_csv(cluster_tsv, sep='\t', low_memory=False)

    # list of enzymes
    enzymes = ['DL-endopeptidase', 'LD-carboxypeptidase', 
               'LD-endopeptidase', 'Glucosaminidase',
               'DD-carboxypeptidase', 'DD-endopeptidase',
               'Amidase', 'Muramidase']
    
    extra_classes = ['Saga', 'UC118']

    clustered_df = pd.DataFrame()
    
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
