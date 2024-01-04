import pandas as pd
import xml.etree.ElementTree as ET

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



