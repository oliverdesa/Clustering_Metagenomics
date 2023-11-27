import pandas as pd

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






