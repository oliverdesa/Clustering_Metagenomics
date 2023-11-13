import pandas as pd
import glob
from pathlib import Path

accession = 'met4'

# Define the path to the directory containing the TSV files
path_to_tsv_files = '/home/oliver/Desktop/MET4/met4_humann/' 

# Use glob to get all the TSV file paths
tsv_files = glob.glob(f'{path_to_tsv_files}*.tsv')

# Read each TSV into a DataFrame and store them in a list
dfs = []
for file in tsv_files:
    df = pd.read_csv(file, sep='\t', index_col='# Gene Family')
    # Rename the abundance column to the patient identifier based on the file name
    df.columns = [Path(file).stem]
    dfs.append(df)

# Combine all DataFrames, filling missing values with 0
combined_df = pd.concat(dfs, axis=1).fillna(0)

# write combined df to a feather file
combined_df.reset_index().to_feather(f'{accession}_combined_data.feather')
