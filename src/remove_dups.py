import pandas as pd

# Load the .tsv file into a pandas DataFrame
filename = "merged_output.tsv"
df = pd.read_csv(filename, sep='\t', index_col=0)

# Add a new column with simplified 'Gene Family' names
df['Simplified Gene Family'] = df.index.to_series().apply(lambda x: x.split('|')[0])

# Drop duplicates based on the new column
df = df.drop_duplicates(subset='Simplified Gene Family', keep='first')

# Optional: Save the cleaned data back to a new .tsv file, without the 'Simplified Gene Family' column
df.drop(columns='Simplified Gene Family').to_csv("cleaned_data.tsv", sep='\t')
