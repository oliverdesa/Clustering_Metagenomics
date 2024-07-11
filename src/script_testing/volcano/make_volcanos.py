import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
import os

def load_data(file_path):
    # Load the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep='\t')
    return df

def parse_data(df):
    # Extract enzyme type from the 'feature' column
    def extract_enzyme_type(feature):
        parts = feature.split('.')
        if len(parts) == 2:
            return parts[0]
        elif len(parts) == 3:
            return '.'.join(parts[:2])
        else:
            return parts[0]
    
    df['enzyme_type'] = df['feature'].apply(extract_enzyme_type)
    return df

def create_volcano_plot(df, enzyme_type, output_dir):
    # Filter data for the given enzyme type
    df_enzyme = df[df['enzyme_type'] == enzyme_type]
    
    # Create a volcano plot
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=df_enzyme, x='coef', y=-np.log10(df_enzyme['qval']),
                    edgecolor='k', s=100)
    plt.axhline(y=-np.log10(0.05), linestyle='--', color='red')
    plt.axvline(x=0, linestyle='--', color='grey')
    plt.xlabel('Coefficient')
    plt.ylabel('-log10(q-value)')
    plt.title(f'DA of {enzyme_type} in Dysbiotic Communities')
    plt.savefig(os.path.join(output_dir, f'volcano_{enzyme_type}.png'), dpi=600)
    plt.close()

def generate_volcano_plots(file_path, output_dir):
    # Load and parse data
    df = load_data(file_path)
    df = parse_data(df)
    
    # Get unique enzyme types
    enzyme_types = df['enzyme_type'].unique()
    
    # Generate volcano plot for each enzyme type
    for enzyme_type in enzyme_types:
        create_volcano_plot(df, enzyme_type, output_dir)

def main():
    parser = argparse.ArgumentParser(description="Generate volcano plots for each enzyme type")
    parser.add_argument('file_path', type=str, help="Path to the TSV file containing the data")
    parser.add_argument('--output_dir', type=str, help="Path to the directory where the volcano plots will be saved", default='.')
    args = parser.parse_args()
    
    # Create output directory if it does not exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    generate_volcano_plots(args.file_path, args.output_dir)

if __name__ == "__main__":
    main()
