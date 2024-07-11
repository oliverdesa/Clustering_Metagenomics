import pandas as pd
import seaborn as sns
import numpy as np
import argparse
import os
from adjustText import adjust_text

import matplotlib.pyplot as plt

def load_data(file_path):
    # Load the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep='\t')
    # Drop rows with value 'X' in the 'feature' column
    df = df[df['feature'] != 'X']
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
    
    # Create a dictionary to map enzyme types to colors
    enzyme_colors = {enzyme_type: color for enzyme_type, color in zip(df['enzyme_type'].unique(), sns.color_palette())}
    
    return df, enzyme_colors

def create_volcano_plot(df, condition, output_dir, enzyme_colors, labels):
    # Filter data for the various diagnosis
    df_condition = df[df['value'] == condition]
    
    # Create a volcano plot
    plt.figure(figsize=(20, 12))
    sns.scatterplot(data=df_condition, x='coef', y=-np.log10(df_condition['qval']),
                    hue='enzyme_type', palette=enzyme_colors, edgecolor='k', s=100)
    plt.axhline(y=-np.log10(0.05), linestyle='--', color='red')
    plt.axvline(x=0, linestyle='--', color='grey')
    plt.xlabel('Coefficient')
    plt.ylabel('-log10(q-value)')
    plt.title(f'DA of PGH clusters in {condition} patients')
    
    if labels:
        # Add labels to points
        texts = []
        for i, row in df_condition.iterrows():
            if row['coef'] > 1e-05 or row['coef'] < -1e-05 or -np.log10(row['qval']) > 6.5:
                texts.append(plt.text(row['coef'], -np.log10(row['qval']), row['feature'], fontsize=8))
        
        # Adjust the position of labels to avoid overlap
        adjust_text(texts)
    
    plt.savefig(os.path.join(output_dir, f'volcano_{condition}.png'), dpi=600)
    plt.close()

def generate_volcano_plots(file_path, output_dir, labels):
    # Load data
    df = load_data(file_path)
    df, enzyme_colors = parse_data(df)
    
    # Get unique conditions
    diagnosis = df['value'].unique()
    
    # Generate volcano plot for each enzyme type
    for condition in diagnosis:
        create_volcano_plot(df, condition, output_dir, enzyme_colors, labels)

def main():
    parser = argparse.ArgumentParser(description="Generate volcano plots for all PGH clusters in each condition.")
    parser.add_argument('file_path', type=str, help="Path to the TSV file containing the data")
    parser.add_argument('--output_dir', type=str, help="Path to the directory where the volcano plots will be saved", default='.')
    parser.add_argument('--labels', type=str, help="Whether or not to label individual points", default=True, choices=['True', 'False'])

    args = parser.parse_args()
    
    # Create output directory if it does not exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    generate_volcano_plots(args.file_path, args.output_dir, args.labels)

if __name__ == "__main__":
    main()
