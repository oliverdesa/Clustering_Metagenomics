import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def load_data(file_path):
    df = pd.read_csv(file_path, sep='\t')
    return df

def calculate_correlations(df, gene_columns, score_column):
    spearman_results = {}
    pearson_results = {}
    
    # Calculate correlations for each gene and store the results
    for gene in gene_columns:
        spearman_corr, spearman_p = spearmanr(df[gene], df[score_column])
        pearson_corr, pearson_p = pearsonr(df[gene], df[score_column])
        
        spearman_results[gene] = (spearman_corr, spearman_p)
        pearson_results[gene] = (pearson_corr, pearson_p)
    
    return spearman_results, pearson_results

def plot_correlations(correlation_results, method):
    genes = list(correlation_results.keys())
    corrs = [correlation_results[gene][0] for gene in genes]
    pvals = [correlation_results[gene][1] for gene in genes]
    
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=genes, y=corrs, size=-np.log10(pvals), legend=False, sizes=(20, 200))
    plt.axhline(0, linestyle='--', color='grey')
    plt.xlabel('Gene')
    plt.ylabel(f'{method} Correlation')
    plt.title(f'{method} Correlation of Gene Counts with Dysbiosis Scores')
    plt.xticks(rotation=90)
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Calculate correlations between gene counts and dysbiosis scores")
    parser.add_argument('file_path', type=str, help="Path to the TSV file containing the data")
    parser.add_argument('--score_column', type=str, default='dysbiosis_score', help="Column name for dysbiosis scores")
    args = parser.parse_args()
    
    # Load the data
    df = load_data(args.file_path)
    
    # Identify gene columns (assuming all columns except the score column are gene counts)
    gene_columns = [col for col in df.columns if col != args.score_column]
    
    # Calculate correlations
    spearman_results, pearson_results = calculate_correlations(df, gene_columns, args.score_column)
    
    # Plot the results
    plot_correlations(spearman_results, 'Spearman')
    plot_correlations(pearson_results, 'Pearson')

if __name__ == "__main__":
    main()