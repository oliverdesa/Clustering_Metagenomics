import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from adjustText import adjust_text

def make_merged_table(abundance_path, dysbiosis_path):
    abundance = pd.read_csv(abundance_path, sep='\t')
    dysbiosis = pd.read_csv(dysbiosis_path, sep='\t', names=['sample_id', 'dysbiosis_score', 'dysbiosis_category'])
    abundance_merged = abundance.merge(dysbiosis, left_on='sample_id', right_on='sample_id')
    abundance_merged.drop(['Unnamed: 0', 'dysbiosis_category'], axis=1, inplace=True)
    return abundance_merged

def calculate_correlations(df, gene_columns, score_column):
    spearman_results = {}
    pearson_results = {}
    
    # Calculate correlations for each gene and store the results
    for gene in gene_columns:
        # Ensure columns are numeric
        df[gene] = pd.to_numeric(df[gene], errors='coerce')
        df[score_column] = pd.to_numeric(df[score_column], errors='coerce')
        
        # Drop NaN values
        valid_rows = df[[gene, score_column]].dropna()
        
        if not valid_rows.empty:
            spearman_corr, spearman_p = spearmanr(valid_rows[gene], valid_rows[score_column])
            pearson_corr, pearson_p = pearsonr(valid_rows[gene], valid_rows[score_column])
            
            spearman_results[gene] = (spearman_corr, spearman_p)
            pearson_results[gene] = (pearson_corr, pearson_p)
    
    return spearman_results, pearson_results

def plot_correlations(correlation_results, method):
    # Filter results for significant correlations (p < 0.05)
    significant_results = {gene: (corr, pval) for gene, (corr, pval) in correlation_results.items() if pval < 0.05}
    
    if significant_results:
        genes = list(significant_results.keys())
        corrs = [significant_results[gene][0] for gene in genes]
        pvals = [significant_results[gene][1] for gene in genes]

        plt.figure(figsize=(20, 12))
        sns.scatterplot(x=range(len(genes)), y=corrs, size=-np.log10(pvals), legend=False, sizes=(20, 200))
        plt.axhline(0, linestyle='--', color='grey')
        plt.xlabel('Gene')
        plt.ylabel(f'{method} Correlation')
        plt.title(f'{method} Correlation of Gene Counts with Dysbiosis Scores')
        plt.xticks([], [])  # Remove x-axis labels

        # Identify top 10 positive and negative correlations
        sorted_indices = np.argsort(corrs)
        top_negative_indices = sorted_indices[:10]
        top_positive_indices = sorted_indices[-10:]

        # Annotate the top negative and positive correlations
        texts = []
        for idx in top_negative_indices:
            texts.append(plt.text(idx, corrs[idx], genes[idx], ha='center', va='center'))
        for idx in top_positive_indices:
            texts.append(plt.text(idx, corrs[idx], genes[idx], ha='center', va='center'))
        
        adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black'))

        plt.savefig(f'{method}_correlation_clusters_dysbiosis', dpi=600)
    else:
        print(f"No significant {method} correlations (p < 0.05) found.")

def main():
    parser = argparse.ArgumentParser(description="Calculate correlations between gene counts and dysbiosis scores")
    parser.add_argument('abundance_file_path', type=str, help="Path to the TSV file containing the data")
    parser.add_argument('dysbiosis_file_path', type=str, help="Path to the TSV file containing the dysbiosis scores")
    parser.add_argument('--score_column', type=str, default='dysbiosis_score', help="Column name for dysbiosis scores")
    args = parser.parse_args()
    
    # Load the data
    df = make_merged_table(args.abundance_file_path, args.dysbiosis_file_path)
    
    # Identify gene columns (assuming all columns except the score column are gene counts)
    gene_columns = [col for col in df.columns if col != args.score_column]
    
    # Calculate correlations
    spearman_results, pearson_results = calculate_correlations(df, gene_columns, args.score_column)
    
    # Plot the results
    plot_correlations(spearman_results, 'Spearman')
    plot_correlations(pearson_results, 'Pearson')

if __name__ == "__main__":
    main()