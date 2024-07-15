library(vegan)
library(ggplot2)
library(dplyr)
library(plotly)
library(pairwiseAdonis)

## Load data ##

# Read in the feather
metadata_ibd <- read.csv("/Volumes/PGH-Backup/ibd_data/metadata/hmp2_metagenomics_metadata.csv")

# Read in the abundance data
abundance_data_ibd <- read.csv("/Volumes/PGH-Backup/ibd_data/humann_second_run/ibd_genefamilies_relab_clustered.tsv", 
                               sep = "\t", header = TRUE)

######################

## Pre-process data ##

# Drop the column named 'X'
if ("X" %in% colnames(abundance_data_ibd)) {
  abundance_data_ibd <- abundance_data_ibd %>% select(-X)
}

labels <- abundance_data_ibd$sample_id

df_numeric <- abundance_data_ibd %>% select(-sample_id)

# Remove empty rows

empty_rows <- apply(df_numeric, 1, function(row) all(is.na(row) | row == 0))

df_numeric_clean <- df_numeric[!empty_rows, ]

any_empty_rows <- any(apply(df_numeric_clean, 1, function(row) all(is.na(row) | row == 0)))
print(any_empty_rows)

######################

## subset for genefamilies ##

# Extract unique gene families with the first two parts of the identifier if they exist
get_gene_family <- function(col_name) {
  parts <- strsplit(col_name, "\\.")[[1]]
  if (length(parts) >= 3) {
    return(paste(parts[1], parts[2], sep = "."))
  } else {
    return(parts[1])
  }
}

gene_families <- unique(sapply(colnames(df_numeric_clean), get_gene_family))

print(gene_families)

######################

# Define PCoA function to caluclate PCoA, plot result and PERMANOVA test ##

perform_pcoa <- function(gene_family) {
  # Subset the dataframe for the specific gene family
  subset_df <- df_numeric_clean %>% select(starts_with(gene_family))
  
  # Remove columns with all zeros
  non_zero_columns <- apply(subset_df, 2, function(col) any(col != 0))
  subset_df <- subset_df[, non_zero_columns]
  
  # Remove rows with all zeros
  non_zero_rows <- apply(subset_df, 1, function(row) any(row != 0))
  subset_df <- subset_df[non_zero_rows, ]
  
  # Check if the subset_df is empty or has only one column (not enough for distance matrix)
  if (ncol(subset_df) < 2 || nrow(subset_df) < 2) {
    print(paste("Not enough data for gene family:", gene_family))
    return(NULL)
  }
  
  print(paste("Subset data for", gene_family, ":"))
  print(head(subset_df))
  
  # Calculate the distance matrix
  distance_matrix <- vegdist(subset_df, method = "bray")
  print("Distance matrix calculated.")
  
  # Perform PCoA
  pcoa_result <- cmdscale(distance_matrix, eig = TRUE, k = 2)
  
  # Check if the PCoA result is valid
  if (is.null(pcoa_result$points)) {
    print(paste("PCoA failed for gene family:", gene_family))
    return(NULL)
  }
  
  print(paste("PCoA result for", gene_family, ":"))
  print(head(pcoa_result$points))
  
  # Calculate the percentage of variance explained
  eig_vals <- pcoa_result$eig
  variance_explained <- eig_vals / sum(eig_vals) * 100
  percent_variance <- round(variance_explained[1:2], 2)
  
  labels_clean <- labels[!empty_rows][non_zero_rows]
  pcoa_df <- data.frame(pcoa_result$points, label = labels_clean)
  colnames(pcoa_df) <- c("PCoA1", "PCoA2", "label")
  
  # Ensure column names match for merging
  colnames(metadata_ibd)[colnames(metadata_ibd) == "External.ID"] <- "label"
  
  # Merge PCoA results with metadata
  merged_df <- merge(pcoa_df, metadata_ibd, by = "label")
  
  # Plot the results
  pcoa_plot <- ggplot(merged_df, aes(x = PCoA1, y = PCoA2, color = diagnosis)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(
      title = paste("PCoA Plot (2D) -", gene_family), 
      x = paste0("PCoA1 (", percent_variance[1], "%)"),
      y = paste0("PCoA2 (", percent_variance[2], "%)")
    ) +
    scale_color_manual(values = c("nonIBD" = "blue", "UC" = "red", "CD" = "green"))
  
  ggsave(paste0("IBD_PCoA_plot_", gene_family, ".png"), plot = pcoa_plot, width = 10, height = 10, dpi = 600)
  
  # Run PERMANOVA test
  permanova_result <- adonis(distance_matrix ~ diagnosis, data = merged_df)
  print(permanova_result$aov.tab)
  
  # Run pairwise PERMANOVA tests
  pairwise_results <- pairwise.adonis(distance_matrix, merged_df$diagnosis)
  print(pairwise_results)
  
  return(list(pcoa_result = pcoa_result, merged_df = merged_df))
}

# Loop over each gene family and perform PCoA
results <- lapply(gene_families, perform_pcoa)

# Remove NULL results
results <- results[!sapply(results, is.null)]

# Feature selection for each gene family
for (gene_family in gene_families) {
  result <- results[[which(gene_families == gene_family)]]
  if (is.null(result)) next
  
  print(paste("Processing gene family:", gene_family))
  pcoa_result <- result$pcoa_result
  feature_scores <- pcoa_result$points
  
  subset_df <- df_numeric_clean %>% select(starts_with(gene_family))
  
  # Ensure that subset_df and feature_scores have the same number of rows
  subset_df <- subset_df[!empty_rows, ]
  
  print("Subset data frame:")
  print(head(subset_df))
  
  print("Feature scores:")
  print(head(feature_scores))
  
  contributions <- apply(subset_df, 2, function(x) cor(x, feature_scores[,1], use = "complete.obs"))
  contributions <- cbind(contributions, apply(subset_df, 2, function(x) cor(x, feature_scores[,2], use = "complete.obs")))
  
  colnames(contributions) <- c("PCoA1", "PCoA2")
  
  # ordered_contributions_PCoA1 <- contributions[order(abs(contributions[,1]), decreasing = TRUE),]
  # ordered_contributions_PCoA2 <- contributions[order(abs(contributions[,2]), decreasing = TRUE),]
  
  # print(paste("Top 10 features for PCoA1 in", gene_family, ":"))
  # print(ordered_contributions_PCoA1[1:10,])
  # 
  # print(paste("Top 10 features for PCoA2 in", gene_family, ":"))
  # print(ordered_contributions_PCoA2[1:10,])
}






















# Calculate the distance matrix
distance_matrix <- vegdist(df_numeric_clean, method = "bray")

# Perform PCoA
pcoa_result <- cmdscale(distance_matrix, eig = TRUE, k = 2)

# Calculate the percentage of variance explained
eig_vals <- pcoa_result$eig
variance_explained <- eig_vals / sum(eig_vals) * 100
percent_variance <- round(variance_explained[1:2], 2)  # Get the first two components

labels_clean <- labels[!empty_rows]
pcoa_df <- data.frame(pcoa_result$points, label = labels_clean)
colnames(pcoa_df) <- c("PCoA1", "PCoA2", "label")

# Ensure column names match for merging
colnames(metadata_ibd)[colnames(metadata_ibd) == "External.ID"] <- "label"

# Merge PCoA results with metadata
merged_df <- merge(pcoa_df, metadata_ibd, by = "label")

######################

## Subset and Perform PCoA for each gene family ##

# Extract unique gene families
gene_families <- unique(sub("(.+?)\\..+", "\\1", colnames(df_numeric_clean)))

print(gene_families))


## Visualize results ##

# 2D plot using ggplot2 with percentage of variance explained in axis labels
pcoa_plot <- ggplot(merged_df, aes(x = PCoA1, y = PCoA2, color = diagnosis)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "PCoA Plot (2D)", 
    x = paste0("PCoA1 (", percent_variance[1], "%)"),
    y = paste0("PCoA2 (", percent_variance[2], "%)")
  ) +
  scale_color_manual(values = c("nonIBD" = "blue", "UC" = "red", "CD" = "green"))

ggsave("IBD_PCoA_plot.png", plot = pcoa_plot, width = 10, height = 10, dpi = 600)

######################

## Statistical analysis ##

# Run PERMANOVA test
permanova_result <- adonis(distance_matrix ~ diagnosis, data = merged_df)

permanova_table <- permanova_result$aov.tab

# View the results
print(permanova_table)

# Run pairwise PERMANOVA tests
pairwise_results <- pairwise.adonis(distance_matrix, merged_df$diagnosis)

# View pairwise comparison results
print(pairwise_results)

######################

## Feature selection ##

# Extract the most important features from the PCoA
feature_scores <- pcoa_result$points

# Calculate the pearson correlation between each column and the scores on the
# first PCoA. complete.obs ensures only complete cases are used in the correlation.
contributions <- apply(df_numeric_clean, 2, function(x) cor(x, feature_scores[,1], use = "complete.obs"))

# Does the same but for PCoA2 and then combines the results into a single matrix
contributions <- cbind(contributions, apply(df_numeric_clean, 2, function(x) cor(x, feature_scores[,2], use = "complete.obs")))

colnames(contributions) <- c("PCoA1", "PCoA2")

# Order contributions by magnitude for PCoA1
ordered_contributions_PCoA1 <- contributions[order(abs(contributions[,1]), decreasing = TRUE),]
# Order contributions by magnitude for PCoA2
ordered_contributions_PCoA2 <- contributions[order(abs(contributions[,2]), decreasing = TRUE),]

# Print the most important features for PCoA1
print("Top 10 features for PCoA1:")
print(ordered_contributions_PCoA1[1:10,])

# Print the most important features for PCoA2
print("Top 10 features for PCoA2:")
print(ordered_contributions_PCoA2[1:10,])



