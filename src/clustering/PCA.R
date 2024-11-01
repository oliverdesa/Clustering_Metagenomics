library(Maaslin2)
library(feather)
library(arrow)
library(tidyverse)
library(vegan)

####################################################################
# PCA script
####################################################################

# Read in the feather
metadata_8156 <- arrow::read_feather("C:\\Users\\odesa\\OneDrive - University of Toronto\\LabWork\\CRC\\LatestDataJan\\DRA008156\\DRA008156_metadata_complete.feather")

# Read in the abundance data
abundance_data_8156 <- read.csv("C:\\Users\\odesa\\Desktop\\Code\\CRC-Final\\data\\clustering\\humann_clustered\\clustered_complete_DRA008156.tsv", 
                                sep = "\t", header = TRUE)

abundance_data_8156 <- abundance_data_8156 %>% column_to_rownames(var = "sample_id")

# Patterns to search for
patterns <- c("^DL.endopeptidase")

# Combine patterns into a single regular expression
pattern <- paste(patterns, collapse = "|")

selected_columns <- grep(pattern, names(abundance_data_8156), value = TRUE)

filtered_df_8156 <- abundance_data_8156[, selected_columns]

# convert rownames to a column
filtered_df_8156$sample_id <- rownames(filtered_df_8156)

merged_df <- merge(
  x = filtered_df_8156,
  y = metadata_8156[, c("Accession", "Group")],  # Select only the desired columns from df2
  by.x = "sample_id",                          # Specify the identifier column in df1
  by.y = "Accession",                   # Specify the identifier column in df2
  all.x = TRUE                           # Optional: to perform a left join
)


# Separate metadata and gene abundance data
metadata <- merged_df[, "Group"]
gene_abundance <- merged_df[, !names(merged_df) %in% "Group"]

# drop sample_id column
gene_abundance <- subset(gene_abundance, select = -sample_id)

# Perform PCA (scaling recommended)
pca <- prcomp(gene_abundance, scale = TRUE)

# Summarize PCA results
summary(pca)

# Visualize the PCA (e.g., first two PCs)
library(ggplot2)
pca_df <- data.frame(pca$x, State = metadata)
ggplot(pca_df, aes(x = PC1, y = PC2, color = State)) +
  geom_point(size = 3) +
  labs(title = "PCA of Gene Abundance Data", x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal()

# Check loadings (contributions of each gene to the PCs)
loadings <- as.data.frame(pca$rotation)

# Identify the top contributing genes for each principal component
# For example, for PC1 (top 10 genes)
top_genes_pc1 <- loadings[order(abs(loadings$PC1), decreasing = TRUE), "PC1"]
top_10_genes_pc1 <- rownames(loadings)[order(abs(loadings$PC1), decreasing = TRUE)][1:10]

# Print the top genes contributing to PC1
print(top_10_genes_pc1)

# Check for rows that are completely empty or contain only zeros
empty_rows <- rowSums(gene_abundance, na.rm = TRUE) == 0

# Identify rows with any missing values
rows_with_na <- apply(gene_abundance, 1, function(x) any(is.na(x)))

# Combine both conditions to exclude problematic rows
rows_to_exclude <- empty_rows | rows_with_na

# Display the problematic rows
print(which(rows_to_exclude))

# Exclude rows with all zeros or NA values
gene_abundance_filtered <- gene_abundance[!rows_to_exclude, ]

# Example Bray-Curtis dissimilarity matrix
# Recalculate the Bray-Curtis dissimilarity matrix after removing problematic rows
distance_matrix <- vegdist(gene_abundance_filtered, method = "bray")

metadata_filtered <- metadata[!rows_to_exclude, ]


# Run PERMANOVA using 'State' as the grouping variable
perm_results <- adonis2(distance_matrix ~ group, data = metadata)

# Display results
print(perm_results)
