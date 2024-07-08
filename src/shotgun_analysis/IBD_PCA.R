library(vegan)
library(ggplot2)
library(dplyr)
library(plotly)
library(pairwiseAdonis)

# Load data
# Read in the feather
metadata_ibd <- read.csv("/Volumes/PGH-Backup/ibd_data/metadata/hmp2_metagenomics_metadata.csv")

# Read in the abundance data
abundance_data_ibd <- read.csv("/Volumes/PGH-Backup/ibd_data/humann_second_run/ibd_genefamilies_relab_clustered.tsv", 
                               sep = "\t", header = TRUE)

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

# Calculate the distance matrix
distance_matrix <- vegdist(df_numeric_clean, method = "bray")

# Perform PCoA
pcoa_result <- cmdscale(distance_matrix, eig = TRUE, k = 2)

labels_clean <- labels[!empty_rows]
pcoa_df <- data.frame(pcoa_result$points, label = labels_clean)
colnames(pcoa_df) <- c("PCoA1", "PCoA2", "label")

# Ensure column names match for merging
colnames(metadata_ibd)[colnames(metadata_ibd) == "External.ID"] <- "label"

# Merge PCoA results with metadata
merged_df <- merge(pcoa_df, metadata_ibd, by = "label")

# 2D plot using ggplot2
ggplot(merged_df, aes(x = PCoA1, y = PCoA2, color = diagnosis)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCoA Plot (2D)", x = "PCoA1", y = "PCoA2") +
  scale_color_manual(values = c("nonIBD" = "blue", "UC" = "red", "CD" = "green"))

# Run PERMANOVA test
permanova_result <- adonis(distance_matrix ~ diagnosis, data = merged_df)

permanova_table <- permanova_result$aov.tab

# View the results
print(permanova_table)

# Run pairwise PERMANOVA tests
pairwise_results <- pairwise.adonis(distance_matrix, merged_df$diagnosis)

# View pairwise comparison results
print(pairwise_results)



