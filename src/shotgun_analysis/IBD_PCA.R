library(vegan)
library(ggplot2)
library(dplyr)

# Load data
# Read in the feather
metadata_ibd <- read.csv("/Volumes/PGH-Backup/ibd_data/hmp2_metagenomics_metadata.csv")

# Read in the abundance data
abundance_data_ibd <- read.csv("/Volumes/PGH-Backup/ibd_data/humann_second_run/ibd_genefamilies_relab_clustered.tsv", 
                               sep = "\t", header = TRUE)

# Drop the column named 'X'
if ("X" %in% colnames(abundance_data_ibd)) {
  abundance_data_ibd <- abundance_data_ibd %>% select(-X)
}

# Find columns with non-numeric values
non_numeric_columns <- abundance_data_ibd %>%
  select_if(~ any(!is.numeric(.)))

# Print the non-numeric columns
print(non_numeric_columns)


distance_matrix <- vegdist(abundance_data_ibd, method = "bray")

# Perform PCoA
pcoa_result <- cmdscale(distance_matrix, eig = TRUE, k = 2)