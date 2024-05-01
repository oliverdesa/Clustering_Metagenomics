library(Maaslin2)
library(feather)
library(arrow)
library(tidyverse)
library(readr)

# To validate cluster effect sizes, run maaslin on individual genes and extract
# coefs from the genes making up one of our significant clusters, then plot
# These coefs compared to the overall cluster.

# Read in the feather
metadata_8156 <- read.csv("/Users/odesa/Library/CloudStorage/OneDrive-UniversityofToronto/LabWork/CRC/LatestDataJan/DRA008156/DRA008156_metadata_complete.tsv",
                          sep = "\t", header = TRUE)

# Read in the abundance data
abundance_data_8156 <- read.csv("/Users/odesa/Library/CloudStorage/OneDrive-UniversityofToronto/LabWork/CRC/LatestDataJan/DRA008156/clean_joined_genefamilies_relab_DRA008156.tsv",
                                sep = "\t", header = TRUE)

# set sample_id as rownames
abundance_data_8156 <- abundance_data_8156 %>% column_to_rownames(var = "sample_id")
metadata_8156 <- metadata_8156 %>% column_to_rownames(var = "Accession")

# Patterns to search for
patterns <- c("^DL.endopeptidase")

# Combine patterns into a single regular expression
pattern <- paste(patterns, collapse = "|")

selected_columns <- grep(pattern, names(abundance_data_8156), value = TRUE)

filtered_df_8156 <- abundance_data_8156[, selected_columns]

# Run Maaslin2
maaslin2_result <- Maaslin2(
  input_data = filtered_df_8156,
  input_metadata = metadata_8156,
  output = "/Users/odesa/Library/CloudStorage/OneDrive-UniversityofToronto/LabWork/CRC/maaslin/DRA008156/raw",
  fixed_effects = c("Group"),
  random_effects = NULL,
  transform = "none",
  reference=c("Group,Healthy")
)
