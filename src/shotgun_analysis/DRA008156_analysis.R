library(Maaslin2)
library(dplyr)
library(ggbiplot)
library(feather)
library(arrow)
library(tidyverse)
library(scales)
library(vegan)
library(pairwiseAdonis)


pooled_abun_df <- function(abun_df) {
  
  classes <- unlist(lapply(colnames(abun_df), function(col_name) {
    strsplit(col_name, "_")[[1]][1]
  }))
  
  pooled_abun_df <- aggregate(t(abun_df), by = list(classes), FUN = sum)
  
  rownames(pooled_abun_df) <- pooled_abun_df$Group.1
  pooled_abun_df$Group.1 <- NULL
  
  pooled_abun_df <- data.frame(t(pooled_abun_df))
  
  colnames(pooled_abun_df) <- unlist(lapply(colnames(pooled_abun_df), function(col_name) {
    sub(".", "-", col_name, fixed = TRUE)
  }))
  
  pooled_abun_df
  
}

## Load data ##

# Read in the abundance data
abundance_DRA008156 <- arrow::read_feather("C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\LatestDataJan\\DRA008156\\clean_joined_genefamilies_relab_DRA008156.feather")

# Read in the metadata data
metadata_DRA008156 <- read.csv("C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\LatestDataJan\\DRA008156\\DRA008156_metadata_complete.tsv", sep = "\t")

# set sample_id as rownames
abundance_DRA008156 <- abundance_DRA008156 %>% column_to_rownames(var = "sample_id")
metadata_DRA008156 <- metadata_DRA008156 %>% column_to_rownames(var = "Accession")

pooled_abun_df_DRA008156 <- pooled_abun_df(abundance_DRA008156)

enzymes = c('Amidase', 'DD-carboxypeptidase', 'DD-endopeptidase', 'DL-endopeptidase', 
            'Glucosaminidase', 'LD-carboxypeptidase', 'LD-endopeptidase', 'Muramidase')

filt_pooled_abun_df_DRA008156 <- pooled_abun_df_DRA008156[, enzymes]


results <- Maaslin2(
  input_data = pooled_abun_df_DRA008156,
  input_metadata = metadata_DRA008156,
  output = "C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\maaslin\\DRA008156\\2_groups\\pooled_NoNorm",
  fixed_effects = c("Group"),
  random_effects = NULL,
  normalization = "none", 
  transform = "none",
  reference=c("Group,Healthy")
)

results <- Maaslin2(
  input_data = pooled_abun_df_DRA008156,
  input_metadata = metadata_DRA008156,
  output = "C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\maaslin\\DRA008156\\2_groups\\pooled_CLR",
  fixed_effects = c("Group"),
  random_effects = NULL,
  normalization = "CLR", 
  transform = "none",
  reference=c("Group,Healthy")
)

# Clustering Results

clustered_abundance_DRA008156 <- arrow::read_feather("C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\LatestDataJan\\DRA008156\\clustered_complete_DRA008156.feather")

clustered_abundance_DRA008156 <- clustered_abundance_DRA008156 %>% column_to_rownames(var = "sample_id")

results <- Maaslin2(
  input_data = clustered_abundance_DRA008156,
  input_metadata = metadata_DRA008156,
  output = "C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\maaslin\\DRA008156\\2_groups\\clustered_NoNorm",
  fixed_effects = c("Group"),
  random_effects = NULL,
  normalization = "none", 
  transform = "none",
  reference=c("Group,Healthy")
)

results <- Maaslin2(
  input_data = clustered_abundance_DRA008156,
  input_metadata = metadata_DRA008156,
  output = "C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\maaslin\\DRA008156\\2_groups\\clustered_CLR",
  fixed_effects = c("Group"),
  random_effects = NULL,
  normalization = "CLR", 
  transform = "none",
  reference=c("Group,Healthy")
)



