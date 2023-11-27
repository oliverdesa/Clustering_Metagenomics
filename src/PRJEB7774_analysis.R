library(Maaslin2)
library(dplyr)
library(ggbiplot)
library(feather)
library(arrow)
library(tidyverse)

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

## Load Data ##

# Read in the feather
metadata_7774 <- arrow::read_feather("C:/Users/odesa/Desktop/CRCFinal/PRJEB7774/PRJEB7774_metadata.feather")

# Read in the abundance data
abundance_data_7774 <- arrow::read_feather("C:/Users/odesa/Desktop/CRCFinal/PRJEB7774/clean_joined_genefamilies_relab_7774.feather")

# set sample_id as rownames
abundance_data_7774 <- abundance_data_7774 %>% column_to_rownames(var = "sample_id")
metadata_7774 <- metadata_7774 %>% column_to_rownames(var = "run_accession")

# pool the data
pooled_abun_df_7774 <- pooled_abun_df(abundance_data_7774)

# drop the UNMAPPED column
pooled_abun_df_7774 <- subset(pooled_abun_df_7774, select = -c(UNMAPPED))

## Run Maaslin2 ##

results <- Maaslin2(
  input_data = pooled_abun_df_7774,
  input_metadata = metadata_7774,
  output = "C:/users/odesa/Desktop/CRCFinal/CRC-Final/data/PRJEB7774_maaslin",
  fixed_effects = c("sample_title"),
  random_effects = NULL,
  normalization = "CLR", 
  transform = "none",
  reference=c("sample_title,Control")
)
