library(Maaslin2)
library(feather)
library(arrow)
library(tidyverse)

## Load Data ##

# Read in the feather
metadata_10878 <- arrow::read_feather("C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\LatestDataJan\\PRJEB7774\\PRJEB7774_metadata.feather")

# Read in the abundance data
abundance_data_10878 <- arrow::read_feather("C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\LatestDataJan\\PRJEB10878\\clean_joined_genefamilies_relab_10878.feather")

# set sample_id as rownames
abundance_data_7774 <- abundance_data_7774 %>% column_to_rownames(var = "sample_id")
metadata_7774 <- metadata_7774 %>% column_to_rownames(var = "run_accession")