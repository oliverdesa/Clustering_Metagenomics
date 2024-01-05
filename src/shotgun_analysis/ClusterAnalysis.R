library(Maaslin2)
library(feather)
library(arrow)
library(tidyverse)

## Load Data ##

# Read in the feather
metadata_7774 <- arrow::read_feather("C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\LatestDataJan\\PRJEB7774\\PRJEB7774_metadata.feather")

# Read in the abundance data
abundance_data_7774 <- arrow::read_feather("C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\LatestDataJan\\PRJEB7774\\clustered_complete_7774.feather")

# set sample_id as rownames
abundance_data_7774 <- abundance_data_7774 %>% column_to_rownames(var = "sample_id")
metadata_7774 <- metadata_7774 %>% column_to_rownames(var = "run_accession")


## Run Maaslin2 ##

results <- Maaslin2(
  input_data = abundance_data_7774,
  input_metadata = metadata_7774,
  output = "C:\\Users\\odesa\\Desktop\\Code\\CRC-Final\\data\\clustering\\Maaslin2Results",
  fixed_effects = c("sample_title"),
  random_effects = NULL,
  normalization = "none", 
  transform = "none",
  reference=c("sample_title,Control")
)

results <- Maaslin2(
  input_data = abundance_data_7774,
  input_metadata = metadata_7774,
  output = "C:\\Users\\odesa\\Desktop\\Code\\CRC-Final\\data\\clustering\\Maaslin2ResultsCLR",
  fixed_effects = c("sample_title"),
  random_effects = NULL,
  normalization = "CLR", 
  transform = "none",
  reference=c("sample_title,Control")
)
