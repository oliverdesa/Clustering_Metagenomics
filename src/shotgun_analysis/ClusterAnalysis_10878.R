library(Maaslin2)
library(feather)
library(arrow)
library(tidyverse)

## Load Data ##

# Read in the feather
metadata_10878 <- read.csv("C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\LatestDataJan\\PRJEB10878\\PRJEB10878Metadata.txt")

# Read in the cluster abundance data
abundance_data_10878 <- arrow::read_feather("C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\LatestDataJan\\PRJEB10878\\clustered_complete_10878.feather")

# set sample_id as rownames
abundance_data_10878 <- abundance_data_10878 %>% column_to_rownames(var = "sample_id")
metadata_10878 <- metadata_10878 %>% column_to_rownames(var = "Run")


## Run Maaslin2 ##

results <- Maaslin2(
  input_data = abundance_data_10878,
  input_metadata = metadata_10878,
  output = "C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\maaslin\\PRJEB10878_cluster\\Maaslin2ResultsNoNorm",
  fixed_effects = c("config"),
  random_effects = NULL,
  normalization = "none", 
  transform = "none",
)

results <- Maaslin2(
  input_data = abundance_data_10878,
  input_metadata = metadata_10878,
  output = "C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\maaslin\\PRJEB10878_cluster\\Maaslin2ResultsCLR",
  fixed_effects = c("config"),
  random_effects = NULL,
  normalization = "CLR", 
  transform = "none",
)
