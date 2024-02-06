library(Maaslin2)
library(feather)
library(arrow)
library(tidyverse)


############ Colorectal cancer data ############

############ DRA008156 ############

## Load Data ##

# Read in the feather
metadata_8156 <- arrow::read_feather("C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\LatestDataJan\\DRA008156\\DRA008156_metadata_complete.feather")

# Read in the abundance data
abundance_data_8156 <- read.csv("C:\\Users\\odesa\\Desktop\\Code\\CRC-Final\\data\\clustering\\humann_clustered\\clustered_complete_DRA008156.tsv", 
                                sep = "\t", header = TRUE)

# set sample_id as rownames
abundance_data_8156 <- abundance_data_8156 %>% column_to_rownames(var = "sample_id")
metadata_8156 <- metadata_8156 %>% column_to_rownames(var = "Accession")

# Patterns to search for
patterns <- c("^DL.endopeptidase", "^LD.carboxypeptidase", "^DD.carboxypeptidase")

# Combine patterns into a single regular expression
pattern <- paste(patterns, collapse = "|")

selected_columns <- grep(pattern, names(abundance_data_8156), value = TRUE)

filtered_df_8156 <- abundance_data_8156[, selected_columns]

## Run Maaslin2 ##

results <- Maaslin2(
  input_data = filtered_df_8156,
  input_metadata = metadata_8156,
  output = "C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\maaslin\\DRA008156\\new_cluster\\no_clr",
  fixed_effects = c("Group"),
  random_effects = NULL,
  normalization = "none", 
  transform = "none",
  reference = c("Group,Healthy")
)

results <- Maaslin2(
  input_data = filtered_df_8156,
  input_metadata = metadata_8156,
  output = "C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\maaslin\\DRA008156\\new_cluster\\clr",
  fixed_effects = c("Group"),
  random_effects = NULL,
  normalization = "CLR", 
  transform = "none",
  reference = c("Group,Healthy")
)


############ PRJEB7774 ############

## Load Data ##

# Read in the feather
metadata_7774 <- arrow::read_feather("C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\LatestDataJan\\PRJEB7774\\PRJEB7774_metadata.feather")

# Read in the abundance data
abundance_data_7774 <- read.csv("C:\\Users\\odesa\\Desktop\\Code\\CRC-Final\\data\\clustering\\humann_clustered\\clustered_complete_PRJEB7774.tsv",
                                sep='\t', header = TRUE)

# set sample_id as rownames
abundance_data_7774 <- abundance_data_7774 %>% column_to_rownames(var = "sample_id")
metadata_7774 <- metadata_7774 %>% column_to_rownames(var = "run_accession")

# Patterns to search for
patterns <- c("^DL.endopeptidase", "^LD.carboxypeptidase", "^DD.carboxypeptidase")

# Combine patterns into a single regular expression
pattern <- paste(patterns, collapse = "|")

selected_columns <- grep(pattern, names(abundance_data_7774), value = TRUE)

filtered_df_7774 <- abundance_data_7774[, selected_columns]


## Run Maaslin2 ##

results <- Maaslin2(
  input_data = filtered_df_7774,
  input_metadata = metadata_7774,
  output = "C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\maaslin\\PRJEB7774\\new_cluster\\no_clr",
  fixed_effects = c("sample_title"),
  random_effects = NULL,
  normalization = "none", 
  transform = "none",
  reference=c("sample_title,Control")
)

results <- Maaslin2(
  input_data = filtered_df_7774,
  input_metadata = metadata_7774,
  output = "C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\maaslin\\PRJEB7774\\new_cluster\\clr",
  fixed_effects = c("sample_title"),
  random_effects = NULL,
  normalization = "CLR", 
  transform = "none",
  reference=c("sample_title,Control")
)


############ PRJEB10878 ############

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






############ ICI Trial Data ############

############ PRJEB22893 ############

## Load Data ##















