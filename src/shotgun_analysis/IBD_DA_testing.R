library(Maaslin2)
library(feather)
library(arrow)
library(tidyverse)


############ Inflammatory Bowel Disease Data ############

## Load Data ##

# Read in the feather
metadata_ibd <- read.csv("/Volumes/PGH-Backup/ibd_data/hmp2_metadata_2018-08-20.csv")

# Read in the abundance data
abundance_data_ibd <- read.csv("/Volumes/PGH-Backup/ibd_data/humann_second_run/ibd_genefamilies_relab_clustered.tsv", 
                                sep = "\t", header = TRUE)

abundance_data_ibd_unique <- abundance_data_ibd[!duplicated(abundance_data_ibd$sample_id), ]

# Temp remove rownames
abundance_data_ibd_unique <- abundance_data_ibd_unique %>% 
  rownames_to_column(var = "temp_rowname") %>% 
  select(-temp_rowname)

# Then, set the 'sample_id' column as row names
abundance_data_ibd_unique <- abundance_data_ibd_unique %>% 
  column_to_rownames(var = "sample_id")

# Now filter metadata to only include rows present in abundance df
abundance_sample_ids <- rownames(abundance_data_ibd_unique)

# Also filter to metagenomic data
filtered_metadata_ibd <- metadata_ibd %>%
  filter(`External.ID` %in% abundance_sample_ids & data_type == "metagenomics")

filtered_metadata_ibd <- filtered_metadata_ibd %>% column_to_rownames(var = "External.ID")

# Patterns to search for
patterns <- c("^DL.endopeptidase")

# Combine patterns into a single regular expression
pattern <- paste(patterns, collapse = "|")

selected_columns <- grep(pattern, names(abundance_data_ibd_unique), value = TRUE)

abundance_data_ibd_unique <- abundance_data_ibd_unique[, selected_columns]

## Run Maaslin2 ##

results <- Maaslin2(
  input_data = abundance_data_ibd_unique,
  input_metadata = filtered_metadata_ibd,
  output = "/Volumes/PGH-Backup/ibd_data/maaslin2/maaslin2_results_ibd",
  fixed_effects = c("diagnosis"),
  random_effects = NULL,
  normalization = "none", 
  transform = "none",
  reference = c("diagnosis,nonIBD")
)






