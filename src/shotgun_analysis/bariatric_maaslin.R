library(Maaslin2)
library(feather)
library(arrow)
library(tidyverse)

############ Bariatric Data ############

## Load Grouped Data ##

# Read in the feather
metadata_bariatric <- arrow::read_feather("E:/bariatric/bariatric_metadata.feather")

metadata_bariatric$unique_sample_id <- paste(metadata_bariatric$sample_id, metadata_bariatric$TimePoint, sep = "_")

metadata_bariatric <- metadata_bariatric %>% column_to_rownames(var = "unique_sample_id")

metadata_bariatric <- metadata_bariatric[,-3]


# Read in the abundance data
grouped_bariatric <- arrow::read_feather("E:/bariatric/grouped_bariatric.feather")

grouped_bariatric <- grouped_bariatric[,-14]

grouped_bariatric$sample_id <- metadata_bariatric$unique_sample_id

grouped_bariatric <- grouped_bariatric %>% column_to_rownames(var = "sample_id")

## Load Cluster Data ##

# Read in the feather
clustered_bariatric <- arrow::read_feather("E:/bariatric/clustered_complete_bariatric.feather")

clustered_bariatric$sample_id <- metadata_bariatric$unique_sample_id

clustered_bariatric <- clustered_bariatric %>% column_to_rownames(var = "sample_id")


## Run Maaslin2 ##

results <- Maaslin2(
  input_data = grouped_bariatric,
  input_metadata = metadata_bariatric,
  output = "E:/bariatric/maaslin/grouped_no_clr",
  fixed_effects = c("TimePoint"),
  normalization = "none", 
  transform = "none",
  reference = c("TimePoint,BL"),
  random_effects = c("sample_id")
)

results <- Maaslin2(
  input_data = grouped_bariatric,
  input_metadata = metadata_bariatric,
  output = "E:/bariatric/maaslin/grouped_clr",
  fixed_effects = c("TimePoint"),
  normalization = "CLR", 
  transform = "none",
  reference = c("TimePoint,BL"),
  random_effects = c("sample_id")
)
