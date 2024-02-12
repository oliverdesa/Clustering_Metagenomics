library(Maaslin2)
library(feather)
library(arrow)
library(tidyverse)
library(MicrobiomeStat)

############ Bariatric Data ############

## Load Grouped Data ##

# Read in the feather
metadata_bariatric <- arrow::read_feather("E:/bariatric/bariatric_metadata.feather")

metadata_bariatric$unique_sample_id <- paste(metadata_bariatric$sample_id, metadata_bariatric$TimePoint, sep = "_")

metadata_bariatric <- metadata_bariatric %>% column_to_rownames(var = "unique_sample_id")

metadata_bariatric <- metadata_bariatric[,-3]

metadata_bariatric$TimePoint <- factor(metadata_bariatric$TimePoint, levels = c("BL", "OR", "1M", "6M"))

metadata_bariatric$TimePoint <- as.numeric(metadata_bariatric$TimePoint)

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
  reference = c("TimePoint,1"),
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


# transpose the grouped_bariatric

# grouped_bariatric <- t(grouped_bariatric)

#' Feature.tab <- grouped_bariatric
#' 
#' metadata_bariatric$unique_id <- row.names(metadata_bariatric)
#' 
#' Meta.dat <- metadata_bariatric
#' 
#' Meta.dat$TimePoint <- relevel(as.factor(Meta.dat$TimePoint), ref="1")
#' 
#' MicrobiomeData <- list(
#'   feature.tab = Feature.tab,
#'   meta.dat = Meta.dat
#' )
#' 
#' test.list <- generate_taxa_trend_test_long(
#'   data.obj = MicrobiomeData,
#'   subject.var = "unique_id",
#'   time.var = "TimePoint",
#'   #group.var
#'   feature.level = "original",
#'   feature.dat.type = "proportion"
#' )
#' 
#' #' Generate the volcano plot
#' plot.list <- generate_taxa_volcano_single(
#'   data.obj = MicrobiomeData,
#'   test.list = test.list,
#'   feature.sig.level = 0.1,
#'   feature.mt.method = "none"
#' )
