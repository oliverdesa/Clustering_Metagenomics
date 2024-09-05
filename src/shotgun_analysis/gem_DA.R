library(Maaslin2)
library(readxl)
library(readr)
library(tidyverse)

#### Read in Data ####

# Read in abundance data
abundance <- readr::read_tsv("/Volumes/PGH-Backup/gem/genefamilies_shallow/gem_clustered_mapped.tsv")

# visualize IDs
abundance['original_id']

# reformat IDs
abundance$original_id <- gsub("_Sdna", "", abundance$original_id)

# Read in metadata
meta <- readr::read_tsv("/Volumes/PGH-Backup/gem/metadata_maaslin_filtered.tsv")

meta['File name of sequencing (HC)']

########################################

# filter metadata
filtered_meta <- meta %>%
  select(`File name of sequencing (HC)`, `cluster`, `diagnosis (0/1)` )

filtered_meta <- filtered_meta %>% distinct(`File name of sequencing (HC)`, .keep_all = TRUE)

abundance_column <- "original_id"
meta_column <- "File name of sequencing (HC)"

########################################

# Find the intersection of the values in the two different columns
matching_samples <- intersect(filtered_meta[[meta_column]], abundance[[abundance_column]])

# Check for duplicates in the matching samples
duplicates_in_meta <- sum(duplicated(filtered_meta$`File name of sequencing (HC)`))
duplicates_in_abundance <- sum(duplicated(abundance$original_id))

print(paste("Duplicates in metadata:", duplicates_in_meta))
print(paste("Duplicates in abundance:", duplicates_in_abundance))

########################################

# Subset both tables to keep only the matching rows
filtered_abundance <- abundance[abundance[[abundance_column]] %in% matching_samples, ]
filtered_meta <- filtered_meta[filtered_meta[[meta_column]] %in% matching_samples, ]

filtered_abundance <- filtered_abundance[, -c(1)]

# Check the dimensions of both filtered tables
dim(filtered_abundance)
dim(filtered_meta)

########################################

#### Check for non-numeric columns ####
non_numeric <- filtered_abundance %>% select_if(~ !is.numeric(.))
print(colnames(non_numeric_columns))

filtered_abundance <- filtered_abundance %>% column_to_rownames(var= "original_id")
filtered_meta <- filtered_meta %>% column_to_rownames(var= "File name of sequencing (HC)")

########################################

#### Re-order to match order ####

# Ensure the rownames are aligned based on matching sample IDs
filtered_meta <- filtered_meta[rownames(filtered_meta) %in% rownames(filtered_abundance), ]
filtered_abundance <- filtered_abundance[rownames(filtered_abundance) %in% rownames(filtered_meta), ]

# Reorder filtered_meta rows to match the order of filtered_abundance
filtered_meta <- filtered_meta[match(rownames(filtered_abundance), rownames(filtered_meta)), ]

# Check that the rownames are now in the same order
all(rownames(filtered_abundance) == rownames(filtered_meta))

########################################

# Check column names
print(colnames(filtered_abundance))

dim(filtered_abundance)
dim(filtered_meta)

names(filtered_meta)[names(filtered_meta) == 'diagnosis (0/1)'] <- 'diagnosis'

filtered_meta <- filtered_meta %>%
  mutate(diagnosis = case_when(
    diagnosis == 0 ~ "healthy",
    diagnosis == 1 ~ "ibd",
    TRUE ~ as.character(diagnosis)
  ))

# Patterns to search for
patterns <- c("^DL.endopeptidase")

# Combine patterns into a single regular expression
pattern <- paste(patterns, collapse = "|")

selected_columns <- grep(pattern, names(filtered_abundance), value = TRUE)

filtered_abundance <- filtered_abundance[, selected_columns]

## Begin Maaslin Analysis ##
results <- Maaslin2(
  input_data = filtered_abundance,
  input_metadata = filtered_meta,
  output = "/Volumes/PGH-Backup/gem/maaslin2",
  fixed_effects = c("diagnosis"),
  random_effects = "none",
  normalization = "none", 
  transform = "none",
  reference = c("diagnosis,healthy")
)









