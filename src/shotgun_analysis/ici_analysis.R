library(Maaslin2)
library(feather)
library(arrow)
library(tidyverse)
library(readr)

############ Immune Checkpoint Inhibitor Therapy Data ############

############ 70966 and 43119 ############

## Load Data ##

# Read in the feather
metadata_baseline_survival <- read.csv("/Users/odesa/OneDrive - University of Toronto/LabWork/ICI/metadata/cleaned_mdat_baseline_survival.csv")
metadata_longitudinal <- read.csv("/Users/odesa/OneDrive - University of Toronto/LabWork/ICI/metadata/cleaned_mdat_longitudinal.csv")

# Read in the abundance data
abundance_data_clustered <- read_tsv("/Users/odesa/OneDrive - University of Toronto/LabWork/ICI/LatestData/clustered_complete_70966_43119.tsv")
abundance_data_grouped <- read_tsv("/Users/odesa/OneDrive - University of Toronto/LabWork/ICI/LatestData/grouped_70966_43119.tsv")

# reformat the accessions
# Splitting the string and keeping only the first element
abundance_data_clustered$sample_id <- sapply(strsplit(abundance_data_clustered$sample_id, "_"), `[`, 1)

abundance_data_grouped$sample <- sapply(strsplit(abundance_data_grouped$sample, "_"), `[`, 1)

# set sample_id as rownames
abundance_data_clustered <- abundance_data_clustered %>% column_to_rownames(var = "sample_id")
abundance_data_grouped <- abundance_data_grouped %>% column_to_rownames(var = "sample")

metadata_baseline_survival <- metadata_baseline_survival %>% column_to_rownames(var = "accession")
metadata_longitudinal <- metadata_longitudinal %>% column_to_rownames(var = "accession")

# Patterns to search for
patterns <- c("^DL-endopeptidase")

# Combine patterns into a single regular expression
pattern <- paste(patterns, collapse = "|")

selected_columns <- grep(pattern, names(abundance_data_clustered), value = TRUE)
abundance_data_clustered <- abundance_data_clustered[, selected_columns]

abundance_data_grouped <- abundance_data_grouped[, "DL-endopeptidase", drop = FALSE]

#### Start MaAsLin2 ####

## colitis associations ##

# probably need to cat all metadata for this

results <- Maaslin2(
  input_data = abundance_data_clustered,
  input_metadata = metadata_baseline_survival,
  output = "/Users/odesa/OneDrive - University of Toronto/LabWork/ICI/maaslin/clustered_survival/colitis",
  fixed_effects = c("colitis"),
  random_effects = NULL,
  normalization = "none", 
  transform = "none",
  reference = c("colitis,no")
) # nothing

results <- Maaslin2(
  input_data = abundance_data_grouped,
  input_metadata = metadata_baseline_survival,
  output = "/Users/odesa/OneDrive - University of Toronto/LabWork/ICI/maaslin/grouped_survival/colitis",
  fixed_effects = c("colitis"),
  random_effects = NULL,
  normalization = "none", 
  transform = "none",
  reference = c("colitis,no")
) # nothing

## PFS12 associations ##

results <- Maaslin2(
  input_data = abundance_data_clustered,
  input_metadata = metadata_baseline_survival,
  output = "/Users/odesa/OneDrive - University of Toronto/LabWork/ICI/maaslin/clustered_survival/PFS12",
  fixed_effects = c("PFS12"),
  random_effects = NULL,
  normalization = "none", 
  transform = "none",
  reference = c("PFS12,no")
) # nothing

results <- Maaslin2(
  input_data = abundance_data_grouped,
  input_metadata = metadata_baseline_survival,
  output = "/Users/odesa/OneDrive - University of Toronto/LabWork/ICI/maaslin/grouped_survival/PFS12",
  fixed_effects = c("PFS12"),
  random_effects = NULL,
  normalization = "none", 
  transform = "none",
  reference = c("PFS12,no")
) # nothing

## antibiotic associations ##

results <- Maaslin2(
  input_data = abundance_data_clustered,
  input_metadata = metadata_baseline_survival,
  output = "/Users/odesa/OneDrive - University of Toronto/LabWork/ICI/maaslin/clustered_survival/antibiotic",
  fixed_effects = c("antibiotics"),
  random_effects = NULL,
  normalization = "none", 
  transform = "none",
  reference = c("antibiotics,no")
) # nothing

results <- Maaslin2(
  input_data = abundance_data_grouped,
  input_metadata = metadata_baseline_survival,
  output = "/Users/odesa/OneDrive - University of Toronto/LabWork/ICI/maaslin/grouped_survival/antibiotic",
  fixed_effects = c("antibiotics"),
  random_effects = NULL,
  normalization = "none", 
  transform = "none",
  reference = c("antibiotics,no")
) # nothing

## Progression associations ##

results <- Maaslin2(
  input_data = abundance_data_clustered,
  input_metadata = metadata_baseline_survival,
  output = "/Users/odesa/OneDrive - University of Toronto/LabWork/ICI/maaslin/clustered_survival/progression",
  fixed_effects = c("Progression"),
  random_effects = NULL,
  normalization = "none", 
  transform = "none",
  reference = c("Progression,0")
) # nothing

results <- Maaslin2(
  input_data = abundance_data_grouped,
  input_metadata = metadata_baseline_survival,
  output = "/Users/odesa/OneDrive - University of Toronto/LabWork/ICI/maaslin/grouped_survival/progression",
  fixed_effects = c("Progression"),
  random_effects = NULL,
  normalization = "none", 
  transform = "none",
  reference = c("Progression,0")
) # nothing

## recist associations ##

results <- Maaslin2(
  input_data = abundance_data_clustered,
  input_metadata = metadata_baseline_survival,
  output = "/Users/odesa/OneDrive - University of Toronto/LabWork/ICI/maaslin/clustered_survival/recist",
  fixed_effects = c("recist"),
  random_effects = NULL,
  normalization = "none", 
  transform = "none",
  reference = c("recist,PD")
) # nothing

results <- Maaslin2(
  input_data = abundance_data_grouped,
  input_metadata = metadata_baseline_survival,
  output = "/Users/odesa/OneDrive - University of Toronto/LabWork/ICI/maaslin/grouped_survival/recist",
  fixed_effects = c("recist"),
  random_effects = NULL,
  normalization = "none", 
  transform = "none",
  reference = c("recist,PD")
) # nothing

#### Longitudinal Analyses ####

results <- Maaslin2(
  input_data = abundance_data_clustered,
  input_metadata = metadata_baseline_survival,
  output = "/Users/odesa/OneDrive - University of Toronto/LabWork/ICI/maaslin/clustered_survival/longi_response",
  fixed_effects = c("visit", "ORR"),
  random_effects = NULL,
  normalization = "none", 
  transform = "none",
  reference = c("visit,0")
) # This doesnt work, ned to subset by visit

# subset by visit

output_base_path <- "/Users/odesa/OneDrive - University of Toronto/LabWork/ICI/maaslin/clustered_survival/longi_response"

# Loop through each visit number
for(visit_number in c(0, 1, 2, 3)) {
  # Subset metadata for the current visit
  metadata_subset <- metadata_longitudinal[metadata_longitudinal$visit == visit_number, ]
  
  # Subset abundance data to match the rownames of the metadata subset
  abundance_subset <- abundance_data_clustered[rownames(abundance_data_clustered) %in% rownames(metadata_subset), ]
  
  # Define output directory for this subset
  output_directory <- paste0(output_base_path, "/visit_", visit_number)
  
  # Run MaAsLin2
  results <- Maaslin2(
    input_data = abundance_subset,
    input_metadata = metadata_subset,
    fixed_effects = c("ORR"),
    random_effects = NULL,
    output = output_directory,
    normalization = "NONE",
    transform = "NONE"
  ) # nothing... Unfortunate
  
}

output_base_path <- "/Users/odesa/OneDrive - University of Toronto/LabWork/ICI/maaslin/clustered_survival/longi_PFS"

# Loop through each visit number
for(visit_number in c(0, 1, 2, 3)) {
  # Subset metadata for the current visit
  metadata_subset <- metadata_longitudinal[metadata_longitudinal$visit == visit_number, ]
  
  # Subset abundance data to match the rownames of the metadata subset
  abundance_subset <- abundance_data_clustered[rownames(abundance_data_clustered) %in% rownames(metadata_subset), ]
  
  # Define output directory for this subset
  output_directory <- paste0(output_base_path, "/visit_", visit_number)
  
  # Run MaAsLin2
  results <- Maaslin2(
    input_data = abundance_subset,
    input_metadata = metadata_subset,
    fixed_effects = c("PFS12"),
    random_effects = NULL,
    output = output_directory,
    normalization = "NONE",
    transform = "NONE"
  ) # Basically nothing. 2 Cluster close to significance, q = 0.13 for both at visit 2
  
}

output_base_path <- "/Users/odesa/OneDrive - University of Toronto/LabWork/ICI/maaslin/clustered_survival/longi_recist"

# Loop through each visit number
for(visit_number in c(0, 1, 2, 3)) {
  # Subset metadata for the current visit
  metadata_subset <- metadata_longitudinal[metadata_longitudinal$visit == visit_number, ]
  
  # Subset abundance data to match the rownames of the metadata subset
  abundance_subset <- abundance_data_clustered[rownames(abundance_data_clustered) %in% rownames(metadata_subset), ]
  
  # Define output directory for this subset
  output_directory <- paste0(output_base_path, "/visit_", visit_number)
  
  # Run MaAsLin2
  results <- Maaslin2(
    input_data = abundance_subset,
    input_metadata = metadata_subset,
    fixed_effects = c("recist"),
    random_effects = NULL,
    output = output_directory,
    normalization = "NONE",
    transform = "NONE",
    reference = c("recist,SD")
  ) # Basically nothing. 2 Cluster close to significance, q = 0.13 for both at visit 2
  
}

output_base_path <- "/Users/odesa/OneDrive - University of Toronto/LabWork/ICI/maaslin/clustered_survival/longi_colitis"

# Loop through each visit number
for(visit_number in c(0, 1, 2, 3)) {
  # Subset metadata for the current visit
  metadata_subset <- metadata_longitudinal[metadata_longitudinal$visit == visit_number, ]
  
  # Subset abundance data to match the rownames of the metadata subset
  abundance_subset <- abundance_data_clustered[rownames(abundance_data_clustered) %in% rownames(metadata_subset), ]
  
  # Define output directory for this subset
  output_directory <- paste0(output_base_path, "/visit_", visit_number)
  
  # Run MaAsLin2
  results <- Maaslin2(
    input_data = abundance_subset,
    input_metadata = metadata_subset,
    fixed_effects = c("colitis"),
    random_effects = NULL,
    output = output_directory,
    normalization = "NONE",
    transform = "NONE"
  ) # Basically nothing. 2 Cluster close to significance, q = 0.13 for both at visit 2
  
}
