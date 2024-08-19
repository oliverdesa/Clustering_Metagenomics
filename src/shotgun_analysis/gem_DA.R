library(Maaslin2)
library(readxl)
library(readr)
library(tidyverse)

# Read in abundance data
abundance <- readr::read_tsv("/Volumes/PGH-Backup/gem/genefamilies/clean_all_shallow_gems_genefamilies_relab_clustered.tsv")

abundance['sample_id']

# Read in metadata
meta <- 