library(Maaslin2)
library(dplyr)
library(ggbiplot)
library(feather)
library(arrow)
library(tidyverse)
library(ggplot2)

pooled_abun_df <- function(abun_df) {
  
  classes <- unlist(lapply(colnames(abun_df), function(col_name) {
    strsplit(col_name, "_")[[1]][1]
  }))
  
  pooled_abun_df <- aggregate(t(abun_df), by = list(classes), FUN = sum)
  
  rownames(pooled_abun_df) <- pooled_abun_df$Group.1
  pooled_abun_df$Group.1 <- NULL
  
  pooled_abun_df <- data.frame(t(pooled_abun_df))
  
  colnames(pooled_abun_df) <- unlist(lapply(colnames(pooled_abun_df), function(col_name) {
    sub(".", "-", col_name, fixed = TRUE)
  }))
  
  pooled_abun_df
  
}

## Load Data ##

# Read in the feather
metadata <- arrow::read_feather("C:/Users/odesa/Desktop/CRCFinal/PRJEB10878/PRJEB10878_metadata.feather")

# Read in the abundance data
abundance_data <- arrow::read_feather("C:/Users/odesa/Desktop/CRCFinal/PRJEB10878/clean_joined_genefamilies_relab.feather")

# read in the DAC abundance data
DAC_abundance_data <- arrow::read_feather("C:/Users/odesa/Desktop/CRCFinal/PRJEB10878/clean_joined_DAC_genefamilies_relab.feather")

## Process data ##

# Set the sample_id column as the index
metadata <- metadata %>% column_to_rownames(var = "Run")
abundance_data <- abundance_data %>% column_to_rownames(var = "sample_id")
DAC_abundance_data <- DAC_abundance_data %>% column_to_rownames(var = "sample_id")

# Try pooling data
pooled_abun_df <- pooled_abun_df(abundance_data)

DAC_pooled_abun_df <- pooled_abun_df(DAC_abundance_data)

# Drop the UNMAPPED column
pooled_abun_df <- subset(pooled_abun_df, select = -c(UNMAPPED))

DAC_pooled_abun_df <- subset(DAC_pooled_abun_df, select = -c(UNMAPPED))


## Run Maaslin2, CLR Normalization ##

results <- Maaslin2(
  input_data = pooled_abun_df,
  input_metadata = metadata,
  output = "C:/users/odesa/Desktop/CRCFinal/CRC-Final/data/PRJEB10878_PGH_maaslin",
  fixed_effects = c("config"),
  random_effects = NULL,
  normalization = "none", 
  transform = "none"
)

# DAC Analysis #

DAC_abundance_data <- arrow::read_feather("C:/Users/odesa/Desktop/CRCFinal/PRJEB10878/clean_joined_DAC_genefamilies_relab.feather")

DAC_abundance_data <- DAC_abundance_data %>% column_to_rownames(var = "sample_id")

DAC_pooled_abun_df <- pooled_abun_df(DAC_abundance_data)

DAC_pooled_abun_df <- subset(DAC_pooled_abun_df, select = -c(UNMAPPED))

DAC_pooled_abun_df$sample_id <- rownames(DAC_pooled_abun_df)

metadata$sample_id <- rownames(metadata)

# add the config column from metadata to DAC_pooled_abun_df on the same sample_id
DAC_pooled_abun_df <- merge(DAC_pooled_abun_df, metadata, by = "sample_id")


DAC_plot <- ggplot(DAC_pooled_abun_df, aes(x = config, y = Diadenylate, color = config)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1.5, aes(color = config)) +
  theme_minimal() +
  labs(x = "Patient", y = "DAC Relative Abundance")

ggsave("C:/users/odesa/Desktop/CRCFinal/CRC-Final/figures/DAC_abundance.png", plot = DAC_plot, width = 10, height = 10, dpi=300)
