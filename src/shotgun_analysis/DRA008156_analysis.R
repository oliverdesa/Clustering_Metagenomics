library(Maaslin2)
library(dplyr)
library(ggbiplot)
library(feather)
library(arrow)
library(tidyverse)
library(scales)
library(vegan)
library(pairwiseAdonis)


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

## Load data ##

# Read in the abundance data
abundance_DRA008156 <- arrow::read_feather("C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\LatestDataJan\\DRA008156\\clean_joined_genefamilies_relab_DRA008156.feather")

# Read in the metadata data
metadata_DRA008156 <- read.csv("C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\LatestDataJan\\DRA008156\\DRA008156_metadata_complete.tsv", sep = "\t")

# set sample_id as rownames
abundance_DRA008156 <- abundance_DRA008156 %>% column_to_rownames(var = "sample_id")
metadata_DRA008156 <- metadata_DRA008156 %>% column_to_rownames(var = "Accession")

pooled_abun_df_DRA008156 <- pooled_abun_df(abundance_DRA008156)

enzymes = c('Amidase', 'DD-carboxypeptidase', 'DD-endopeptidase', 'DL-endopeptidase', 
            'Glucosaminidase', 'LD-carboxypeptidase', 'LD-endopeptidase', 'Muramidase')

filt_pooled_abun_df_DRA008156 <- pooled_abun_df_DRA008156[, enzymes]


results <- Maaslin2(
  input_data = pooled_abun_df_DRA008156,
  input_metadata = metadata_DRA008156,
  output = "C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\maaslin\\DRA008156\\2_groups\\pooled_NoNorm",
  fixed_effects = c("Group"),
  random_effects = NULL,
  normalization = "none", 
  transform = "none",
  reference=c("Group,Healthy")
)

results <- Maaslin2(
  input_data = pooled_abun_df_DRA008156,
  input_metadata = metadata_DRA008156,
  output = "C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\maaslin\\DRA008156\\2_groups\\pooled_CLR",
  fixed_effects = c("Group"),
  random_effects = NULL,
  normalization = "CLR", 
  transform = "none",
  reference=c("Group,Healthy")
)

# Clustering Results

clustered_abundance_DRA008156 <- arrow::read_feather("C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\LatestDataJan\\DRA008156\\clustered_complete_DRA008156.feather")

clustered_abundance_DRA008156 <- clustered_abundance_DRA008156 %>% column_to_rownames(var = "sample_id")

results <- Maaslin2(
  input_data = clustered_abundance_DRA008156,
  input_metadata = metadata_DRA008156,
  output = "C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\maaslin\\DRA008156\\2_groups\\clustered_NoNorm",
  fixed_effects = c("Group"),
  random_effects = NULL,
  normalization = "none", 
  transform = "none",
  reference=c("Group,Healthy")
)

results <- Maaslin2(
  input_data = clustered_abundance_DRA008156,
  input_metadata = metadata_DRA008156,
  output = "C:\\Users\\odesa\\OneDrive - University of Toronto\\CRC\\maaslin\\DRA008156\\2_groups\\clustered_CLR",
  fixed_effects = c("Group"),
  random_effects = NULL,
  normalization = "CLR", 
  transform = "none",
  reference=c("Group,Healthy")
)

# make indices of metadata_7774 into column
metadata_DRA008156$run_accession <- rownames(metadata_DRA008156)
filt_pooled_abun_df_DRA008156$run_accession <- rownames(filt_pooled_abun_df_DRA008156)

# Add sample_title to pooled_abun_df_7774 where run_accession matches
filt_pooled_abun_df_DRA008156$sample_title <- metadata_DRA008156$Group[match(filt_pooled_abun_df_DRA008156$run_accession, 
                                                                             metadata_DRA008156$run_accession)]

# convert to long df format
long_clr_pooled<- filt_pooled_abun_df_DRA008156 %>% 
  pivot_longer(
    cols = c(Amidase:Muramidase), names_to = "Enzyme", 
    values_to = "Relative Abundance")

# Plot the differential abundances and save to a tiff for pub

my_plot <- ggplot(long_clr_pooled, aes(x = Enzyme, y = `Relative Abundance`, color = sample_title)) +
  geom_point(alpha = 0.8, size = 2, stroke = 0.3, position = position_jitterdodge(jitter.height = 0,
                                                                                  jitter.width = 0.15)) +
  stat_summary(aes(group = sample_title), fun.min = function(x) {quantile(x, 0.25)}, fun.max = function(x) {quantile(x, 0.75)},
               size = 0.2, width = 0.3, geom = "errorbar", position = position_dodge(width = 0.8), color='black') +
  stat_summary(aes(group = sample_title), fun.min = median, fun.max = median,
               size = 0.2, width = 0.6, geom = "errorbar", position = position_dodge(width = 0.8), color='black') +
  scale_color_manual(values = c("Cancer" = "#cf5154", "Healthy" = "#4e88aa")) +
  labs(x = "Enzyme", y = "Relative Abundance") +
  theme_minimal() +
  coord_flip() +
  theme_linedraw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit("0.35", "cm"),
        legend.position = "bottom",
        legend.margin = margin(-0.25, 0, 0, -0.2, unit = "cm"),
        legend.justification = "left") +
  guides(color = guide_legend(title = "Condition"))+

  ggsave("C:\\Users\\odesa\\Desktop\\Code\\CRC-Final\\Figures\\DRRA008156_DA_pooled.tiff",
         width = 10, height = 10, units = "in", device = 'tiff', dpi = 295)

#print(my_plot)


  
