library(Maaslin2)
library(dplyr)
library(ggbiplot)
library(feather)
library(arrow)
library(tidyverse)

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

czm_clr <- function(count_df) {
  
  set.seed(101)
  
  clr <- t(apply(count_df, 1, function(x){log(x) - mean(log(x))}))
  
  clr
  
}

## Load Data ##

# Read in the feather
metadata_7774 <- arrow::read_feather("C:/Users/odesa/Desktop/CRCFinal/PRJEB7774/PRJEB7774_metadata.feather")

# Read in the abundance data
abundance_data_7774 <- arrow::read_feather("C:/Users/odesa/Desktop/CRCFinal/PRJEB7774/clean_joined_genefamilies_relab_7774.feather")

# set sample_id as rownames
abundance_data_7774 <- abundance_data_7774 %>% column_to_rownames(var = "sample_id")
metadata_7774 <- metadata_7774 %>% column_to_rownames(var = "run_accession")

# pool the data
pooled_abun_df_7774 <- pooled_abun_df(abundance_data_7774)

# drop the UNMAPPED column
pooled_abun_df_7774 <- subset(pooled_abun_df_7774, select = -c(UNMAPPED))

## Run Maaslin2 ##

results <- Maaslin2(
  input_data = pooled_abun_df_7774,
  input_metadata = metadata_7774,
  output = "C:/users/odesa/Desktop/CRCFinal/CRC-Final/data/PRJEB7774_maaslin",
  fixed_effects = c("sample_title"),
  random_effects = NULL,
  normalization = "CLR", 
  transform = "none",
  reference=c("sample_title,Control")
)
## Plotting ##

# List of columns needed
enzymes = c('Amidase', 'DD-carboxypeptidase', 'DD-endopeptidase', 'DL-endopeptidase', 
            'Glucosaminidase', 'LD-carboxypeptidase', 'LD-endopeptidase', 'Muramidase')

filt_pooled_abun_df_7774 <- pooled_abun_df_7774[, enzymes]

# Apply CLR Normalization and change to DF
clr_pooled_abun_df_7774 <- czm_clr(filt_pooled_abun_df_7774)
clr_pooled_abun_df_7774 <- as.data.frame(clr_pooled_abun_df_7774)

# make indices of metadata_7774 into column
metadata_7774$run_accession <- rownames(metadata_7774)
clr_pooled_abun_df_7774$run_accession <- rownames(clr_pooled_abun_df_7774)

# Add sample_title to pooled_abun_df_7774 where run_accession matches
clr_pooled_abun_df_7774$sample_title <- metadata_7774$sample_title[match(clr_pooled_abun_df_7774$run_accession, 
                                                                         metadata_7774$run_accession)]

clr_pooled_abun_df_7774 <- subset(clr_pooled_abun_df_7774, select = -c(run_accession))

# convert to long df format
long_clr_pooled<- clr_pooled_abun_df_7774 %>% 
  pivot_longer(
    cols = c(Amidase:Muramidase), names_to = "Enzyme", 
    values_to = "CLR Relative Abundance")

# Plot the differential abundances and save to a tiff for pub

my_plot <- ggplot(long_clr_pooled, aes(x = Enzyme, y = `CLR Relative Abundance`, color = sample_title)) +
  geom_point(alpha = 0.8, size = 2, stroke = 0.3, position = position_jitterdodge(jitter.height = 0,
                                                                                  jitter.width = 0.15)) +
  stat_summary(aes(group = sample_title),fun.min = function(x) {quantile(x, 0.25)}, fun.max = function(x) {quantile(x, 0.75)},
               size = 0.2, width = 0.3, geom = "errorbar", position = position_dodge(width = 0.8), color='black') +
  stat_summary(aes(group = sample_title), fun.min = median, fun.max = median,
               size = 0.2, width = 0.6, geom = "errorbar", position = position_dodge(width = 0.8), color='black') +
  scale_color_manual(values = c("Carcinoma" = "#cf5154", "Control" = "#4e88aa", "Adenoma" = "#D7B2E5")) + # Define custom colors
  labs(x = "Enzyme", y = "CLR Relative Abundance") + # Label axes
  theme_minimal() + # Minimal theme
  coord_flip() + # Flip axes
  
  theme_linedraw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 14), # Base text size
        axis.title = element_text(size = 16), # Title for axes
        axis.text = element_text(size = 12), # Text for axis ticks
        legend.title = element_text(size = 14), # Text for legend title
        legend.text = element_text(size = 12),
        legend.key.size = unit("0.35", "cm"),
        legend.position = "bottom",
        legend.margin = margin(-0.25, 0, 0, -0.2, unit = "cm"),
        legend.justification = "left") +
  guides(color = guide_legend(title = "Condition")) +# Add legend
  scale_y_continuous(limits = c(-5, 5))

ggsave("C:/users/odesa/Desktop/CRCFinal/Figures/PRJEB7774_clr_abun.tif",
       width = 10, height = 10, units = "in", device = 'tiff', dpi = 600)

# write new feathers of pooled data

arrow::write_feather(clr_pooled_abun_df_7774, "C:/users/odesa/Desktop/CRCFinal/PRJEB7774/CLR_PRJEB7774_pooled_abun.feather")

