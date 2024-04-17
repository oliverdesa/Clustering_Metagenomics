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

czm_clr <- function(count_df) {
  
  set.seed(101)
  
  clr <- t(apply(count_df, 1, function(x){log(x) - mean(log(x))}))
  
  clr
  
}

## Load Data ##

# Read in the feather
metadata_7774 <- arrow::read_feather("E:/CRC/PRJEB7774/humann/PRJEB7774_metadata.feather")

# Read in the abundance data
abundance_data_7774 <- arrow::read_feather("E:/CRC/PRJEB7774/humann/new_combined/clean_joined_genefamilies_relab_7774.feather")

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

filt_pooled_abun_df_7774 <- as.data.frame(filt_pooled_abun_df_7774)

# make indices of metadata_7774 into column
metadata_7774$run_accession <- rownames(metadata_7774)
clr_pooled_abun_df_7774$run_accession <- rownames(clr_pooled_abun_df_7774)

filt_pooled_abun_df_7774$run_accession <- rownames(filt_pooled_abun_df_7774)

# Add sample_title to pooled_abun_df_7774 where run_accession matches
clr_pooled_abun_df_7774$sample_title <- metadata_7774$sample_title[match(clr_pooled_abun_df_7774$run_accession, 
                                                                         metadata_7774$run_accession)]

filt_pooled_abun_df_7774$sample_title <- metadata_7774$sample_title[match(filt_pooled_abun_df_7774$run_accession,
                                                                         metadata_7774$run_accession)]

clr_pooled_abun_df_7774 <- subset(clr_pooled_abun_df_7774, select = -c(run_accession))

filt_pooled_abun_df_7774 <- subset(filt_pooled_abun_df_7774, select = -c(run_accession))

# convert to long df format
long_clr_pooled<- filt_pooled_abun_df_7774 %>% 
  pivot_longer(
    cols = c(Amidase:Muramidase), names_to = "Enzyme", 
    values_to = "Relative Abundance")

# Plot the differential abundances and save to a tiff for pub

my_plot <- ggplot(long_clr_pooled, aes(x = Enzyme, y = `Relative Abundance`, color = sample_title)) +
  geom_point(alpha = 0.8, size = 2, stroke = 0.3, position = position_jitterdodge(jitter.height = 0,
                                                                                  jitter.width = 0.15)) +
  #stat_summary(aes(group = sample_title), fun.min = function(x) {quantile(x, 0.25)}, fun.max = function(x) {quantile(x, 0.75)},
  
  stat_summary(aes(group = sample_title),fun.min = function(x) {quantile(x, 0.25)}, fun.max = function(x) {quantile(x, 0.75)},
               size = 0.2, width = 0.3, geom = "errorbar", position = position_dodge(width = 0.8), color='black') +
  stat_summary(aes(group = sample_title), fun.min = median, fun.max = median,
               size = 0.2, width = 0.6, geom = "errorbar", position = position_dodge(width = 0.8), color='black') +
  scale_color_manual(values = c("Carcinoma" = "#cf5154", "Control" = "#4e88aa", "Adenoma" = "#D7B2E5")) + # Define custom colors
  labs(x = "Enzyme", y = "Relative Abundance") + # Label axes
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
  guides(color = guide_legend(title = "Condition"))# Add legend
  # scale_y_continuous(limits = c(-5, 5))

ggsave(plot = my_plot, "C:/users/odesa/Desktop/PRJEB7774_clr_abun.tif",
       width = 10, height = 10, units = "in", device = 'tiff', dpi = 600)

# write new feathers of pooled data
# Want 2 variations of the data for ML training, relab all, CLR - UNMAPPED 

# add rownames to column
clr_pooled_abun_df_7774$sample_id <- rownames(clr_pooled_abun_df_7774)
pooled_abun_df_7774$sample_id <- rownames(pooled_abun_df_7774)

arrow::write_feather(clr_pooled_abun_df_7774, "C:/users/odesa/Desktop/CRCFinal/PRJEB7774/CLR_PRJEB7774_pooled_abun.feather")
arrow::write_feather(pooled_abun_df_7774, "C:/users/odesa/Desktop/CRCFinal/PRJEB7774/relab_PRJEB7774_pooled_abun.feather")


### PCA Analysis ###

relab_pooled <- arrow::read_feather("C:/users/odesa/Desktop/CRCFinal/PRJEB7774/relab_PRJEB7774_pooled_abun.feather")

# Rename the 'run_accession' column in metadata_table to 'sample_id'
names(metadata_7774)[names(metadata_7774) == "run_accession"] <- "sample_id"

# add labels
relab_pooled_labeled <- merge(relab_pooled, metadata_7774[, c("sample_id", "sample_title")], by = "sample_id")

# drop unwanted columns
relab_pooled_labeled <- subset(relab_pooled_labeled, select = -c(sample_id, UNMAPPED, UC118))

# Separate the labels from the data
data <- relab_pooled_labeled[, -ncol(relab_pooled_labeled)]
labels <- relab_pooled_labeled[, ncol(relab_pooled_labeled)]

# Z-score normalization
data_normalized <- as.data.frame(scale(data))

# PCA
pca_result <- prcomp(data_normalized, center = TRUE, scale. = TRUE)

# Create a data frame for plotting
pca_data <- as.data.frame(pca_result$x)
pca_data$label <- labels

# Plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = label)) +
    geom_point() +
    stat_ellipse(type = "t", level = 0.95) +  # Adds ellipses for each group
    theme_minimal() +
    annotate("text", x = Inf, y = Inf, label = annotation_text, hjust = 1, vjust = 1, size = 3.5) +
    theme(plot.margin = margin(1, 1, 1, 1, "cm"))

pca_plot <- pca_plot + 
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white", colour = "white"),
        legend.background = element_rect(fill = "white"),
        plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black"))

ggsave("PCA_plot.png", plot = pca_plot, width = 10, height = 10, units = "in")


# Permanova
distance_matrix <- vegdist(data_normalized, method = "euclidean")
permanova_result <- adonis2(distance_matrix ~ labels)
print(permanova_result)

pairwise_results <- pairwise.adonis(distance_matrix, labels)
print(pairwise_results)

annotation_text <- "Carcinoma vs Control: R2 = 0.0447, p.adj = 0.003"



