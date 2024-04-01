library(feather)
library(ggplot2)
library(tidyverse)
library(Boruta)
library(scales)
library(vegan)
library(pairwiseAdonis)

czm_clr <- function(count_df) {
  
  set.seed(101)
  
  clr <- t(apply(count_df, 1, function(x){log(x) - mean(log(x))}))
  
  clr
  
}

boruta <- function(count_df, meta_df) {
  
  set.seed(101)
  Boruta(count_df,
         as.factor(meta_df[rownames(count_df), "config"]),
         maxRuns = 1000)
  
}

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

## Load the Data
abundance_data <- read_feather('E:/CRC/PRJEB10878/clean_joined_genefamilies_relab_10878.feather')
abundance_data <- abundance_data %>% column_to_rownames(var = "sample_id")

pooled_data <- pooled_abun_df(abundance_data)

# add sample_id back into the data frame
pooled_data$SampleIdentifier <- rownames(pooled_data)


metadata <- read.csv('E:/CRC/PRJEB10878/PRJEB10878_metadata.csv')
metadata <- as.data.frame(metadata)

# add the 'config' column from metadata to the pooled_data df
pooled_data$config <- metadata[match(pooled_data$SampleIdentifier, metadata$Run),]$config

# Save the identifier columns
identifiers_and_config <- pooled_data[, c("SampleIdentifier", "config")]

# List of columns needed
enzymes = c('Amidase', 'DD-carboxypeptidase', 'DD-endopeptidase', 'DL-endopeptidase', 
            'Glucosaminidase', 'LD-carboxypeptidase', 'LD-endopeptidase', 'Muramidase')

#make new DF with only these columns
filtered_pooled_data <- pooled_data[enzymes]

# Apply CLR Normalization and change to DF
CLR_pooled_data <- czm_clr(filtered_pooled_data)
CLR_pooled_data <- as.data.frame(CLR_pooled_data)

filtered_pooled_data <- as.data.frame(filtered_pooled_data)

# Add config column back into data frame
filtered_pooled_data <- cbind(identifiers_and_config, filtered_pooled_data)


# Conver to long format DF
long_df <- filtered_pooled_data %>%
  pivot_longer(
    cols = -c(SampleIdentifier, config),
    names_to = "Enzyme",
    values_to = "Relative Abundance"
  )

# Plot the differential abundances and save to a tiff for pub

my_plot <- ggplot(long_df, aes(x = Enzyme, y = `Relative Abundance`, color = config)) +
    geom_point(alpha = 0.8, size = 2, stroke = 0.3, position = position_jitterdodge(jitter.height = 0,
                                                                                    jitter.width = 0.15)) +
    stat_summary(aes(group = config),fun.min = function(x) {quantile(x, 0.25)}, fun.max = function(x) {quantile(x, 0.75)},
                 size = 0.2, width = 0.3, geom = "errorbar", position = position_dodge(width = 0.8), color='black') +
    stat_summary(aes(group = config), fun.min = median, fun.max = median,
                 size = 0.2, width = 0.6, geom = "errorbar", position = position_dodge(width = 0.8), color='black') +
    scale_color_manual(values = c("case" = "#cf5154", "control" = "#4e88aa")) + # Define custom colors
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
    #scale_y_continuous(limits = c(-5, 5))

  ggsave("C:/users/odesa/Desktop/PRJEB10878_enzyme_abundance.tif", plot = my_plot, width = 10, height = 10, units = "in", 
       dpi = 600, device = "tiff")

# ## Apply and plot Boruta ##
# 
# # Apply Boruta
# boruta <- boruta(CLR_pooled_data, metadata)
# 
# # Plot Boruta
# imp <- attStats(boruta)
# 
# # Create a dataframe for plotting
# imp_df <- data.frame(
#   Attribute = rownames(imp),
#   Importance = imp$meanImp,
#   Decision = imp$decision
# )
# 
# print(imp_df)
# 
# 
# # Plot using ggplot2
# ggplot(imp_df, aes(x = reorder(Attribute, Importance), y = Importance, fill = Decision)) +
#   geom_bar(stat = "identity") +
#   coord_flip() + # Flip coordinates to make it a horizontal bar plot
#   scale_fill_manual(values = c("Rejected" = "red", "Confirmed" = "green", "Tentative" = "yellow")) +
#   theme_minimal() +
#   labs(title = "Feature Importance from Boruta Analysis", x = "Features", y = "Importance (Z-Score)")


## PCA Analysis ##

names(metadata)[names(metadata) == "Run"] <- "SampleIdentifier"

# add labels
relab_pooled_labeled <- merge(pooled_data, metadata[, c("SampleIdentifier", "config")], by = "SampleIdentifier")

# drop unwanted columns
relab_pooled_labeled <- subset(relab_pooled_labeled, select = -c(SampleIdentifier, UNMAPPED, UC118))

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

ggsave("PCA_plot_10878.png", plot = pca_plot, width = 10, height = 10, units = "in")


# Permanova
distance_matrix <- vegdist(data_normalized, method = "euclidean")
permanova_result <- adonis2(distance_matrix ~ labels)
print(permanova_result)

pairwise_results <- pairwise.adonis(distance_matrix, labels)
print(pairwise_results)

annotation_text <- "Case vs Control R2 = 0.00886, p value = 0.322"
