library(feather)
library(ggplot2)
library(tidyverse)
library(Boruta)

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

## Load the Data
pooled_data <- read_feather('pooled_abundance_data.feather')
pooled_data <- as.data.frame(pooled_data)  

metadata <- read_feather('metadf.feather')
metadata <- as.data.frame(metadata)

# add the 'config' column from metadata to the pooled_data df
pooled_data$config <- metadata[match(pooled_data$SampleIdentifier, metadata$SampleIdentifier),]$config

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

# Add config column back into data frame
CLR_pooled_data <- cbind(identifiers_and_config, CLR_pooled_data)


# Conver to long format DF
long_df <- CLR_pooled_data %>%
  pivot_longer(
    cols = -c(SampleIdentifier, config),
    names_to = "Enzyme",
    values_to = "CLR Relative Abundance"
  )

# Plot the differential abundances and save to a tiff for pub

my_plot <- ggplot(long_df, aes(x = Enzyme, y = `CLR Relative Abundance`, color = config)) +
    geom_point(alpha = 0.8, size = 2, stroke = 0.3, position = position_jitterdodge(jitter.height = 0,
                                                                                    jitter.width = 0.15)) +
    stat_summary(aes(group = config),fun.min = function(x) {quantile(x, 0.25)}, fun.max = function(x) {quantile(x, 0.75)},
                 size = 0.2, width = 0.3, geom = "errorbar", position = position_dodge(width = 0.8), color='black') +
    stat_summary(aes(group = config), fun.min = median, fun.max = median,
                 size = 0.2, width = 0.6, geom = "errorbar", position = position_dodge(width = 0.8), color='black') +
    scale_color_manual(values = c("case" = "#cf5154", "control" = "#4e88aa")) + # Define custom colors
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

## Apply and plot Boruta ##

# Apply Boruta
boruta <- boruta(CLR_pooled_data, metadata)

# Plot Boruta
imp <- attStats(boruta)

# Create a dataframe for plotting
imp_df <- data.frame(
  Attribute = rownames(imp),
  Importance = imp$meanImp,
  Decision = imp$decision
)

print(imp_df)


# Plot using ggplot2
ggplot(imp_df, aes(x = reorder(Attribute, Importance), y = Importance, fill = Decision)) +
  geom_bar(stat = "identity") +
  coord_flip() + # Flip coordinates to make it a horizontal bar plot
  scale_fill_manual(values = c("Rejected" = "red", "Confirmed" = "green", "Tentative" = "yellow")) +
  theme_minimal() +
  labs(title = "Feature Importance from Boruta Analysis", x = "Features", y = "Importance (Z-Score)")
