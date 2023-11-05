library(feather)
library(ggplot2)
library(tidyverse)

czm_clr <- function(count_df) {
  
  set.seed(101)
  
  clr <- t(apply(count_df, 1, function(x){log(x) - mean(log(x))}))
  
  clr
  
}

pooled_data <- read_feather('pooled_abundance_data.feather')
pooled_data <- as.data.frame(pooled_data)  

metadata <- read_feather('metadf.feather')
metadata <- as.data.frame(metadata)

# add the 'config' column from metadata to the pooled_data df, make sure they align
# properly on the SampleIdentifier column present in both dfs
pooled_data$config <- metadata[match(pooled_data$SampleIdentifier, metadata$SampleIdentifier),]$config

identifiers_and_config <- pooled_data[, c("SampleIdentifier", "config")]

# List of columns needed
enzymes = c('Amidase', 'DD-carboxypeptidase', 'DD-endopeptidase', 'DL-endopeptidase', 
            'Glucosaminidase', 'LD-carboxypeptidase', 'LD-endopeptidase', 'Muramidase')

#make new DF with only these columns
filtered_pooled_data <- pooled_data[enzymes]

CLR_pooled_data <- czm_clr(filtered_pooled_data)

CLR_pooled_data <- as.data.frame(CLR_pooled_data)

# Add config column back into data frame
CLR_pooled_data <- cbind(identifiers_and_config, CLR_pooled_data)


long_df <- CLR_pooled_data %>%
  pivot_longer(
    cols = -c(SampleIdentifier, config),
    names_to = "Enzyme",
    values_to = "CLR Relative Abundance"
  )

# Convert the 'CLR Relative Abundance' to numeric, coercing non-numeric to NAs
long_df$`CLR Relative Abundance` <- as.numeric(as.character(long_df$`CLR Relative Abundance`))

# Check again after conversion
str(long_df$`CLR Relative Abundance`)

summary(long_df$`CLR Relative Abundance`)

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

ggsave("PRJEB10878_diff_abundance.tiff", plot = my_plot, width = 11, height = 8.5, units = "in", dpi = 300, type = "cairo")

