library(feather)
library(ggplot2)
library(tidyverse)

pooled_data <- read_feather('pooled_abundance_data.feather')
pooled_data <- as.data.frame(pooled_data)  

metadata <- read_feather('metadf.feather')
metadata <- as.data.frame(metadata)

# add the 'config' column from metadata to the pooled_data df, make sure they align
# properly on the SampleIdentifier column present in both dfs
pooled_data$config <- metadata[match(pooled_data$SampleIdentifier, metadata$SampleIdentifier),]$config

long_df <- pooled_data %>%
  pivot_longer(
    cols = starts_with("Amidase"):starts_with("UC118"),
    names_to = "Enzyme",
    values_to = "Value"
  )

# Define new labels with the known sample sizes
new_labels <- c(case = "Case (n=74)", control = "Control (n=53)")

# Function to create a plot
plot_function <- function(data, enzyme_name) {
  p <- ggplot(data, aes(x = config, y = Value)) +
    geom_boxplot(aes(fill = config)) +
    scale_fill_manual(values = c("control" = "grey", "case" = "blue")) +
    geom_jitter(width = 0.2, size = 1, color = 'black') +
    labs(
      title = paste(enzyme_name), 
      x = "", 
      y = "Reads per Kilobase (RPK)"
    ) +
    scale_x_discrete(labels = new_labels) +
    theme_minimal(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  ggsave(paste0("Enzyme_", enzyme_name, ".png"), dpi = 600, plot = p, width = 10, height = 6)
}

# Loop through each enzyme and plot
enzyme_names <- unique(long_df$Enzyme)

for (enzyme in enzyme_names) {
  enzyme_data <- filter(long_df, Enzyme == enzyme)
  plot_function(enzyme_data, enzyme)
}

