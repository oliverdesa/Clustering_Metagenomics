library(metafor)
library(ggplot2)
library(dplyr)

####################################################################
# This R script is to carry out the meta-analysis of MaAsLin2 results
# across the CRC trials and generate outputted plots
####################################################################


# Tools like funnel plots for assessing publication bias and statistics like I² 
# for quantifying heterogeneity can provide additional insights into the reliability
# of your meta-analysis results.


# Read the files
# data10878 <- read.csv("C:/Users/odesa/OneDrive - University of Toronto/CRC/maaslin/PRJEB10878/new_cluster/no_clr/all_results.tsv", sep="\t")
# data7774 <- read.csv("C:/Users/odesa/OneDrive - University of Toronto/CRC/maaslin/PRJEB7774/new_cluster/no_clr/all_results.tsv", sep="\t")
# data8156 <- read.csv("C:/Users/odesa/OneDrive - University of Toronto/CRC/maaslin/DRA008156/new_cluster/no_clr/all_results.tsv", sep="\t")

data10878 <- read.csv("/Volumes/PGH-Backup/CRC/PRJEB10878/maaslin/carboxys_and_dles/all_results.tsv", sep="\t")
data7774 <- read.csv("/Volumes/PGH-Backup/CRC/PRJEB7774/maaslin/carboxys_and_dles/all_results.tsv", sep="\t")
data8156 <- read.csv("/Volumes/PGH-Backup/CRC/DRA008156/maaslin2/carboxys_and_dles/all_results.tsv", sep="\t")

# Add a dataset identifier
data10878$Dataset <- 'PRJEB10878'
data7774$Dataset <- 'PRJEB7774'
data8156$Dataset <- 'DRA008156'

data7774 <- data7774[!grepl("Adenoma", data7774$value), ]

# Select relevant columns
data10878 <- subset(data10878, select=c(Dataset, feature, coef, stderr, pval, qval))
data7774 <- subset(data7774, select=c(Dataset, feature, coef, stderr, pval, qval))
data8156 <- subset(data8156, select=c(Dataset, feature, coef, stderr, pval, qval))



combined_data <- rbind(data10878, data7774, data8156)

# Subsetting to include only features that contain "DL.e"

patterns <- c("^DD.c")

# Loop through enzymes, subset df
for (pattern in patterns) {
  
  # subset data for current enzyme
  subset_data <- combined_data[grepl(pattern, combined_data$feature), ]
  
  # Print the number of rows
  print(paste("Number of rows for", pattern, ":", nrow(subset_data))) # ok
  
  # save and print unqiue features
  unique_features <- unique(subset_data$feature)
  print(unique_features) # ok
  
  # Write the unique features to a file
  # remove special characters
  safe_pattern <- gsub("[^a-zA-Z0-9]", "_", pattern)
  
  #generate file name
  file_name <- paste0("/Volumes/PGH-Backup/CRC/ids_", safe_pattern, "_ids.txt")
  
  # output features file
  writeLines(unique_features, file_name)
  
  # Placeholder for meta-analysis results
  meta_results <- list()
  
  for (feature in unique_features) {
    
    # Make sure the condition is correctly applied
    feature_data <- subset_data[subset_data$feature == feature, ]
    
    # Need to ensure that correct features being subsetted
    print(paste("Processing feature:", feature, "with", nrow(feature_data), "rows"))
    
    # Perform the meta-analysis
    meta_analysis <- rma(yi = coef, sei = stderr, data = feature_data, method = "REML")
    
    # Store the result
    meta_results[[feature]] <- meta_analysis
  }
  
  print(meta_results)
  
  # Placeholder for extracted data
  meta_data <- data.frame(Feature = character(), EffectSize = numeric(), 
                          LowerCI = numeric(), UpperCI = numeric(), pval = numeric())
  
  for (feature in names(meta_results)) {
    res <- meta_results[[feature]]
    significance <- ifelse(res$pval < 0.05, "Red", 
                           ifelse(res$pval < 0.1, "Orange", "Not Significant"))
    meta_data <- rbind(meta_data, data.frame(Feature = feature,
                                             EffectSize = res$beta,
                                             LowerCI = res$ci.lb,
                                             UpperCI = res$ci.ub,
                                             pval = res$pval,
                                             Significance = significance))
  }
  
  # Assuming meta_data is prepared as above
  meta_data$Feature <- factor(meta_data$Feature, levels = meta_data$Feature[order(meta_data$EffectSize)])
  
  filtered_data <- meta_data %>%
    filter(EffectSize < 0, 
           LowerCI < 0, 
           pval < 0.1) 
  
  plot_title <- paste0("Forest Plot of Significant Features for", pattern)
  
  # Plot the forest plot
  ggplot(meta_data, aes(y = EffectSize, x = reorder(Feature, EffectSize))) +
    geom_point(aes(color = Significance), size = 3) + # Use Significance for color
    geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.3) + # Vertical error bars with wider spacing
    # geom_text(data = filtered_data, aes(label = paste0("p=", round(pval, 3)), color = Significance), 
              # vjust = -6.5, angle = 45, hjust = 1, size = 3.5) + # Adjust position, angle, and size of p-value annotations
    scale_color_manual(values = c("Red" = "red", "Orange" = "orange", "Not Significant" = "blue")) + # Define colors
    labs(y = "MaAsLin2 Coefficient", x = "Feature", title = plot_title) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1), # Increase font size for better readability
          axis.text.y = element_text(size = 10, margin = margin(r = 10)), # Add more space between y-axis labels
          plot.title = element_text(size = 14, face = "bold"), # Make title more prominent
          legend.position = "none") 
  
  # Adjust plot width and height for more spacing
  plot_width <- 10 + 0.5 * nrow(meta_data) 
  plot_height <- 25
  
  # Save the plot with increased dimensions
  figure_name <- paste0("/Volumes/PGH-Backup/CRC/meta_analysis", safe_pattern, ".tif")
  
  ggsave(figure_name, height = plot_height, width = plot_width, units = "cm", dpi = 600, device = "tiff")
  
}
