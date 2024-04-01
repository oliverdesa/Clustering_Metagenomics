library(metafor)
library(ggplot2)
library(dplyr)


# Tools like funnel plots for assessing publication bias and statistics like I² 
# for quantifying heterogeneity can provide additional insights into the reliability
# of your meta-analysis results.


# Read the files
data10878 <- read.csv("C:/Users/odesa/OneDrive - University of Toronto/CRC/maaslin/PRJEB10878/new_cluster/no_clr/all_results.tsv", sep="\t")
data7774 <- read.csv("C:/Users/odesa/OneDrive - University of Toronto/CRC/maaslin/PRJEB7774/new_cluster/no_clr/all_results.tsv", sep="\t")
data8156 <- read.csv("C:/Users/odesa/OneDrive - University of Toronto/CRC/maaslin/DRA008156/new_cluster/no_clr/all_results.tsv", sep="\t")

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
subset_data <- combined_data[grepl("DL.e", combined_data$feature), ]

# write.csv(subset_data, "subset_data.csv")

unique_features <- unique(subset_data$feature)

print(unique_features)

writeLines(unique_features, "DLE_ids.txt")

# Placeholder for meta-analysis results
meta_results <- list()

for (feature in unique_features) {
  
  # Make sure the condition is correctly applied
  feature_data <- subset_data[subset_data$feature == feature, ]
  
  # Need to ensrue that correct features being subsetted
  print(paste("Processing feature:", feature, "with", nrow(feature_data), "rows"))
  
  # Perform the meta-analysis
  meta_analysis <- rma(yi = coef, sei = stderr, data = feature_data, method = "REML")
  
  # Store the result
  meta_results[[feature]] <- meta_analysis
}
# Placeholder for extracted data
meta_data <- data.frame(Feature = character(), EffectSize = numeric(), 
                        LowerCI = numeric(), UpperCI = numeric(), pval = numeric())

# Extract data from each meta-analysis result
for (feature in names(meta_results)) {
  res <- meta_results[[feature]]
  meta_data <- rbind(meta_data, data.frame(Feature = feature,
                                           EffectSize = res$beta,
                                           LowerCI = res$ci.lb,
                                           UpperCI = res$ci.ub,
                                           pval = res$pval))
}

dev.off()

print(meta_data)

# Assuming meta_data is prepared as above
meta_data$Feature <- factor(meta_data$Feature, levels = meta_data$Feature[order(meta_data$EffectSize)])

filtered_data <- meta_data %>%
  filter(EffectSize < 0, 
         LowerCI < 0, 
         pval < 0.1) 

# Assuming 'meta_data' contains your full dataset and 'filtered_data' contains the significant features
ggplot(meta_data, aes(x = EffectSize, y = reorder(Feature, EffectSize))) +
  geom_point(aes(color = Feature %in% filtered_data$Feature), size = 3) + # Highlight significant features
  geom_errorbarh(aes(xmin = LowerCI, xmax = UpperCI), height = 0.2) +
  geom_text(data = filtered_data, aes(label = paste0("p=", round(pval, 3))), hjust = 1.5, color = "red") + # Annotate with p-values
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) + # Color significant features differently
  labs(x = "Effect Size", y = "Feature", title = "Forest Plot of Significant Features") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        legend.position = "none") 

  ggsave("C:/users/odesa/Desktop/meta_analysis.tif", width = 20, height = 20, units = "cm", dpi = 600,
         device = "tiff")



print(filtered_data)
