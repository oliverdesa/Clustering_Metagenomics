library(Maaslin2)
library(dplyr)
library(ggbiplot)
library(feather)

## Load + Clean Data ##

abundance_data <- read.table('cleaned_data.tsv', header=TRUE, sep='\t', 
                             stringsAsFactors = FALSE, check.names = FALSE)

# Replace the _Abundance-RPKs part of each column name with an empty string
colnames(abundance_data) <- gsub("_Abundance-RPKs", "", colnames(abundance_data))

original_column_names <- colnames(abundance_data)

abundance_data <- t(abundance_data)

# Sub 0s for low num
abundance_data[abundance_data==0] <- 1e-6

# Make first row colum nnames
colnames(abundance_data) <- abundance_data[1,]
abundance_data <- abundance_data[-1, ]

abundance_data = data.frame(abundance_data)

# convert values from string to foat
abundance_data <- mutate_all(abundance_data, function(x) as.numeric(as.character(x)))

metadf <- read.csv('MetaData.csv')

# Set the first column as row names
rownames(metadf) <- metadf[, 1]
# Drop the first column
metadf <- metadf[, -1]

## Run Maaslin2, log transformation ##

results <- Maaslin2(
  input_data = abundance_data,
  input_metadata = metadf,
  output = "C:/users/odesa/Desktop/CRC",
  fixed_effects = c("config"),
  random_effects = NULL,
  normalization = "NONE", 
  transform = "LOG"
)

# Add a column containing the sample identifiers
metadf$SampleIdentifier <- rownames(metadf)
abundance_data$SampleIdentifier <- rownames(abundance_data)

# Write the DFs to feather to be used in other scripts
write_feather(metadf, "metadf.feather")
write_feather(abundance_data, 'abundance.feather')


## Run PCA + Permanova ##

new_abundance_data <- read_feather('abundance.feather')
new_metadf <- read_feather('metadf.feather')

# From new_metadf remove the row with SampleIdentifier ERR1018218
new_metadf <- new_metadf[-which(new_metadf$SampleIdentifier == 'ERR1018218'), ]

new_abundance_data <- as.data.frame(new_abundance_data)
rownames(new_abundance_data) <- new_abundance_data$SampleIdentifier
new_abundance_data$SampleIdentifier <- NULL

# Adjust the factor levels to ensure "control" is plotted first
new_metadf$config <- factor(new_metadf$config, levels = c("control", "case"))


# Standardize the data
df_standardized <- scale(new_abundance_data)

# Compute PCA
pca_result <- prcomp(df_standardized, center = TRUE, scale. = TRUE)

# Plot the first two principal components with adjusted alpha for transparency
ggbiplot(pca_result, obs.scale = 1, var.scale = 1, 
         groups = new_metadf$config, ellipse = TRUE, 
         circle = TRUE, alpha = 0.5) +   # Add alpha for transparency
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')

# To see a summary of the variance explained by each principal component
summary(pca_result)



