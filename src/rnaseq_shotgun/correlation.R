library(tidyverse)
library(readr)
library(Maaslin2)
library(ggplot2)
library(gridExtra)
library(car)
library(GiANT)

# functions

# Sample function to filter mgx dataframe based on RNAseq samples
find_closest_mgx <- function(rnaseq_df, mgx_df) {
  # Split the 'unique_values' into patient ID and day for both dataframes
  rnaseq_df$Patient_ID <- sapply(strsplit(rnaseq_df$unique_id, "_"), `[`, 1)
  rnaseq_df$Day <- as.numeric(sapply(strsplit(rnaseq_df$unique_id, "_"), `[`, 2))
  
  mgx_df$Patient_ID <- sapply(strsplit(mgx_df$unique_id, "_"), `[`, 1)
  mgx_df$Day <- as.numeric(sapply(strsplit(mgx_df$unique_id, "_"), `[`, 2))
  
  # Store the original column names to ensure consistent rbind
  original_columns <- names(mgx_df)
  
  # Initialize an empty dataframe to store the filtered mgx data
  filtered_mgx <- mgx_df[0, ]  # Create an empty dataframe with the same structure as mgx_df
  
  # Loop through each sample in the RNAseq dataframe
  for (i in 1:nrow(rnaseq_df)) {
    # Get the current patient ID and day from RNAseq
    patient_id <- rnaseq_df$Patient_ID[i]
    target_day <- rnaseq_df$Day[i]
    
    # Subset mgx data for the same patient ID
    mgx_subset <- mgx_df[mgx_df$Patient_ID == patient_id, ]
    
    if (nrow(mgx_subset) > 0) {
      # Find the exact match for patient ID and day
      exact_match <- mgx_subset[mgx_subset$Day == target_day, ]
      
      if (nrow(exact_match) == 1) {
        # If exact match exists, select that row and keep only the original columns
        filtered_mgx <- rbind(filtered_mgx, exact_match[, original_columns])
      } else {
        # If no exact match, find the closest day
        mgx_subset$Day_Diff <- abs(mgx_subset$Day - target_day)
        closest_match <- mgx_subset[which.min(mgx_subset$Day_Diff), ]
        # Select the closest match and keep only the original columns
        filtered_mgx <- rbind(filtered_mgx, closest_match[, original_columns])
      }
    }
  }
  
  # Return the filtered mgx dataframe with only the original columns
  return(filtered_mgx)
}

# Function to clean junk from rna data
clean_rna_cols <- function(df){
  df <- df %>% select(-"biopsy_location", -"week_num_metagenomics",
                      -"External ID_x", -"sample_id", -"External ID_y")
  #df <- df %>% column_to_rownames(var = "unique_id")
  return(df)
}

# Function to parse a GMT file
parse_gmt <- function(file_path) {
  # Initialize an empty list to store gene sets
  gene_sets <- list()
  
  # Read the GMT file line by line
  lines <- readLines(file_path)
  
  # Iterate over each line to parse the gene set information
  for (line in lines) {
    # Split the line by tabs
    parts <- strsplit(line, "\t")[[1]]
    
    # The first item is the gene set name, the second is a description, the rest are genes
    gene_set_name <- parts[1]
    genes <- parts[3:length(parts)]  # Genes start from the third element
    
    # Add to the list
    gene_sets[[gene_set_name]] <- genes
  }
  
  return(gene_sets)
}

# Load the data
mgx <- read_tsv("/Volumes/PGH-Backup/ibd_data/humann_second_run/ibd_genefamilies_relab_clustered_for_rnaseq.tsv")
rnaseq <- read_tsv("/Volumes/PGH-Backup/ibd_data/rnaseq/tmm_normalized_counts_for_mgx.tsv")

# Clean the data
# Drop junk column
mgx <- mgx %>% select(-"Unnamed: 0")

# Seperate by biopsy location and drop duplicates
rna_colon <- rnaseq %>% filter(grepl("Sigmoid Colon", biopsy_location))
rna_ileum <- rnaseq %>% filter(grepl("Ileum", biopsy_location))
rna_rectum <- rnaseq %>% filter(grepl("Rectum", biopsy_location))

rna_colon <- rna_colon[!duplicated(rna_colon$unique_id), ]
rna_ileum <- rna_ileum[!duplicated(rna_ileum$unique_id), ]
rna_rectum <- rna_rectum[!duplicated(rna_rectum$unique_id), ]

# Remove junk from RNA datasets
rna_colon_clean <- clean_rna_cols(rna_colon)
rna_ileum_clean <- clean_rna_cols(rna_ileum)
rna_rectum_clean <- clean_rna_cols(rna_rectum)

# Filter mgx data to only include samples that have rnaseq data
filtered_mgx_colon <- find_closest_mgx(rna_colon_clean, mgx)
filtered_mgx_ileum <- find_closest_mgx(rna_ileum_clean, mgx)
filtered_mgx_rectum <- find_closest_mgx(rna_rectum_clean, mgx)

# Remove duplicates from rnaseq data and mgx
rna_colon_clean <- rna_colon_clean[!duplicated(rna_colon_clean$`Participant ID`), ]
rna_ileum_clean <- rna_ileum_clean[!duplicated(rna_ileum_clean$`Participant ID`), ]
rna_rectum_clean <- rna_rectum_clean[!duplicated(rna_rectum_clean$`Participant ID`), ]

filtered_mgx_colon <- filtered_mgx_colon[!duplicated(filtered_mgx_colon$`Participant ID`), ]
filtered_mgx_ileum <- filtered_mgx_ileum[!duplicated(filtered_mgx_ileum$`Participant ID`), ]
filtered_mgx_rectum <- filtered_mgx_rectum[!duplicated(filtered_mgx_rectum$`Participant ID`), ]

# Add covariates to rnaseq data for random effects controlling
metadata <- read_tsv("/Volumes/PGH-Backup/ibd_data/metadata/hmp2_rnaseq_metadata.tsv")

#remove dupes
metadata <- metadata[!duplicated(metadata$Participant.ID), ]

# Rename columns to match
colnames(rna_colon_clean)[colnames(rna_colon_clean) == "Participant ID"] <- "Participant.ID"
colnames(rna_ileum_clean)[colnames(rna_ileum_clean) == "Participant ID"] <- "Participant.ID"
colnames(rna_rectum_clean)[colnames(rna_rectum_clean) == "Participant ID"] <- "Participant.ID"

# Select only diagnosis and patient ID
metadata <- metadata %>% select("Participant.ID", "diagnosis")

# merge covariates to RNA data
merged_rna_colon <- merge(rna_colon_clean, metadata, by = "Participant.ID")
merged_rna_ileum <- merge(rna_ileum_clean, metadata, by = "Participant.ID")
merged_rna_rectum <- merge(rna_rectum_clean, metadata, by = "Participant.ID")

dim(merged_rna_rectum)

# Make Participant ID column rownames
filtered_mgx_colon <- filtered_mgx_colon %>% column_to_rownames(var = "Participant ID")
filtered_mgx_ileum <- filtered_mgx_ileum %>% column_to_rownames(var = "Participant ID")
filtered_mgx_rectum <- filtered_mgx_rectum %>% column_to_rownames(var = "Participant ID")

merged_rna_colon <- merged_rna_colon %>% column_to_rownames(var = "Participant.ID")
merged_rna_ileum <- merged_rna_ileum %>% column_to_rownames(var = "Participant.ID")
merged_rna_rectum <- merged_rna_rectum %>% column_to_rownames(var = "Participant.ID")

# Drop junk from mgx data
filtered_mgx_colon <- filtered_mgx_colon %>% select(-"sample_id", -"Patient_ID",
                                                    -"week_num", -"External ID",
                                                    -"unique_id", -"Day")

filtered_mgx_ileum <- filtered_mgx_ileum %>% select(-"sample_id", -"Patient_ID",
                                                    -"week_num", -"External ID",
                                                    -"unique_id", -"Day")

filtered_mgx_rectum <- filtered_mgx_rectum %>% select(-"sample_id", -"Patient_ID",
                                                      -"week_num", -"External ID",
                                                      -"unique_id", -"Day")



print(head(names(filtered_mgx_rectum), 100))

# Gene set to search for in Maaslin2
# inflammatory_markers <- c("IL1B", "IL6", "IL8", "TNF", "CCL2", "CXCL10", "CCL5",
#                           "SAA1", "CRP", "NOD2", "TLR2", "TLR4", "NLRP3", "CASP1",
#                           "IL1B", "CXCL1", "LFNA4", "IFNB", "IFNG", "IL17A")

inflammatory_geneset <- parse_gmt("/Volumes/PGH-Backup/ibd_data/rnaseq/HALLMARK_INFLAMMATORY_RESPONSE.v2024.1.Hs.gmt")
inflammatory_geneset <- unname(inflammatory_geneset[[1]])
print(inflammatory_geneset)


# Define the Maaslin2 parameters

results <- Maaslin2(
  input_data = filtered_mgx_colon,
  input_metadata = merged_rna_colon,
  output = "/Volumes/PGH-Backup/ibd_data/shotgun_rnaseq_maaslin/colon",
  fixed_effects = inflammatory_markers,
  random_effects = c("diagnosis"),
  normalization = "none",
  transform = "none",
  reference = "none"
)

results <- Maaslin2(
  input_data = filtered_mgx_ileum,
  input_metadata = merged_rna_ileum,
  output = "/Volumes/PGH-Backup/ibd_data/shotgun_rnaseq_maaslin/ileum",
  fixed_effects = inflammatory_markers,
  random_effects = c("diagnosis"),
  normalization = "none",
  transform = "none",
  reference = "none"
)

results <- Maaslin2(
  input_data = filtered_mgx_rectum,
  input_metadata = merged_rna_rectum,
  output = "/Volumes/PGH-Backup/ibd_data/shotgun_rnaseq_maaslin/rectum",
  fixed_effects = inflammatory_markers,
  random_effects = c("diagnosis"),
  normalization = "none",
  transform = "none",
  reference = "none"
)

any(is.na(merged_rna_rectum))

merged_rna_rectum$diagnosis

# Most significant results from rectum which is good. 

#########################################################

# Look at Nod2 subset of genes
print(merged_rna_rectum$Age)

merged_rna_rectum_nod2 <- merged_rna_rectum[, colnames(merged_rna_rectum) %in% nod2_genes]

metadata_columns <- merged_rna_rectum[, c("diagnosis", "Age")]

merged_rna_rectum_nod2 <- cbind(merged_rna_rectum_nod2, metadata_columns)

print(merged_rna_rectum_nod2)

if (!is.character(nod2_genes)) {
  nod2_genes <- as.character(nod2_genes)
}

results_nod2 <- Maaslin2(
  input_data = filtered_mgx_rectum,
  input_metadata = merged_rna_rectum_nod2,
  output = "/Volumes/PGH-Backup/ibd_data/shotgun_rnaseq_maaslin/rectum_nod2",
  fixed_effects = nod2_genes,
  random_effects = c("diagnosis", "Age"),
  normalization = "none",
  transform = "none",
  reference = "none"
)

print(nod2_genes)

#########################################################

# PCA correlations

print(tail(names(filtered_mgx_rectum), 10))

filtered_mgx_rectum <- filtered_mgx_rectum %>%
  select(-"sample_id", -"Patient_ID", -"week_num", -"External ID", -"unique_id", -"Day")
print(tail(names(filtered_mgx_rectum), 10))

# Remove columns with all zeros
non_zero_columns <- apply(filtered_mgx_rectum, 2, function(col) any(col != 0))
filtered_mgx_rectum <- filtered_mgx_rectum[, non_zero_columns]

# Remove rows with all zeros
non_zero_rows <- apply(filtered_mgx_rectum, 1, function(row) any(row != 0))
filtered_mgx_rectum <- filtered_mgx_rectum[non_zero_rows, ]

# choose numeric columns
filtered_mgx_rectum <- filtered_mgx_rectum %>% select_if(is.numeric)

mgx_scaled <- scale(filtered_mgx_rectum)


# Step 2: Identify and remove zero-variance columns after scaling
constant_columns <- apply(mgx_scaled, 2, var) == 0

print(constant_columns)

# Drop the constant columns
mgx_filtered <- mgx_scaled[, !constant_columns]

# Check the dimensions of the filtered dataframe
dim(mgx_filtered)

pca_result <- prcomp(mgx_scaled, center = TRUE)

summary(pca_result)

df_pca <- as.data.frame(pca_result$x)

df_pca <- df_pca %>%
  rownames_to_column(var = "Participant ID")

inflammatory_geneset <- trimws(inflammatory_geneset)

inflammatory_geneset <- sapply(inflammatory_geneset, function(x) iconv(x, to = "ASCII//TRANSLIT"))

inflammatory_geneset <- unname(inflammatory_geneset)

print(inflammatory_geneset)

print(typeof(inflammatory_geneset))

inflammatory_subset <- inflammatory_geneset[1:100]

# inflam_test <- inflammatory_geneset[1:20]
# 
# print(inflam_test)
# 
# inflam_test <- trimws(inflam_test)

inflammatory_geneset <- c(
  "ABCA1", "ABI1", "ACVR1B", "ACVR2A", "ADGRE1", "ADM", "ADORA2B", "ADRM1",
  "AHR", "APLNR", "AQP9", "ATP2A2", "ATP2B1", "ATP2C1", "AXL", "BDKRB1",
  "BEST1", "BST2", "BTG2", "C3AR1", "C5AR1", "CALCRL", "CCL17", "CCL2",
  "CCL20", "CCL22", "CCL24", "CCL5", "CCL7", "CCR7", "CCRL2", "CD14",
  "CD40", "CD48", "CD55", "CD69", "CD70", "CD82", "CDKN1A", "CHST2",
  "CLEC5A", "CMKLR1", "CSF1", "CSF3", "CSF3R", "CX3CL1", "CXCL10", "CXCL11",
  "CXCL6", "CXCL8", "CXCL9", "CXCR6", "CYBB", "DCBLD2", "EBI3", "EDN1",
  "EIF2AK2", "EMP3", "EREG", "F3", "FFAR2", "FPR1", "FZD5", "GABBR1",
  "GCH1", "GNA15", "GNAI3", "GP1BA", "GPC3", "GPR132", "GPR183", "HAS2",
  "HBEGF", "HIF1A", "HPN", "HRH1", "ICAM1", "ICAM4", "ICOSLG", "IFITM1",
  "IFNAR1", "IFNGR2", "IL10", "IL10RA", "IL12B", "IL15", "IL15RA", "IL18",
  "IL18R1", "IL18RAP", "IL1A", "IL1B", "IL1R1", "IL2RB", "IL4R", "IL6",
  "IL7R", "INHBA", "IRAK2", "IRF1", "IRF7", "ITGA5", "ITGB3", "ITGB8",
  "KCNA3", "KCNJ2", "KCNMB2", "KIF1B", "KLF6", "LAMP3", "LCK", "LCP2",
  "LDLR", "LIF", "LPAR1", "LTA", "LY6E", "LYN", "MARCO", "MEFV",
  "MEP1A", "MET", "MMP14", "MSR1", "MXD1", "MYC", "NAMPT", "NDP",
  "NFKB1", "NFKBIA", "NLRP3", "NMI", "NMUR1", "NOD2", "NPFFR2", "OLR1",
  "OPRK1", "OSM", "OSMR", "P2RX4", "P2RX7", "P2RY2", "PCDH7", "PDE4B",
  "PDPN", "PIK3R5", "PLAUR", "PROK2", "PSEN1", "PTAFR", "PTGER2", "PTGER4",
  "PTGIR", "PTPRE", "PVR", "RAF1", "RASGRP1", "RELA", "RGS1", "RGS16",
  "RHOG", "RIPK2", "RNF144B", "ROS1", "RTP4", "SCARF1", "SCN1B", "SELE",
  "SELENOS", "SELL", "SEMA4D", "SERPINE1", "SGMS2", "SLAMF1", "SLC11A2", "SLC1A2",
  "SLC28A2", "SLC31A1", "SLC31A2", "SLC4A4", "SLC7A1", "SLC7A2", "SPHK1", "SRI",
  "STAB1", "TACR1", "TACR3", "TAPBP", "TIMP1", "TLR1", "TLR2", "TLR3",
  "TNFAIP6", "TNFRSF1B", "TNFRSF9", "TNFSF10", "TNFSF15", "TNFSF9", "TPBG", "VIP"
)

inflam_test <- inflammatory_geneset[1:75]


# Visualizing the first two principal components
ggplot(df_pca, aes(x = PC1, y = PC2, label =`Participant ID`)) +
  geom_point() +
  geom_text(vjust = 1, hjust = 1) +
  theme_minimal() +
  labs(title = "PCA of Metagenomic Data", x = "PC1", y = "PC2")

  ggsave("/Volumes/PGH-Backup/ibd_data/shotgun_rnaseq_maaslin/figures/rectum_mgx_pca.png",
         width = 10, height = 8, units = "in", dpi = 600)
  
# set column to rownames
df_pca <- df_pca %>% column_to_rownames(var = "Participant ID")

# Correlation between PCA components and RNAseq data
results_nod2 <- Maaslin2(
  input_data = df_pca,
  input_metadata = merged_rna_rectum,
  output = "/Volumes/PGH-Backup/ibd_data/shotgun_rnaseq_maaslin/rectum_pca_inflam",
  fixed_effects = inflam_test,
  random_effects = c("diagnosis"),
  normalization = "none",
  transform = "none",
  reference = "none"
)

#########################################################

# Generate a loading plot

# Step 1: Extract the loadings (contributions of each gene to the PC)
loadings <- pca_result$rotation[, 1]  # Extract loadings for PC1 (change "1" to any other PC if needed)

# Step 2: Create a data frame for plotting
loadings_df <- data.frame(Gene = rownames(pca_result$rotation), 
                          Loading = loadings)

# Step 3: Order by absolute loading value (most contributing genes)
loadings_df <- loadings_df[order(abs(loadings_df$Loading), decreasing = TRUE), ]



# Step 4: Plot the top contributing genes (you can adjust the number to focus on top 10 or 20 genes)
ggplot(loadings_df[1:20, ], aes(x = reorder(Gene, Loading), y = Loading)) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +  # Flip coordinates for better visibility
  theme_minimal() +
  labs(title = "Top Genes Contributing to PC1",
       x = "Genes",
       y = "Loading Value")

# Adjusted loadings plot

# Function to generate a bar plot for top positive and negative contributors of a given PC
plot_top_pos_neg_contributors <- function(pca_result, pc_number, top_n = 10) {
  # Step 1: Extract the loadings (contributions of each gene to the PC)
  loadings <- pca_result$rotation[, pc_number]  # Extract loadings for the given PC
  
  # Step 2: Create a data frame for plotting
  loadings_df <- data.frame(Gene = rownames(pca_result$rotation), 
                            Loading = loadings)
  
  # Step 3: Separate into positive and negative loadings
  top_pos_loadings <- loadings_df %>%
    filter(Loading > 0) %>%
    arrange(desc(Loading)) %>%
    head(top_n)
  
  top_neg_loadings <- loadings_df %>%
    filter(Loading < 0) %>%
    arrange(Loading) %>%
    head(top_n)
  
  # Combine positive and negative loadings for plotting
  top_loadings <- rbind(top_pos_loadings, top_neg_loadings)
  
  # Step 4: Plot the top positive and negative contributing genes
  plot <- ggplot(top_loadings, aes(x = reorder(Gene, Loading), y = Loading, fill = Loading > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # Flip coordinates for better visibility
    theme_minimal() +
    scale_fill_manual(values = c("blue", "red"), labels = c("Negative", "Positive")) +
    labs(title = paste0("PC", pc_number),
         x = "Cluster",
         y = "Loading Value",
         fill = "Loading Direction")
  
  return(plot)
}

# PCs to visualize
pcs_to_plot <- c(9, 8, 6, 25, 15)

# Generate the plots for each PC
plots <- lapply(pcs_to_plot, function(pc) plot_top_pos_neg_contributors(pca_result, pc))

# Arrange all plots in a grid for easier comparison (adjust ncol as needed)
combined_plot <- grid.arrange(grobs = plots, ncol = 2)

ggsave("/Volumes/PGH-Backup/ibd_data/shotgun_rnaseq_maaslin/figures/rectum_mgx_pca_loadings.png",
       plot = combined_plot, width = 10, height = 8, units = "in", dpi = 600)
















