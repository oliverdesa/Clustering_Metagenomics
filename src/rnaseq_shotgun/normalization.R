library(DESeq2)
library(tidyverse)
library(edgeR)

rnaseq <- readr::read_tsv("/Volumes/PGH-Backup/ibd_data/rnaseq/host_tx_counts.tsv")

metadata <- read.csv("/Volumes/PGH-Backup/ibd_data/metadata/hmp2_metadata_2018-08-20.csv")

metadata_rnaseq <- metadata %>%
  filter(data_type == "host_transcriptomics")

# Set patient IDs as row names
metadata_rnaseq <- metadata_rnaseq %>% column_to_rownames(var = "External.ID")

# Set transcripts as rownames
rnaseq <- rnaseq %>% column_to_rownames(var = "Transcript")

common_samples <- intersect(colnames(rnaseq), rownames(metadata_rnaseq))

rnaseq_filtered <- rnaseq[, common_samples]
metadata_rnaseq_filtered <- metadata_rnaseq[common_samples, ]

dim(rnaseq_filtered)
dim(metadata_rnaseq_filtered)

# Create DESeq2 dataset
# control for repeated measures, diagnosis and lcoation of biopsy
# Issue is some of these are linear, address later
# dds <- DESeqDataSetFromMatrix(countData = rnaseq_filtered,
#                               colData = metadata_rnaseq_filtered,
#                               design = ~ Participant.ID + diagnosis + biopsy_location)

dds <- DESeqDataSetFromMatrix(countData = rnaseq_filtered,
                              colData = metadata_rnaseq_filtered,
                              design = ~ 1)

# Check library sizes (total counts per sample)
library_sizes <- colSums(counts(dds))

# Identify samples with zero total counts
zero_count_samples <- names(library_sizes[library_sizes == 0])
print(zero_count_samples)

# If there are any zero-count samples, remove them
if (length(zero_count_samples) > 0) {
  # Remove zero-count samples from the count matrix
  rnaseq_filtered_nonzero <- rnaseq_filtered[, !(colnames(rnaseq_filtered) %in% zero_count_samples)]
  
  # Remove zero-count samples from the metadata
  metadata_rnaseq_filtered_nonzero <- metadata_rnaseq_filtered[!(rownames(metadata_rnaseq_filtered) %in% zero_count_samples), ]
  
  # Update the DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(countData = rnaseq_filtered_nonzero,
                                colData = metadata_rnaseq_filtered_nonzero,
                                design = ~ 1)
} else {
  rnaseq_filtered_nonzero <- rnaseq_filtered
  metadata_rnaseq_filtered_nonzero <- metadata_rnaseq_filtered
}

# LEts try EdgeR
dge <- DGEList(counts = counts(dds))

keep <- filterByExpr(dge, group = metadata_rnaseq_filtered$diagnosis)
dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]

# TMM Normalization
dge_filtered <- calcNormFactors(dge_filtered)

## Write normalized counts tabel to csv

# Extract normalized counts as a matrix
normalized_counts <- cpm(dge_filtered, normalized.lib.sizes = TRUE)

# Convert to data frame and add gene names as a column
normalized_counts_df <- as.data.frame(normalized_counts)
normalized_counts_df$Gene <- rownames(normalized_counts_df)

# Write the normalized counts to a TSV file
write.table(normalized_counts_df, "/Volumes/PGH-Backup/ibd_data/normalized_counts.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

write.table(metadata_rnaseq_filtered, "/Volumes/PGH-Backup/ibd_data/hmp2_rnaseq_metadata.tsv", sep = "\t", row.names = FALSE)



