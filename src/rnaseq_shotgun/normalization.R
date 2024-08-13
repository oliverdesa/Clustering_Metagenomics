library(DESeq2)
library(tidyverse)

rnaseq <- readr::read_tsv("/Volumes/PGH-Backup/ibd_data/rnaseq/host_tx_counts.tsv")

metadata <- read.csv("/Volumes/PGH-Backup/ibd_data/metadata/hmp2_metadata_2018-08-20.csv")

metadata_rnaseq <- metadata %>%
                    filter(data_type == "host_transcriptomics")

dim(rnaseq)
dim(metadata_rnaseq)

colnames(rnaseq)
rownames(metadata_rnaseq)

metadata_rnaseq <- metadata_rnaseq %>% column_to_rownames(var = "External.ID")

rnaseq <- rnaseq %>% column_to_rownames(var = "Transcript")

common_samples <- intersect(colnames(rnaseq), rownames(metadata_rnaseq))

rnaseq_filtered <- rnaseq[, common_samples]

dim(rnaseq_filtered)
dim(metadata_rnaseq)

# Create DESeq2 dataset
# control for repeated measures, diagnosis and lcoation of biopsy
dds <- DESeqDataSetFromMatrix(countData = rnaseq,
                              colData = metadata_rnaseq,
                              design = ~ Participant.ID + diagnosis + biopsy_location)

