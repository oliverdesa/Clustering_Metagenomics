library(ggmsa)
library(ggplot2)
library(Biostrings)


alignment_dir <- "/Users/odesa/Desktop/transient_files/aligned_fastas"

# Get a list of all FASTA files in the directory
alignment_files <- list.files(alignment_dir, pattern = "\\.faa$", full.names = TRUE)

# Loop through each alignment file
for (alignment_file in alignment_files) {
  # Extract the file name without extension for the output file
  file_name <- tools::file_path_sans_ext(basename(alignment_file))
  output_file <- file.path(alignment_dir, paste0(file_name, "_msa.png"))
  
  # Print message
  cat("Processing file:", alignment_file, "...\n")
  
  # Read the alignment file to determine the sequence length
  alignment <- Biostrings::readAAStringSet(alignment_file)  # Use readAAStringSet for protein sequences
  seq_length <- max(width(alignment))  # Determine the length of the longest sequence
  
  # Visualize the entire MSA and save to a PNG file
  msa_plot <- ggmsa(alignment_file, start = 1, end = seq_length, char_width = 0.5, color = "Clustal")
  ggsave(filename = output_file, plot = msa_plot, width = 10, height = 8, dpi = 600)
  
  # Print completion message
  cat("Saved MSA plot to:", output_file, "\n")
}
