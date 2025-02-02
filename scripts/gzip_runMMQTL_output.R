## Kailash BP, 2025
## gzip runMMQTL test_peak_1 output files

# Arguments
library(optparse)
library(tidyverse)

# Define options
option_list <- list(
  make_option(c('--chunk_meta'), help = 'chunk meta file', default = "results/example/example"),
  make_option(c('--eQTL_number'), help = 'number of eQTL peaks', default = 1),
  make_option(c('--output_file'), help = 'chunk output success dummy file', default = "results/example/example")
)

# Parse arguments
option.parser <- OptionParser(option_list = option_list)
opt <- parse_args(option.parser)

# Assign parsed arguments to variables
chunk_meta <- opt$chunk_meta
eQTL_number <- as.numeric(opt$eQTL_number)
output_file <- opt$output_file
mmQTL_tmp_folder <- dirname(chunk_meta)
fname <- basename(chunk_meta)

# Debug messages
message(" * Starting gzipping runMMQTL output files ...")
message("   - Chunk meta file: ", chunk_meta)
message("   - Temp folder: ", mmQTL_tmp_folder)
message("   - Output dummy file: ", output_file)

# Validate eQTL_number
if (eQTL_number < 1) {
  stop("Error: eQTL_number must be at least 1.")
}

# Read the chunk meta file
chunk_meta <- read_tsv(chunk_meta, col_names = c("chr", "start", "end", "feature"))

# Check if features are present
if (nrow(chunk_meta) == 0) {
  message(" * No features found in the chunk meta file. Exiting.")
  writeLines("No features to gzip. Process completed.", output_file)
  quit(status = 0)
}

# Debug: Check chunk_meta content
message(" * Number of features to process: ", nrow(chunk_meta))

# Function to find and gzip relevant files for each feature
gzip_files <- function(feature) {
  # Loop through each statistical signal index
  for (i in seq_len(eQTL_number)) {
    # Dynamically construct file paths
    stat_file <- file.path(mmQTL_tmp_folder, feature, paste0("test_peak_", i, "_statistical_signal"))
    summary_file <- file.path(mmQTL_tmp_folder, feature, paste0("test_peak_", i, "_statistical_signal_summary_correlation"))
    
    # Check and gzip the stat_file if it exists
    if (file.exists(stat_file)) {
      message("   - Gzipping: ", stat_file)
      system(paste("gzip", shQuote(stat_file)))
    } else {
      message("   - WARNING: File not found: ", stat_file)
    }
    
    # Check and gzip the summary_file if it exists
    if (file.exists(summary_file)) {
      message("   - Gzipping: ", summary_file)
      system(paste("gzip", shQuote(summary_file)))
    } else {
      message("   - WARNING: File not found: ", summary_file)
    }
  }
}


# Apply gzipping for each feature
walk(chunk_meta$feature, gzip_files)

# Write a dummy success file
writeLines("Successfully gzipped mmQTL output files", output_file)
message(" * Gzipping completed successfully for ", fname)