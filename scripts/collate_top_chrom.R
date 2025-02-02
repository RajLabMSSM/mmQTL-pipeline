## Collate top associations for multiple peaks

# args
# --out - the path to output file
# positional arguments for input files

library(optparse)
library(tidyverse)

option_list <- list(
  make_option(c('--output_peak_1', '-o'), help = 'prefix of out file peak 1', default = "results/example/example"),
   make_option(c('--output_prefix', '-o'), help = 'prefix of out file', default = "results/example/example"),
   make_option(c('--eQTL_number'), help = 'Number of eQTL peaks', default = 1, type = "integer")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser, positional_arguments = TRUE)

output_peak_1 <- opt$output_peak_1 # Dummy variable, to ensure input for fullCollate.
output_prefix <- opt$output_prefix

inputs <- opt$args
eQTL_number <- as.numeric(opt$eQTL_number)

# Loop through each eQTL peak and process separately
for (peak in seq_len(eQTL_number)) {
  message("Processing eQTL peak ", peak)
  
  # Filter files corresponding to the current peak
  peak_files <- inputs[str_detect(inputs, paste0("peak_", peak, "_"))]
  
  if (length(peak_files) == 0) {
    message("No input files found for peak ", peak, ". Skipping...")
    next
  }
  
  # Read and combine input files
  res <- map_df(peak_files, read_tsv)
  
  # Calculate q-values using the Random_FDR column
  res$qval <- qvalue::qvalue(res$Random_FDR)$qval
  
  # Sort by q-value
  res <- arrange(res, qval)
  
  # Write output to a separate file for the current peak
  output_file <- paste0(output_prefix, "_peak_", peak, "_top_assoc.tsv")
  write_tsv(res, output_file)
  
  message("Top associations for peak ", peak, " written to: ", output_file)
}
