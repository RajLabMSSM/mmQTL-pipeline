## Jack Humphrey, Kailash BP 2025
## Collate top associations per peak

library(optparse)
library(tidyverse)

option_list <- list(
   make_option(c('--output_file', '-o'), help = 'name of out file', default = "results/example/example"),
   make_option(c('--prefix', '-p'), help = 'prefix of chr specific top assoc files', default = "results/example/example"),
   make_option(c('--eQTL_number'), help = 'Number of eQTL peaks', default = 1, type = "integer")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser, positional_arguments = TRUE)
output_file <- opt$options$output_file
prefix <- opt$options$prefix
peak <- as.numeric(opt$options$eQTL_number)
top_pattern <- paste0("^chr[0-9]+_", "peak_", peak, "_top_assoc\\.tsv\\.gz$")
print(prefix)
inputs <- list.files(prefix, pattern = top_pattern, recursive = FALSE, full.names = TRUE)

message(" * Processing top associations for eQTL peak ", peak)

# Check for input files
if (length(inputs) == 0) {
  message(" * WARNING: No top files found for peak ", peak, ". Writing an empty file.")
  write_tsv(tibble(), output_file)  # Write an empty output file and exit
  quit(save = "no")
}

# Read and combine input files
res <- map_df(inputs, read_tsv)

if (nrow(res) > 1) {
  res$Random_Bonf_FDR <- p.adjust(res$Random_bonf, method = "fdr")
} else {
  res$Random_bonf_fdr <- NA
}

# Sort by Random_Bonf_FDR and write output
res <- arrange(res, Random_Bonf_FDR)
write_tsv(res, output_file)

message(" * Top associations for peak ", peak, " written to: ", output_file)
message(" * Collation complete.")