## Jack Humphrey, Kailash BP 2025
## Collate top associations per peak

library(optparse)
library(tidyverse)

option_list <- list(
   make_option(c('--output_file', '-o'), help = 'name of out file', default = "results/example/example"),
   make_option(c('--prefix', '-p'), help = 'prefix of chr specific top assoc files', default = "results/example/example"),
   make_option(c('--data_key'), help = 'TSV with dataset names used for column renaming', default = ""),
   make_option(c('--QTL_number'), help = 'Number of QTL peaks', default = 1, type = "integer"),
   make_option(c('--QTL_type'), help = 'cis or trans', default = "cis")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser, positional_arguments = TRUE)
output_file <- opt$options$output_file
prefix <- opt$options$prefix
peak <- as.numeric(opt$options$QTL_number)
top_pattern <- paste0("^chr[0-9]+_", "peak_", peak, "_top_assoc\\.tsv\\.gz$")
QTL_type <- opt$options$QTL_type

if (!(QTL_type %in% c("cis", "trans"))) {
  stop("QTL_type must be either 'cis' or 'trans'")
}

print(prefix)
inputs <- list.files(prefix, pattern = top_pattern, recursive = FALSE, full.names = TRUE)

message(" * Processing top associations for QTL peak ", peak)

# Check for input files
if (length(inputs) == 0) {
  message(" * WARNING: No top files found for peak ", peak, ". Writing an empty file.")
  write_tsv(tibble(), output_file)  # Write an empty output file and exit
  quit(save = "no")
}

# Read and combine input files
res <- map_df(inputs, read_tsv)

# Read the data_key file
data_key <- read_tsv(opt$options$data_key)

rename_columns <- function(res, data_key){
  cohorts <- as.list(data_key$dataset)
  naming <- function(x){
    y <- paste0(x,c(".beta",".sd",".z"))
    return(y)
  }
  cohort_names <- unlist(lapply(cohorts, naming))
  res_names <- c("feature","variant_id","chr","pos","ref","alt","Allele", cohort_names,
                 "fixed_beta","fixed_sd","fixed_z","Random_Z","Fixed_P","Random_P",
                 "Fixed_bonf","Random_bonf","Fixed_FDR","Random_FDR","qval")
  
  colnames(res) <- res_names
  return(res)
  
}

res <- rename_columns(res, data_key)

if (QTL_type == "cis") {
   if (nrow(res) > 1) {
     res$Random_Bonf_FDR <- p.adjust(res$Random_bonf, method = "fdr")
   } else {
     res$Random_Bonf_FDR <- NA
   }


   # Sort by Random_Bonf_FDR and write output
   res <- arrange(res, Random_Bonf_FDR)
   write_tsv(res, output_file)
}

if (QTL_type == "trans") {
   res <- res %>% group_by(feature) %>% slice_min(Random_P, with_ties = FALSE) %>% ungroup() %>% arrange(Random_P)
   write_tsv(res, output_file)
}

message(" * Top associations for peak ", peak, " written to: ", output_file)
message(" * Collation complete.")
