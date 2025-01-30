## Kailash BP, 2025
## Make the trans phenotype matrix

# The trans phenotype matrix copies each row of gene in harmonised.tsv file into 22 rows
# - redundant but necessary to test for QTLs for each feature in each chromosome

# arguments

library(optparse)
library(purrr)
library(dplyr)
library(readr)
library(tidyverse)

# Define options
option_list <- list(
  make_option(c('--pheno'), help = 'path to input harmonised phenotype files', default = ""),
  make_option(c('--pheno_meta_trans'), help = 'path to phenotype metadata file of trans features to test', default = ""),
  make_option(c('--outFolder'), help = 'folder to write results to', default = "results/example/example")
)

# Parse arguments
option.parser <- OptionParser(option_list = option_list)
opt <- parse_args(option.parser)

# Debug: Print options
print(opt)

# Assign parsed arguments to variables
pheno <- opt$pheno
pheno_meta_trans <- opt$pheno_meta_trans
outFolder <- opt$outFolder

outFileName <- str_replace(basename(pheno), ".tsv" , ".trans.tsv")
            
# Debug: Print arguments
message(" * Processing trans phenotypes...")
message("   - Harmonized phenotype file: ", pheno)
message("   - Trans phenotype-metadata file: ", pheno_meta_trans)
message("   - Output file name: ", outFileName)
message("   - Output folder: ", outFolder)

pheno_meta_trans <- read_tsv(pheno_meta_trans)

chromosomes <- sapply(1:22, function(chr){return(paste0("chr", chr))})

# Column 1 is feature name in pheno file
# Column 4 is feature name in phenotype metadata file
pheno_file <- read_tsv(pheno) %>%
  filter(.[[1]] %in% pheno_meta_trans$feature)

map_dfr(chromosomes, ~ pheno_file %>%
          mutate("{names(.)[1]}" := paste0(.data[[names(.)[1]]], "_", .x))) %>%
  write_tsv(paste0(outFolder, (outFileName)))
