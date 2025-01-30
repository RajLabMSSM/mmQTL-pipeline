## Kailash BP, 2025
## Make the trans phenotype metadata file

# Output will have each gene show up 22 times as gene_chr1, gene_chr2, ..., gene_chr22
# The way trans is being run is cis-QTL per chromosome with every feature being as long as the chromosome itself

# arguments

library(optparse)
library(purrr)
library(tidyverse)

# Define options
option_list <- list(
  make_option(c('--outFolder'), help = 'folder to write results to', default = "results/example/example"),
  make_option(c('--pheno_meta'), help = 'path to phenotype metadata file of all features', default = ""),
  make_option(c('--pheno_meta_harmonized'), help = 'path to harmonized phenotype metadata file of all features', default = ""),
  make_option(c('--pheno_meta_trans_to_test'), help = 'path to phenotype metadata file of all features to test', default = "")
)

# Parse arguments
option.parser <- OptionParser(option_list = option_list)
opt <- parse_args(option.parser)

# Debug: Print parsed options
print(opt)

# Debug: Print parsed options
print(names(opt))

# Assign parsed arguments to variables
outFolder <- opt$outFolder
pheno_meta <- opt$pheno_meta
pheno_meta_trans_to_test <- opt$pheno_meta_trans_to_test
pheno_meta_harmonized <- opt$pheno_meta_harmonized

# Debug: Print arguments
message(" * Processing trans phenotype metadata files ...")
message("   - Outfolder: ", outFolder)
message("   - Original phenotype metadata file: ", pheno_meta)
message("   - Harmonized phenotype metadata file: ", pheno_meta_harmonized)
message("   - Phenotype metadata file of features to test for trans: ", pheno_meta_trans_to_test)

pheno_meta <- read_tsv(pheno_meta) # To get max CHR START END 

pheno_meta_harmonized <- read_tsv(pheno_meta_harmonized, col_names = F) # To test only for harmonized features
colnames(pheno_meta_harmonized) <- c("chr", "start", "end", "feature")

pheno_meta_trans_to_test <- read_tsv(pheno_meta_trans_to_test) %>% # To test only for selected features
  filter(feature %in% pheno_meta_harmonized$feature)

pheno_meta_trans_to_test <- pheno_meta_trans_to_test %>%
  dplyr::rename("feature_start" = "start", "feature_end" = "end")

# Retaining only autosomal features
chromosomes <- sapply(1:22, function(chr){return(paste0("chr", chr))})

chr_max_start_end <- pheno_meta %>%
  group_by(chr) %>%
  mutate(start = min(start), end = max(end)) %>%
  distinct(chr, start, end) %>%
  dplyr::rename("chr_start" = "start", "chr_end" = "end") %>% 
  filter(chr %in% chromosomes)

# Processing chromosomes for phenotype metadata
map(chromosomes, ~ {
  pheno_meta_trans_to_test %>%
    mutate(chr = .x, feature = paste0(feature, "_", .x)) %>%
    left_join(chr_max_start_end, by = "chr") %>%
    dplyr::rename(start = chr_start, end = chr_end) %>%
    select(chr, start, end, feature) %>% 
    write_tsv(paste0(outFolder, "pheno_meta_", .x, ".trans.genes.tsv"))
})

phen_mdata_chr_files <- list.files(path = outFolder, pattern = "pheno_meta_chr([1-9]|1[0-9]|2[0-2])", full.names = TRUE)
# Merging all chromosome-specific phenotype metadata into one file
map_dfr(phen_mdata_chr_files, read_tsv) %>%
  write_tsv(paste0(outFolder, "phenotype_metadata_trans.tsv"), col_names = FALSE)

# Check if files exist before attempting to remove
if (length(phen_mdata_chr_files) > 0) {
  file.remove(phen_mdata_chr_files)
} else {
  message("No matching files found to remove.")
}
