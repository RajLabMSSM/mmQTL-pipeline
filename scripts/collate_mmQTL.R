## Collate MMQTL results
## Kailash BP, Jack Humphrey 2025

## for each feature in a chromosome
## read in result file
## match on SNP coordinates and sort 
## convert Random_Z to P-value
## Bonferroni adjust
## write top association to a file
## write sorted associations to a file - no header
## concatenate all the pieces into a single chromosome file

library(optparse)

option_list <- list(
   make_option(c('--prefix'), help = 'stem of out file', default = "results/example/example"),
   make_option(c('--chrom'), help = 'the chromosome', default = ""),
   make_option(c('--metadata'), help = 'phenotype metadata file', default = ""),
   make_option(c('--geno'), help = 'path to genotype folder', default = ""),
   make_option(c('--eQTL_number'), help = 'Number of eQTL peaks', default = 1),
   make_option(c('--QTL_type'), help = 'cis or trans', default = "cis")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

print(opt)

prefix <- opt$prefix
chrom <- opt$chrom
pheno_meta <- opt$metadata
peak <- as.numeric(opt$eQTL_number)
QTL_type <- opt$QTL_type
geno_folder <- opt$geno # Get genotype SNP coordinates
top_file <- paste0(prefix, chrom, "_peak_", peak, "_top_assoc.tsv.gz")
print(paste0("Output file path : ", top_file))

library(readr)
library(purrr)
library(dplyr)

message(" * Collating for chromosome: ", chrom, ", eQTL peak: ", peak)

meta <- read_tsv(pheno_meta, col_names = FALSE)
names(meta)[1:4] <- c("chr", "start", "end", "feature") 
meta_loc <- meta[ meta$chr %in% chrom, ]

## helper functions
add_P <- function(data){
  data$Fixed_P <- 2*pnorm(q=abs(data$fixed_z), lower.tail=FALSE)
  data$Random_P <- 2*pnorm(q=abs(data$Random_Z), lower.tail=FALSE)
  data$Fixed_bonf <- p.adjust(data$Fixed_P, method = "bonferroni")
  data$Random_bonf <- p.adjust(data$Random_P, method = "bonferroni")
  data$Fixed_FDR <- p.adjust(data$Fixed_P, method = "BH" )
  data$Random_FDR <- p.adjust(data$Random_P, method = "BH" )
  return(data)
}

remove_files <- function(file_list) {for (file in file_list) {
  # Check for the original file and its corresponding summary correlation file
  summary_file <- sub("statistical_signal\\.gz$", "statistical_signal_summary_correlation.gz", file)
  # Check and remove the original file
  if (file.exists(file)) {file.remove(file)}
  # Check and remove the corresponding summary correlation file
  if (file.exists(summary_file)) {file.remove(summary_file)}
}
}

# Dynamically generate pattern for current peak
peak_pattern <- paste0("test_peak_", peak, "_statistical_signal\\.gz$")
message("Processing files for peak: ", peak)

if (chrom %in% meta$chr == FALSE) {
  message("WARNING: Chromosome ", chrom, " not found in metadata. Exiting.")
  file.create(top_file)  # Write empty file and exit
  quit(save = "no") # quit R without saving workspace
}

all_files <- list.files(prefix, pattern = peak_pattern, recursive = TRUE, full.names = FALSE)
files_loc <- all_files[dirname(all_files) %in% meta_loc$feature]

if (length(files_loc) == 0) {
  message("WARNING: No files found for peak ", peak, " on chromosome ", chrom)
  file.create(top_file)  # Write empty file and exit
  next
}

features_loc <- dirname(files_loc)
files_loc <- paste0(prefix, files_loc) # Full file paths now
names(files_loc) <- features_loc

message(" * ", length(files_loc), " files found for peak ", peak)

## get SNP coordinates from the coordinate BIM files for each dataset
bim_files <- list.files(geno_folder, pattern = paste0("*", chrom, ".bim"), recursive = TRUE, full.names = TRUE)
## read in and get distinct rows
bim_res <- map_df(bim_files, read_tsv, col_names = c("chr", "Variant", "n", "pos", "alt", "ref")) %>% distinct()
bim_res$chr <- paste0("chr", bim_res$chr)
bim_res$n <- NULL

# get number of datasets
n_datasets <- length(bim_files)

# some features do not use all datasets in their meta-analysis
# and thus have missing columns which messes up concatenation
## prepare columns according to the number of datasets
data_cols <- paste0(rep(c("beta_", "sd_", "z_"), n_datasets), rep(paste0("tissue_", 0:(n_datasets - 1)), each = 3))
all_cols <- c("Variant", "Allele", data_cols, "fixed_beta", "fixed_sd", "fixed_z", "Random_Z")
dummy_cols <- setNames(data.frame(matrix(ncol = length(all_cols), nrow = 1)), all_cols)

top_assoc <- list()

# Process each feature
for (feature in features_loc) {
  message(" * Processing feature: ", feature)
  
  d <- read_tsv(files_loc[feature], col_types = list(Allele = "c"))
  
  # ignore empty or malformed files
  if (nrow(d) == 0 || !"fixed_z" %in% names(d)) next
  
  # standardise columns
  d <- bind_rows(dummy_cols, d)[2:nrow(d), ]
  d$feature <- feature
  # add in SNP info
  d <- d %>% left_join(bim_res, by = "Variant") %>% arrange(pos) %>%
    select(feature, variant_id = Variant, chr, pos, ref, alt, everything())
  
  # add P values
  d <- add_P(d)
  
  # Write per-feature nominal associations
  out_file <- paste0(prefix, chrom, "_", feature, "_peak_", peak, "_all_nominal.tsv.gz")
  
  if (QTL_type == "cis") {
    write_tsv(d, out_file, col_names = FALSE)
  } else if (QTL_type == "trans") {
    write_tsv(d %>% select(feature, variant_id, chr, pos, ref, alt, Random_P, Random_Z), out_file, col_names = FALSE)
  }
 
  # Store top association
  top <- arrange(d, Random_P) %>% head(1)
  top_assoc[[feature]] <- top
}

# Combine and write top associations for current peak
if (length(top_assoc) == 0) {
  message("No top associations found for peak ", peak)
  file.create(top_file)
} else {
  # Combine top associations and write to file
  top_res <- bind_rows(top_assoc)
  if (QTL_type == "cis") {
    write_tsv(top_res, top_file)
  } else if (QTL_type == "trans") {
    write_tsv(top_res %>% select(feature, variant_id, chr, pos, ref, alt, Random_P, Random_Z), top_file)
  }
}

# going to keep these files for now.
# if (length(files_loc) > 0) {
#   remove_files(files_loc)
# }

message(" * Collation complete for chromosome: ", chrom, ", peak: ", peak)
