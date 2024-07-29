## Prepare phenotype matrix for edQTL mapping
# Winston H Cuddleston
# 2024
library(tidyverse)
library(optparse)

## Arguments
# - key - data frame matching sample_id (phenotype) to participant_id (genotype)
# - pheno_matrix - a matrix of samples by features (editing sites)
# - pheno_meta - a data frame of feature ID, chr, start and end, this is made by edqtl_collate_phenotype_metadata.R
# - fraction - decimal, what proportion of samples should the coverage threshold be greater than (0.5 default)
# - threshold - a number, what coverage threshold to exclude features on (editing ratio 0.05 default)

option_list <- list(
    make_option(c('--key'), help = 'sample key', default = ''),
    make_option(c('--pheno_matrix'), help = 'editing ratio matrix', default = ''),
    make_option(c('--pheno_meta'), help = 'first version of phenotype metadata', default = ''),
    make_option(c('--threshold'), help = 'minimum editing ratio threshold', default = 0.05),
    make_option(c('--fraction'), help = 'minimum proportion of samples threshold is met in', default = 0.5),
    make_option(c('--prefix'), help = 'stem of out file', default = '')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

sample_key <- opt$key
pheno_matrix <- opt$pheno_matrix
pheno_meta <- opt$pheno_meta
min_threshold <- opt$threshold
min_fraction <- opt$fraction
prefix <- opt$prefix

# read in matrix and meta
meta <- read_tsv(pheno_meta)
sk <- read_tsv(sample_key)

# process metadata
chroms <- paste0("chr", 1:22)
meta <- meta[ meta$chr %in% chroms,]
row.names(meta) <- meta$feature 
message( " * ", nrow(meta), " autosomal features kept " )

# process raw phenotypes to matrix prepared for normalization
process_pheno <- function(mat){
    mat <- read_tsv(mat)
    message(" * read in ", nrow(mat), " features from ", ncol(mat), " samples ")
    # subset out samples only present in sample_key
    stopifnot( all( sk$sample_id %in% names(mat) ) )
    mat <- mat[, c(names(mat)[1], sk$sample_id) ]
    message(" * ", ncol(mat), " samples kept from sample key " )
    # rename columns from samples to donors (participant_id in sample key)
    names(mat) <- c("feature", sk$participant_id )
    mat <- column_to_rownames(mat, var = "feature" )
    # subset out features found in metadata and reorder
    meta_loc <- meta %>% filter(feature %in% row.names(mat)) 
    mat <- mat[ meta_loc$feature, ]
    return(mat)
}

pheno <- process_pheno(pheno_matrix)
message(" * ", nrow(pheno), " features present in metadata" )

# filter out lowly edited sites with high missingness - as defined by flags:
## threshold (editing ratio cutoff, default 5%) and
## fraction (how many samples editing ratio must be met in, default 50%)

# messages along the way help gather info about how many editing sites were dropped at various points

message(" * filtering RNA editing sites ")
features_clean <- rowSums(pheno > min_threshold, na.rm = TRUE) >= min_fraction * ncol(pheno)
pheno <- pheno[ features_clean, ]
meta <- meta[ row.names(pheno), ]
message(" * ", nrow(pheno), " features pass thresholds" ) # how many editing sites had editing rate >5% in >50% of samples?

miss_rate <- sum( is.na(pheno) ) / (nrow(pheno) * ncol(pheno) )
message( " * missing data rate: ", miss_rate ) # how many NAs are still left in the matrix with the default parameters?

# NA values occur when there is no coverage of that editing site in a sample
# mmQTL can't handle 0s so we use mean imputation
# only impute rows with < 25% missing data
missing_data <- rowSums(is.na(pheno) ) <= 0.25 * ncol(pheno)
message(" * removing ", sum(!missing_data), " rows with > 25% missingness") # how many editing sites had to be dropped before mean imputation because >=75% of samples were NA?
pheno <- pheno[ missing_data,]
na_entries <- which(is.na(pheno), arr.ind=TRUE)
pheno[na_entries] <- rowMeans(pheno, na.rm=TRUE)[na_entries[,"row"]]

# we also remove sites with 0 variance - not relevant for QTL mapping
row_sd <- apply(pheno, MARGIN = 1, sd)
message(" * removing ", sum(row_sd == 0), " features with 0 SD" ) # how many editing sites had 0 variance across samples?
pheno <- pheno[ row_sd > 0,] 
meta <- meta[ row.names(pheno), ]
message( " * ", nrow(pheno), " features kept") # how many editing sites are finally left in this cohort?


# scale and center to mean=0 & SD=1 
message(" * scaling and centering ")
pheno <- as.data.frame(t(scale(t(pheno) )))

# quantile normalise
message(" * quantile normalize individuals" )
pheno_q <- as.data.frame(preprocessCore::normalize.quantiles(data.matrix(pheno)))
colnames(pheno_q) <- colnames(pheno)
row.names(pheno_q) <- row.names(pheno)
pheno <- pheno_q

# write out phenotype and metadata
out_file <- paste0(prefix, ".edqtl.phenotype_matrix.tsv.gz")

pheno <- rownames_to_column(pheno,var = "feature")
message(" * writing out")
write_tsv(pheno, file = out_file)

