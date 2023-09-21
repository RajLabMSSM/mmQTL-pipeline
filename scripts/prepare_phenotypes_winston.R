## Prepare phenotypes for QTL mapping
# Jack Humphrey
# 2021
library(tidyverse)
library(optparse)

# arguments

# - sample_key - data frame matching sample_id to participant_id
# - pheno_matrix - a matrix of samples by features, ideally TPM normalised, with feature ID
# - pheno_meta - a data frame of feature ID, group ID, chr, start and end
# - junctions - whether the data is from junctions
# - group - whether to divide features by a grouping id column in the pheno_meta
# - prop_samples - decimal, what proportion of samples should the coverage threshold be greater than (0.5)?
# - cov_threshold - a number, what coverage threshold to exclude features on (TPM 1)
# - counts - if this is added, then use counts as the phenotype and filter on TPM

option_list <- list(
    make_option(c('--key'), help = 'sample key', default = ""),
    make_option(c('--pheno_matrix'), help = 'phenotype count matrix', default = ""),
    make_option(c('--pheno_meta'), help = 'phenotype metadata', default = ""),
    make_option(c('--group'), help = "whether to divide features by a group_id column", default = FALSE, action = "store_true"),
    make_option(c('--edqtl'), help = "RNA editing ratios as phenotypes", default = TRUE, action = "store_true"),
    make_option(c('--threshold'), help = 'minimum threshold', default = 1),
    make_option(c('--fraction'), help = 'minimum fraction of samples', default = 0.5),
    make_option(c('--prefix'), help = 'stem of out file', default = "results/example/example"),
    make_option(c('--counts'), help = "path to count matrix", default = "")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

sample_key <- opt$key
pheno_matrix <- opt$pheno_matrix
pheno_meta <- opt$pheno_meta
group <- opt$group
edqtl <- opt$edqtl
min_threshold <- opt$threshold
min_fraction <- opt$fraction
counts <- opt$counts
if( counts == ""){ counts <- NA}
prefix <- opt$prefix

# functions
# read in matrix and meta

meta <- read_tsv(pheno_meta)
sk <- read_tsv(sample_key)

# process metadata
# remove features from chromosomes not present in VCF - TODO
chroms <- paste0("chr", 1:22)
blacklist <- c("chr2:206122585:A:G","chr2:37104484:A:G")
meta <- meta[ meta$chr %in% chroms,]
row.names(meta) <- meta$feature 
message( " * ", nrow(meta), " autosomal features kept " )


process_pheno <- function(mat){
    mat <- read_tsv(mat)
    message(" * read in ", nrow(mat), " features from ", ncol(mat), " samples ")
    # subset out samples only present in sample_key
    stopifnot( all( sk$sample_id %in% names(mat) ) )
    if( "group" %in% names(mat) ){mat$group <- NULL }
    # if single feature
    mat <- mat[, c(names(mat)[1], sk$sample_id) ]
    message(" * ", ncol(mat), " samples kept from sample key " )
    # rename columns from samples to donors (participant_id in sample key)
    names(mat) <- c("feature", sk$participant_id )
    mat <- column_to_rownames(mat, var = "feature" )
    # subset out features found in metadata and reorder
    meta_loc <- meta %>% filter(feature %in% row.names(mat)) %>% filter(!feature %in% blacklist) 
    mat <- mat[ meta_loc$feature, ]
    return(mat)
}


pheno <- process_pheno(pheno_matrix)

if( !is.na(counts) ) {
    message( " * reading additional counts matrix")
    counts_df <- process_pheno(counts)
}


message(" * ", nrow(pheno), " features present in metadata" )


# apply missingness thresholds
if( edqtl == TRUE){
  message(" * filtering RNA editing sites ")
  # filter lowly edited sites with high missingness
  features_clean <- rowSums(pheno > min_threshold, na.rm = TRUE) >= min_fraction * ncol(pheno) # editing ratio cutoff 5% and missingness cutoff 50%
  pheno <- pheno[ features_clean, ]
  meta <- meta[ row.names(pheno), ]
  message(" * ", nrow(pheno), " features pass thresholds" )

  miss_rate <- sum( is.na(pheno) ) / (nrow(pheno) * ncol(pheno) )
  message( " * missing data rate: ", miss_rate )
  # set divide by 0 errors to 0 - bad approach
  # NA values occur when there is no coverage of the editing site in a sample - impute with the mean
  # only impute rows with < 25% missing data
  missing_data <- rowSums(is.na(pheno) ) <= 0.25 * ncol(pheno)
  message(" * removing ", sum(!missing_data), " rows with > 25% missingness")
  pheno <- pheno[ missing_data,]
  # then impute missing entries with row mean
  na_entries <- which(is.na(pheno), arr.ind=TRUE)
  pheno[na_entries] <- rowMeans(pheno, na.rm=TRUE)[na_entries[,"row"]]
  # remove sites with no variance
  row_sd <- apply(pheno, MARGIN = 1, sd)
  message(" * removing ", sum(row_sd == 0), " features with 0 SD" )
  pheno <- pheno[ row_sd > 0,] 
  meta <- meta[ row.names(pheno), ]
  message( " * ", nrow(pheno), " features kept") 
}

message(" * scaling and centering ")
# scale and centre to means of 0 and SD of 1 
pheno <- as.data.frame(t(scale(t(pheno) )))

# quantile normalise
message(" * quantile normalize individuals" )
pheno_q <- as.data.frame(preprocessCore::normalize.quantiles(data.matrix(pheno)))
colnames(pheno_q) <- colnames(pheno)
row.names(pheno_q) <- row.names(pheno)

pheno <- pheno_q
# write out phenotype and metadata
out_file <- paste0(prefix, "_pheno.tsv.gz")
out_meta <- paste0(prefix, "_meta.tsv")

pheno <- rownames_to_column(pheno,var = "feature")
message(" * writing out")
write_tsv(pheno, file = out_file)
write_tsv(meta, file = out_meta)

