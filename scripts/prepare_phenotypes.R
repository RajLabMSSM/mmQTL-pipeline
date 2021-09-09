
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

option_list <- list(
    make_option(c('--key'), help = 'sample key', default = ""),
    make_option(c('--pheno_matrix'), help = 'phenotype count matrix', default = ""),
    make_option(c('--pheno_meta'), help = 'phenotype metadata', default = ""),
    make_option(c('--group'), help = "whether to divide features by a group_id column", default = FALSE, action = "store_true"),
    make_option(c('--threshold'), help = 'minimum threshold', default = 1),
    make_option(c('--fraction'), help = 'minimum fraction of samples', default = 0.5),
    make_option(c('--prefix'), help = 'stem of out file', default = "results/example/example")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

sample_key <- opt$key
pheno_matrix <- opt$pheno_matrix
pheno_meta <- opt$pheno_meta
group <- opt$group
min_threshold <- opt$threshold
min_fraction <- opt$fraction

prefix <- opt$prefix

# functions

# read in matrix and meta
pheno <- read_tsv(pheno_matrix)
meta <- read_tsv(pheno_meta)

row.names(meta) <- meta$feature 

message(" * read in ", nrow(pheno), " features from ", ncol(pheno), " samples ")

sk <- read_tsv(sample_key)

# subset out samples only present in sample_key
stopifnot( all( sk$sample_id %in% names(pheno) ) )
pheno <- pheno[, c(names(pheno)[1], sk$sample_id) ]

message(" * ", ncol(pheno), " samples kept from sample key " )


# rename columns from samples to donors (participant_id in sample key)
names(pheno) <- c("feature", sk$participant_id )
#print(head(pheno) )
pheno <- column_to_rownames(pheno, var = "feature" )

# subset out features found in metadata and reorder
pheno <- pheno[ meta$feature, ]

message(" * ", nrow(pheno), " features present in metadata" )

# apply missingness thresholds

features_clean <- rowSums(pheno >= min_threshold) > min_fraction * ncol(pheno)

pheno <- pheno[ features_clean, ]

message(" * ", nrow(pheno), " features pass missingness thresholds" )

meta <- meta[ row.names(pheno), ]

# if requested, divide each feature by the sum of the group (for tuQTLs, SUPPA and txrevise)
if( group == TRUE){
    # 
}

# log normalise matrix
pheno <- log2(pheno + 1) 
# scale and centre to means of 0 and SD of 1 
pheno <- as.data.frame(t(scale(t(pheno) )))

# clean up meta
row.names(meta) <- meta$feature
meta <- meta[ row.names(pheno), ]


# write out phenotype and metadata
out_file <- paste0(prefix, "_pheno.tsv.gz")
out_meta <- paste0(prefix, "_meta.tsv")

pheno <- rownames_to_column(pheno,var = "feature")

write_tsv(pheno, file = out_file)

write_tsv(meta, file = out_meta)

