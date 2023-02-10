
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
    meta_loc <- filter(meta,feature %in% row.names(mat) ) 
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
# right now assume non-grouped phenotypes are gene expression
if( group == FALSE ){
    
    features_clean <- rowSums(pheno >= min_threshold) > min_fraction * ncol(pheno)

if( is.na(counts) ){
    message(" * using TPM matrix")
    # use TPMs as phenotype. Filter and log normalise
    pheno <- pheno[ features_clean, ]
    meta <- meta[ row.names(pheno), ]
    # log normalise matrix
    pheno <- log2(pheno + 1) 

    }else{
    message(" * using count matrix" )
    # use counts as phenotype. Filter on TPM expression
    pheno <- counts_df[ features_clean,] 
    stopifnot(nrow(pheno) >= 0) 
    print(head(pheno) )
    save.image("debug.RData")
    # then perform TMM normalisation
    # round counts if estimated
    pheno <- floor(pheno)
    norm <- edgeR::calcNormFactors(pheno, method = "TMM") #normalized counts with TMM
    dge <- edgeR::DGEList(counts=pheno, norm.factors = norm)  
    v <- limma::voom(dge)
    pheno <- as.data.frame(v$E)
    }
message(" * ", nrow(pheno), " features pass missingness thresholds" )

}

# grouped features - eg transcripts
# apply low expression filter
# remove any singleton features  
# impute missing values 
if( group == TRUE){
    message(" * grouping ")
    # filter low expression
    features_clean <- rowSums(pheno > min_threshold) >= min_fraction * ncol(pheno)
    pheno <- pheno[ features_clean, ]
    meta <- meta[ row.names(pheno), ]
    message(" * ", nrow(pheno), " features pass expression thresholds" )
 
    # tally groups in meta - any group that now only appears once should be removed from pheno matrix and metadata
    group_tally <- group_by(meta, group) %>% tally()
    singletons <- group_tally$group[ group_tally$n == 1]
    meta <- meta[ !meta$group %in% singletons, ] 
    pheno <- pheno[meta$feature,]    
    message( " * ", length(singletons), " singleton groups removed" )

    pheno_split <- split(pheno, meta$group )
    pheno <- map_df( pheno_split, ~{
            df <- sweep(.x, MARGIN = 2, STATS =  colSums(.x), FUN = "/")
            return(as.data.frame(df))
    })
    miss_rate <- sum( is.na(pheno) ) / (nrow(pheno) * ncol(pheno) )
    message( " * missing data rate: ", miss_rate )
    # set divide by 0 errors to 0 - bad approach
    # NA values occurs when all transcripts are 0 in a sample - in this case impute with the mean
    # only impute rows with < 10% missing data
    missing_data <- rowSums(is.na(pheno) ) <= 0.1 * ncol(pheno)
    message(" * removing ", sum(!missing_data), " rows with > 10% missingness")
    pheno <- pheno[ missing_data,]
    # then impute missing entries with row mean
    na_entries <- which(is.na(pheno), arr.ind=TRUE)
    pheno[na_entries] <- rowMeans(pheno, na.rm=TRUE)[na_entries[,"row"]]

    # a few weird entries where transcripts are duplicated - all PSI are 0.5
    row_sd <- apply(pheno, MARGIN = 1, sd)
    message(" * removing ", sum(row_sd == 0), " features with 0 SD" )
    pheno <- pheno[ row_sd > 0,] 
    meta <- meta[ row.names(pheno), ]
    message( " * ", nrow(pheno), " features kept") 
}

message(" * scaling and centering ")
# scale and centre to means of 0 and SD of 1 
pheno <- as.data.frame(t(scale(t(pheno) )))

#message(" * ", nrow(pheno), " features pass missingness thresholds" )
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

