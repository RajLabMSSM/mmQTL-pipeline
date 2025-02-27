
## Prepare phenotypes for QTL mapping
# Beomjin Jang
# 2023
library(dplyr)
library(data.table)
library(tidyverse)
library(matrixStats)
library(optparse)
library(edgeR)


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
    make_option(c('--threshold'), help = 'minimum threshold', default = 10),
    make_option(c('--fraction'), help = 'minimum fraction of samples', default = 0),
    make_option(c('--prefix'), help = 'stem of out file', default = "results/example/example"),
    make_option(c('--counts'), help = "path to count matrix", default = "")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

# path_dir <- '/sc/arion/projects/bigbrain/data/ROSMAP/analysis/mmQTL-pipeline/'
# sample_key <- '/sc/arion/projects/bigbrain/data/ROSMAP/analysis/QTL-mapping-pipeline/input_files/ROSMAP_sample_key_EUR_Microglia.txt'
# pheno_matrix <- '/sc/arion/projects/bigbrain/data/ROSMAP/analysis/rds_cell/rosmap_Microglia_RNA_mean.tsv'
# pheno_meta <- '/sc/arion/projects/bigbrain/gencode.v38.annotation_gene_pheno_meta.tsv'
# group <- 'FALSE'
# min_threshold <- 10
# min_fraction <- 0
# counts <- '/sc/arion/projects/bigbrain/data/ROSMAP/analysis/rds_cell/rosmap_Microglia_End_counts.tsv'
# prefix <- 

sample_key <- opt$key
pheno_matrix <- opt$pheno_matrix
pheno_meta <- opt$pheno_meta
group <- opt$group
min_threshold <- opt$threshold
min_fraction <- opt$fraction
min_fraction <- 2
counts <- opt$counts
if( counts == ""){ counts <- NA}
prefix <- opt$prefix

# Processing GTF
# message(" * Processing GTF ")
# annotation_gtf <- '/sc/arion/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v30.annotation.gtf'
# GENCODE <- rtracklayer::import(annotation_gtf)
# gtf <- as.data.frame(GENCODE)
# gtf <- gtf %>% dplyr::filter(type=='gene') %>% dplyr::select( gene_id, gene_name)
# names(gtf) <- c('gene_id', 'transcript_id')
# 
# rm(GENCODE)

# functions
# read in matrix and meta
message(" * Load Meta and Sample Key ")
meta <- read_tsv(pheno_meta)
sk <- read_tsv(sample_key)

# process metadata
# remove features from chromosomes not present in VCF - TODO
message(" * process metadata ")
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

message(" * Load Pheno ")
pheno <- process_pheno(pheno_matrix)
stopifnot(nrow(pheno) >= 0)

message(" * Load Counts ")
counts_df <- process_pheno(counts)
stopifnot(nrow(counts_df) >= 0)


message(" * ", nrow(pheno), " features present in metadata" )


# apply missingness thresholds
# right now assume non-grouped phenotypes are gene expression
message(" * Filtering & Normalization ")
# Gene Expression or Transcript Expression
if( group == FALSE ){
    remove.list =''
    for (i in 1:dim(counts_df)[1]) {
        if (rowCounts(counts_df[i,1:ncol(counts_df)]>0) < as.numeric(min_threshold)) {
            remove.list <- append(remove.list, i)
        }
    }
    remove.list <- remove.list[-1]
    
    dim(counts_df)
    length(remove.list)
    print(remove.list)
    counts_df <- counts_df[-c(as.numeric(remove.list)),]
    stopifnot(nrow(counts_df) >= 0)
    #genes <- rownames(counts_df)
    #length(genes)
    dim(pheno)
    
    pheno <- pheno[-c(as.numeric(remove.list)),]
    #rownames(pheno) <- genes
    print(dim(pheno))
    message(" * ", nrow(counts_df), " features pass missingness thresholds" )
    
    print("CPM Normalization & Invsersetranformation")
    y <- DGEList(counts = counts_df)
    # keep <- filterByExpr(y)
    # y <- y[keep,,keep.lib.sizes=F]
    
    message(" * ", nrow(y$counts), " features pass missingness thresholds" )
    
    y <- calcNormFactors(y, method = "TMM")
    v <- voom(y)
    
    logcpm <- v$E
    
    # pngfile <- sprintf("%s/%s.hist.png", path_dir, prefix)
    # png(pngfile)
    # hist(apply(logcpm, 1, mean), 100)
    # dev.off()
    
    mean_logcpm <- apply(logcpm, 1, mean)
    logcpm <- logcpm[mean_logcpm > 2.0,]
    
    pheno <- logcpm %>%as.data.frame()
    meta <- meta[ row.names(pheno), ]
    message(" * ", nrow(pheno), " features pass missingness thresholds" )

}

# grouped features - eg transcript usage
# apply low expression filter
# remove any singleton features  
# impute missing values 
if( group == TRUE){
    message(" * grouping ")
    # filter low expression
    stop("Group TRUE is not supported in sc-eQTL")

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

