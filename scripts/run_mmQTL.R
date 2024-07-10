
options(echo=TRUE)
## RUN MMQTL

# Jack Humphrey
# 2021
# wraps around the MMQTL binary 
# for massive parallelisation

## Arguments:
# all the paths to the input file lists needed
# path to phenotype metadata
# chunk number (n) and chunk number (i of n)
library(optparse)
option_list <- list(
    make_option(c('--chrom'), help = "which chromosome to run on", default = ""),
    make_option(c('--pheno_meta'), help = 'phenotype metadata', default = ""),
    make_option(c('--chunk_total', '-n'), help = 'total number of chunks', default = 10),
    make_option(c('--chunk', '-i'), help = 'current chunk number', default = 1),
    make_option(c('--prefix'), help = 'stem of out file', default = "results/example/example"),
    make_option(c('--geno_file'), help = "path to geno list file" ),
    make_option(c('--pheno_file'), help = "path to pheno list file"),
    make_option(c('--grm_file'), help = "path to GRM list file"),
    make_option(c('--cov_file'), help = "path to covariate_file"),
    make_option(c('--mmQTL'), help = "full path to MMQTL executable", default = "MMQTL25")
)

absPath <- function(path){
    if(is.null(path) ){ return()}
    if(file.exists(path) ){
        normalizePath(path)
    }
}
option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

chrom <- opt$chrom
meta_file <- absPath(opt$pheno_meta)
geno_file <- absPath(opt$geno_file)
pheno_file <- absPath(opt$pheno_file)
cov_file <- absPath(opt$cov_file)
grm_file <- absPath(opt$grm_file)
i_chunk <- as.numeric(opt$chunk)
n_chunk <- as.numeric(opt$chunk_total)
prefix <- absPath(opt$prefix)
mmqtl_bin <- opt$mmQTL

cis_window <- 1e6

if( !dir.exists(prefix) ){
    dir.create(prefix, recursive = TRUE)
}

# Actions:
#save.image("debug.RData")
# split the list of features into a given number of chunks (n)
split_chunks <- function(meta, i, n){
   bins <- dplyr::ntile(1:nrow(meta), n)
   
   chunk <- meta[ bins == i, ] 
   print(chunk)
   return(chunk)
}

# run commands from within outFolder
setwd(prefix)
# get out the i'th chunk of features
# make a system call to MMQTL for each feature
# handle exceptions somehow
# for j'th entry in meta_loc, make system call to MMQTL

run_mmQTL <- function(meta_loc, j){
     
    feature_j <- meta_loc[j, ]$feature

    cmd1 <- paste0( mmqtl_bin,
        " -b ", 
        " -P ", pheno_file,
        " -Z ", geno_file,
        " -R ", grm_file,
        #" -C ", cov_file,
        " -a ", meta_chunk_file,
        " -A random ",
        " -gene ", feature_j,
        " --threads 4" ,
        #" --out ", prefix,## adding this kills the analysis
        " --primary_only ", ## TODO - make optional
        " -V ", format(cis_window, scientific = FALSE, digits = 1),
        " -Han " # Apply Han & Eskin method to adjust results of random-effect model
        # MMQTL24 -b  -P  pheno_file.txt   -Z  geno_file.txt   -R GRM_file.txt -a feature_annotation.bed  -A random   -gene  gene_name 
    )
    print(cmd1)
    
    system(cmd1)
    
    # move outputs to correct folder
    #cmd2 <- paste0( "mv ", feature_j, " ", prefix, "/" )
    #print(cmd2)
    #system(cmd2)
    
}

# read in phenotype metadata
meta <- readr::read_tsv(meta_file, col_names = c("chr", "start", "end", "feature") )

# keep only the specified chromosome
stopifnot(chrom %in% meta$chr)
meta_loc <- meta[ meta$chr == chrom,] 

head(meta_loc)

meta_loc <- split_chunks(meta_loc,i_chunk , n_chunk)

message( " * ", nrow(meta_loc), " features in chunk ", i_chunk, " of ", n_chunk )

#stopifnot( nrow(meta_loc) > 0 )

# write out chunked metadata and use in mmQTL script as metadata - may save time
meta_chunk_file <- paste0(prefix, "/", chrom, "_chunk_", i_chunk, "_meta.tsv" )
readr::write_tsv(meta_loc, meta_chunk_file, col_names = FALSE)


out_file <- paste0(prefix, "/", chrom, "_chunk_", i_chunk, "_output.txt" )

if( nrow(meta_loc) > 0){
# for the chunk of features, iterate through 
    for( j in 1:nrow(meta_loc) ){

        run_mmQTL(meta_loc, j)

    }
}

writeLines( c("success"), out_file)

