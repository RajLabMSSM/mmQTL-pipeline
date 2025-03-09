
options(echo=TRUE)
## RUN MMQTL

# Kailash BP & Jack Humphrey
# 2025
# wraps around the MMQTL binary 
# for massive parallelisation

## Arguments:
# all the paths to the input file lists needed
# path to phenotype metadata
# chunk number (n) and chunk number (i of n)

library(parallel)
library(optparse)
option_list <- list(
    make_option(c('--chrom'), help = "which chromosome to run on", default = ""),
    make_option(c('--threads'), help = "number of threads for parallelization", default = 4),
    make_option(c('--pheno_meta'), help = 'phenotype metadata', default = ""),
    make_option(c('--chunk_total', '-n'), help = 'total number of chunks', default = 10),
    make_option(c('--chunk', '-i'), help = 'current chunk number', default = 1),
    make_option(c('--prefix'), help = 'stem of out file', default = "results/example/example"),
    make_option(c('--geno_file'), help = "path to geno list file" ),
    make_option(c('--pheno_file'), help = "path to pheno list file"),
    make_option(c('--grm_file'), help = "path to GRM list file"),
    make_option(c('--cov_file'), help = "path to covariate_file"),
    make_option(c('--eQTL_number'), help = "0 for primary QTL, 1 for primary + secondary QTL, 3 for ... you guessed it right :)",  default = 0),
    make_option(c('--mmQTL'), help = "full path to MMQTL executable", default = "MMQTL26a")
)

absPath <- function(path){
    if(is.null(path) ){ return()}
    if(file.exists(path) ){
        normalizePath(path)
    }
}
option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

num_cores <- as.numeric(opt$threads)  # Capture the number of threads

# Print the number of cores being used
message("Using ", num_cores, " cores for parallel processing")

chrom <- opt$chrom
meta_file <- absPath(opt$pheno_meta)
geno_file <- absPath(opt$geno_file)
pheno_file <- absPath(opt$pheno_file)
cov_file <- absPath(opt$cov_file)
grm_file <- absPath(opt$grm_file)
eQTL_number <- opt$eQTL_number
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

# See https://github.com/jxzb1988/MMQTL/blob/master/src_1.5.0/MeQTLPolyG.cpp for options

run_mmQTL <- function(meta_loc, j){
     
    feature_j <- meta_loc[j, ]$feature

    cmd1 <- paste0( mmqtl_bin,
        " -b ", 
        " -P ", pheno_file,
        " -Z ", geno_file,
        " -R ", grm_file,
        " -a ", meta_chunk_file,
        " -A random ",
        " -gene ", feature_j,
        " --eQTL_number ", eQTL_number,
        " -V ", format(cis_window, scientific = FALSE, digits = 1),
        " --Han" # Apply Han & Eskin method to adjust results of random-effect model
    )
    
    tryCatch({
      print(cmd1)  # Print the command to the screen before executing it
      system(cmd1)
    }, error = function(e) {
      message("Error running MMQTL for feature: ", feature_j)
    })
}

# Read phenotype metadata
meta <- readr::read_tsv(meta_file, col_names = c("chr", "start", "end", "feature"))

# Define output paths
out_file <- paste0(prefix, "/", chrom, "_chunk_", i_chunk, "_output.txt")
meta_chunk_file <- paste0(prefix, "/", chrom, "_chunk_", i_chunk, "_meta.tsv")

# Check if chromosome is present in the metadata
if (chrom %in% meta$chr == FALSE) {
  warning_message <- "WARNING: CHROM not in meta file. Either meta file is wrong, or you are not testing all chromosomes."
  writeLines(c(warning_message), out_file)
  readr::write_tsv(data.frame(), meta_chunk_file, col_names = FALSE)  # Write empty metadata file
  quit(save = "no")  # Stop execution, quit script, don't save workspace
}

# Filter metadata for the specified chromosome
meta_loc <- meta[meta$chr == chrom,]

# If there are no features to process, handle it gracefully
if (nrow(meta_loc) == 0) {
  warning_message <- paste("WARNING: No features to process for chunk", i_chunk)
  writeLines(c(warning_message), out_file)
  readr::write_tsv(data.frame(), meta_chunk_file, col_names = FALSE)  # Write empty metadata file
  quit(save = "no")  # Stop execution, quit script, don't save workspace
}

# Split metadata into chunks and write the chunk to a file
meta_loc <- split_chunks(meta_loc, i_chunk, n_chunk)
readr::write_tsv(meta_loc, meta_chunk_file, col_names = FALSE)

# Log the number of features in the chunk
message(" * ", nrow(meta_loc), " features in chunk ", i_chunk, " of ", n_chunk)

# If there are features, process them
if (nrow(meta_loc) > 0) {
  # Run MMQTL in parallel using `mclapply`
  mclapply(1:nrow(meta_loc), function(j) {
    run_mmQTL(meta_loc, j)}, mc.cores = num_cores)  # Use the allocated number of cores
}

# Write success message to the output file
writeLines(c("success"), out_file)
