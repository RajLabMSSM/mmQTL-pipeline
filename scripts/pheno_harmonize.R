## harmonise features across datasets

# arguments

# prefix - where to write out harmonised matrices

# positional arguments - paths to each matrix
library(optparse)

option_list <- list(
   make_option(c('--prefix'), help = 'stem of out file', default = "results/example/example"),
   make_option(c('--metadata'), help = 'phenotype metadata file', default = ""),
   make_option(c('--leafcutter'),action="store_true", default=FALSE)
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser, positional_arguments = TRUE)

print(opt)

prefix <- opt$options$prefix

metadata <- opt$options$metadata
inputs <- opt$args
leafcutter <- opt$options$leafcutter

is.null(metadata)
print(metadata)
print(inputs)

library(tidyverse)
message(" * ", length(inputs), " datasets loaded" )

print(inputs)

res <- map( inputs, read_tsv)

names(res) <- paste0( prefix,  gsub(".tsv.gz", ".harmonised.tsv", basename(inputs) ) )

print(names(res) )
save.image("debug.RData")

if( !leafcutter){
 shared_features <- map(res, "feature" ) %>% reduce(intersect)
   
    meta <- read_tsv(metadata)
    meta <- meta[ meta$feature %in% shared_features, ]
    shared_features <- meta$feature
}

if( leafcutter ){
    res <- map( res ,~{ .x$feature <- gsub(":clu_[0-9]+_[+-]", "", .x$feature); .x })
    shared_features <- map(res, "feature" ) %>% reduce(intersect)

    meta <- as.data.frame(stringr::str_split_fixed(shared_features, ":", 4) )
    print(head(meta) )
    names(meta) <- c("chr", "junc_start", "junc_end", "gene" )
    meta$junc_start <- as.numeric(meta$junc_start)
    meta$junc_end <- as.numeric(meta$junc_end)
    meta$end <- ceiling( meta$junc_start + ( (meta$junc_end - meta$junc_start )/2 ) )
    meta$start <- meta$end - 1
    meta$feature <- shared_features
    meta <- select(meta, chr, start, end, feature, group = gene)
}

message(" * ", length(shared_features) , " shared features between the datasets" )


res <- map(res, ~{ 
    row.names(.x) <- .x$feature
    .x <- .x[ shared_features, ]
    names(.x)[1] <- "Gene"
    return(.x)
    })

walk2( res, names(res), ~{ write_tsv(.x, .y) })

write_tsv(meta, paste0(prefix, "phenotype_metadata.tsv"), col_names = FALSE)

