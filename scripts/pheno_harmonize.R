## harmonise features across datasets

# arguments

# prefix - where to write out harmonised matrices

# positional arguments - paths to each matrix
library(optparse)

option_list <- list(
   make_option(c('--prefix'), help = 'stem of out file', default = "results/example/example"),
   make_option(c('--metadata'), help = 'phenotype metadata file', default = ""),
   make_option(c('--mode'),help = 'how to harmonise - normal (genes, transcripts), SUPPA, leafcutter, txrevise', default="normal"),
   make_option(c('--min_datasets'), help = ' the minimum number of datasets to keep a feature', default = 3 )
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser, positional_arguments = TRUE)

print(opt)

prefix <- opt$options$prefix

metadata <- opt$options$metadata
inputs <- opt$args
mode <- opt$options$mode
min_datasets <- opt$options$min_datasets

is.null(metadata)
print(metadata)
print(inputs)

library(tidyverse)
message(" * ", length(inputs), " datasets loaded" )

print(inputs)

res <- map( inputs, read_tsv)

names(res) <- paste0( prefix,  gsub(".tsv.gz", ".harmonised.tsv", basename(inputs) ) )

print(names(res) )
#save.image("harmonise_debug.RData")

# for leafcutter, throw out cluster ID numbers
# NOW keep those in
#if( mode == "leafcutter" ){
#    res <- map( res ,~{ .x$feature <- gsub(":clu_[0-9]+_[+-]", "", .x$feature); .x })
#}

# remove meddlesome ";" from feature IDs in SUPPA
if( mode == "SUPPA" ){
    res <- map( res ,~{ .x$feature <- gsub(";", ":", .x$feature); .x })
}

# combine features, keep features present in at least 3 datasets - the minimum for mmQTL anyway
all_features <- unlist(map(res, "feature" ))  
feature_tally <- enframe(table(all_features) )
feature_tally <- feature_tally[ feature_tally$value >= min_datasets,]
 
shared_features <- feature_tally$name

if( mode == "normal"){

 meta <- read_tsv(metadata)
    meta <- meta[ meta$feature %in% shared_features, ]
    shared_features <- meta$feature
}

if( mode == "leafcutter" ){
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


if( mode == "SUPPA" ){
   feature_split <- stringr::str_split_fixed(shared_features, ":|-", 8)
    meta <- data.frame(
        chr = feature_split[,3],
        start = as.numeric(feature_split[,5]),
        end = as.numeric(feature_split[,6]),
        feature = shared_features
    ) 


}

message(" * ", length(shared_features) , " shared features between the datasets" )


res <- map(res, ~{ 
    row.names(.x) <- .x$feature
    df <- data.frame( Gene = shared_features, stringsAsFactors = FALSE)
    df <- left_join(df, .x, by = c("Gene" =  names(.x)[1] ) )
    df <- tidyr::drop_na(df)
    #df[ is.na(df) ] <- 0
    #.x <- .x[ shared_features, ]
    #names(.x)[1] <- "Gene"
    #.x <- tidyr::drop_na(.x) # remove NA rows from missing features
    #return(.x)
    return(df)
    })

walk2( res, names(res), ~{ write_tsv(.x, .y) })

write_tsv(meta, paste0(prefix, "phenotype_metadata.tsv"), col_names = FALSE)

