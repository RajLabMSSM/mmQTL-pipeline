## harmonise features across datasets

# arguments

# prefix - where to write out harmonised matrices

# positional arguments - paths to each matrix
library(optparse)

option_list <- list(
   make_option(c('--prefix'), help = 'stem of out file', default = "results/example/example"),
   make_option(c('--metadata'), help = 'phenotype metadata file', default = "")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser, positional_arguments = TRUE)

print(opt)

prefix <- opt$options$prefix

metadata <- opt$options$metadata
inputs <- opt$args

is.null(metadata)
print(metadata)
print(inputs)

library(purrr)
library(readr)

message(" * ", length(inputs), " datasets loaded" )

print(metadata)
meta <- read_tsv(metadata)

print(inputs)

res <- map( inputs, read_tsv)

names(res) <- paste0( prefix,  gsub(".tsv.gz", ".harmonised.tsv.gz", basename(inputs) ) )

print(names(res) )
shared_features <- map(res, "feature" ) %>% reduce(intersect)

meta <- meta[ meta$feature %in% shared_features, ]

message(" * ", length(shared_features) , " shared features between the datasets" )


res <- map(res, ~{ 
    row.names(.x) <- .x$feature
    .x <- .x[ meta$feature, ]
    return(.x)
    })

walk2( res, names(res), ~{ write_tsv(.x, .y) })

write_tsv(meta, paste0(prefix, "phenotype_metadata.tsv") )
