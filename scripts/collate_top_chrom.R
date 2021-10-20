## collate top associations

# args
# --out - the path to output file
# positional arguments for input files
library(optparse)

option_list <- list(
   make_option(c('--output', '-o'), help = 'stem of out file', default = "results/example/example")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser, positional_arguments = TRUE)


output <- opt$options$output

inputs <- opt$args

library(tidyverse)
res <- map_df( inputs, read_tsv)
res$qval <- qvalue::qvalue(res$Random_FDR)$qval
res <- arrange(res, qval)
write_tsv(res, output)
