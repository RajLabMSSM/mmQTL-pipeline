library(tidyverse)
library(optparse)

option_list <- list(make_option(c('--PEER_cov'), help = '', default = ''),
                    make_option(c('--known_cov'), help = '', default = ''),
                    make_option(c('--prefix'), help = '', default = ''))

option.parser <- OptionParser(option_list = option_list)
opt <- parse_args(option.parser)

peer_covariate_file <- opt$PEER_cov
known_covariate_file <- opt$known_cov
prefix <- opt$prefix
outFile <- paste0(opt$prefix, "_covariates.txt")

get_all_covs <- function(peer,known_cov){
  df_peer <- read_tsv(peer)
  nrow(df_peer)
  df_cov <- read_tsv(known_cov)
  
  if(nrow(df_peer) > 0){
    message( " * merging covars with PEER covariates" )
    covariates <- df_cov %>% select(all_of(names(df_peer)))
    merged_covs <- rbind(df_peer,covariates)
  } else {
    message(" * only covars")
    merged_covs <- df_cov
  }
  
  return(merged_covs)
}

dd <- get_all_covs(peer_covariate_file,known_covariate_file)

write_tsv(dd, file = outFile)
