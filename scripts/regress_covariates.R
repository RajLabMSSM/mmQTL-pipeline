# Regress covariates from phenotype file
## Jack Humphrey
## 2021
library(optparse)

option_list <- list(
    make_option(c('--pheno'), help = 'phenotype count matrix', default = ""),
    make_option(c('--cov'), help = 'covariate matrix', default = ""),
    make_option(c('--out'), help = 'path to output file', default = "")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

pheno_file <- opt$pheno
cov_file <- opt$cov
out_file <- opt$out

library(tidyverse)

pheno <- read_tsv(pheno_file)
# read_tsv was assuming InferredCov was Infinite, for some stupid reason
cov <- read.table(cov_file, header=TRUE)

save.image("debug.RData")

# transpose covariate matrix
cov_t <- 
    pivot_longer(cov, names_to= "sample", values_to = "value", !ID ) %>% 
    pivot_wider(names_from = "ID", values_from = "value") %>% 
    column_to_rownames(var = "sample")

pheno_r <- column_to_rownames(pheno, "feature") 

## regress covariates from phenotypes
pheno_norm <- 
  limma::removeBatchEffect(pheno_r, covariates = cov_t )

#save.image("debug.RData")

# transpose back
pheno_df <- pheno_norm %>% as.data.frame() %>% rownames_to_column( var = "feature")

write_tsv(pheno_df, out_file )


