# Regress covariates from phenotype file
## Jack Humphrey
## 2021


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
cov <- read_tsv(cov_file)

# transpose covariate matrix
cov_t <- 
    pivot_longer(cov, names_to= "sample", values_to = "value", !ID ) %>% 
    pivot_wider(names_from = "ID", values_from = "value") %>% 
    column_to_rownames(var = "sample")

# transpose phenotype matrix
pheno_t <- 
    pivot_longer(pheno, names_to= "sample", values_to = "value", !ID ) %>%
    pivot_wider(names_from = "Gene", values_from = "value") %>% 
    column_to_rownames(var = "sample")


## regress covariates from phenotypes
norm_t <- 
  limma::removeBatchEffect(pheno_t, covariates = cov_t )

save.image("debug.RData")

# transpose back
norm <- norm %>%


write_tsv(reg, out_file )


