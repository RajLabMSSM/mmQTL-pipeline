# Regress covariates from phenotype file
## Jack Humphrey
## 2021

## Aline modified October 2024 to handle factor as covariates specified in the code => change if you have different covariates
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

# Convert relevant columns to factors
if("age" %in% colnames(cov_t)){
  cov_t$age <- as.numeric(cov_t$age)
}
if("biological_sex" %in% colnames(cov_t)){
  cov_t$biological_sex <- as.factor(cov_t$biological_sex)
}
if("ADAR" %in% colnames(cov_t)){
  cov_t$ADAR <- as.numeric(cov_t$ADAR)
}
if("ADARB1" %in% colnames(cov_t)){
  cov_t$ADARB1 <- as.numeric(cov_t$ADARB1)
}

# remove the variable with only one factor level
cov_t <- cov_t %>%
    select_if(~ length(unique(.)) > 1)

pheno_r <- column_to_rownames(pheno, "feature")

# handle factors
design_matrix <- model.matrix(~ . - 1, data = cov_t, na.action = na.pass)

## regress covariates from phenotypes
pheno_norm <-
  limma::removeBatchEffect(pheno_r, covariates = design_matrix)

#save.image("debug.RData")

# transpose back
pheno_df <- pheno_norm %>% as.data.frame() %>% rownames_to_column( var = "feature")

write_tsv(pheno_df, out_file )
