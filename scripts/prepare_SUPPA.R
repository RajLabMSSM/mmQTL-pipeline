## Prepare SUPPA event PSIs for mmQTL

library(optparse)
library(readr)
## read in arguments
# sk - sample key 
# tx - transcript TPM matrix
# prefix - the output prefix
# ioe - the SUPPA event IOE file


option_list <- list(
    make_option(c('--key'), help = 'sample key', default = ""),
    make_option(c('--events'), help = 'SUPPA events in IOE format', default = ""),
    make_option(c('--pheno'), help = 'the transcript matrix', default = ""),
    make_option(c('--prefix'), help = 'stem of out file', default = "results/example/example")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

sk <- opt$key
tx <- opt$pheno
events <- opt$events
prefix <- opt$prefix



tx_out <- paste0( prefix, "_transcripts_formatted.tsv")

psi_in <- paste0( prefix, ".psi")
pheno_out <- paste0( prefix, ".SUPPA.phenotype_matrix.tsv.gz")
meta_out <- paste0( prefix, ".SUPPA.phenotype_meta.tsv.gz")
# read in tx matrix
# remove samples not present in sample key
# rename samples with the participant ID
# write out in SUPPA format - transcript IDs as rownames, columns start with participant IDs

tx_matrix <- as.data.frame(readr::read_tsv(tx, col_names = TRUE))
samples <- readr::read_tsv(sk, col_names = TRUE)

row.names(tx_matrix) <- tx_matrix$transcript_id
tx_matrix$transcript_id <- NULL

tx_matrix <- tx_matrix[ ,samples$sample_id] 
names(tx_matrix) <- samples$participant_id

write.table(tx_matrix, file = tx_out,  sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE) 

save.image("debug.RData")

message(" * ", nrow(tx_matrix), " transcripts for testing")

message(" * running SUPPA")
## run SUPPA command 
cmd <- paste0('ml anaconda3; CONDA_BASE=$(conda info --base); source $CONDA_BASE/etc/profile.d/conda.sh; ml purge; conda activate isoseq-pipeline; suppa.py psiPerEvent -i ', events, ' -e ', tx_out, ' -o ', prefix )
#"suppa.py psiPerEvent -i {input.events} -e {input.tpm} -o {params.prefix};"

system(cmd)

## read in PSI matrix
# add feature  column, drop Nan rows
# how to normalise PSI?
# scale and centre

pheno <- read.table(psi_in, sep = "\t", header= TRUE, check.names = FALSE)

message(" * ", nrow(pheno), " events read in" )

# scale and centre to means of 0 and SD of 1 
pheno <- as.data.frame(t(scale(t(pheno) )))

# remove rows with NaN
pheno <- tidyr::drop_na(pheno)

# put in first column
pheno <- tibble::rownames_to_column(pheno, var = "feature")

feature_split <- stringr::str_split_fixed(pheno$feature, ":|-", 7)
meta <- data.frame(
    chr = feature_split[,2],
    start = as.numeric(feature_split[,4]),
    end = as.numeric(feature_split[,5]),
    feature = pheno$feature
    )

message(" * ", nrow(pheno), " events kept for QTL mapping" )

write_tsv(pheno, pheno_out)
write_tsv(meta, meta_out)

