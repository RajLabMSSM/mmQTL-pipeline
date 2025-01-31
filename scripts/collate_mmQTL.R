## Collate MMQTL results

## for each feature in a chromosome
## read in result file
## match on SNP coordinates and sort 
## convert Random_Z to P-value
## Bonferroni adjust
## write top association to a file
## write sorted associations to a file - no header
## concatenate all the pieces into a single chromosome file

library(optparse)

option_list <- list(
   make_option(c('--prefix'), help = 'stem of out file', default = "results/example/example"),
   make_option(c('--metadata'), help = 'phenotype metadata file', default = ""),
   make_option(c('--chrom'), help = 'the chromosome', default = ""),
   make_option(c('--geno'), help = 'path to genotype folder', default = "")
   # add data key here to get order of datasets
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

print(opt)

prefix <- opt$prefix
chrom <- opt$chrom
pheno_meta <- opt$metadata
top_file <- paste0(prefix, chrom, "_top_assoc.tsv" )

library(readr)
library(purrr)
library(dplyr)

message(" * collating ", chrom )

meta <- read_tsv(pheno_meta, col_names = FALSE)
names(meta)[1:4] <- c("chr", "start", "end", "feature") 

if (chrom %in% meta$chr == FALSE) {
  message("WARNING: Chromosome ", chrom, " not found in metadata. Exiting.")
  write_tsv(data.frame(), top_file)  # Write empty file and exit
  quit(save = "no")
}

meta_loc <- meta[ meta$chr %in% chrom, ]

all_files <- list.files(prefix, pattern = "test_peak_1_statistical_signal\\.gz$", recursive = TRUE, full.names = FALSE)

files_loc <- all_files[ dirname(all_files) %in% meta_loc$feature ] 

features_loc <- dirname(files_loc)
files_loc <- paste0( prefix, files_loc)

names(files_loc) <- features_loc

if (length(files_loc) == 0) {
  message("WARNING: No files found to collate for chromosome ", chrom)
  write_tsv(data.frame(), top_file)  # Write empty file and exit
  quit(save = "no")
}

message( " * ", length(files_loc), " to collate" )

geno_folder <- opt$geno
## get SNP coordinates from the coordinate BIM files for each dataset
bim_files <- list.files( geno_folder, pattern = paste0("*", chrom, ".bim" ), recursive = TRUE, full.names = TRUE )
## read in and get distinct rows
bim_res <- purrr::map_df(bim_files, read_tsv, col_names = c("chr", "Variant","n", "pos", "ref", "alt") ) %>% distinct()

bim_res$chr <- paste0("chr", bim_res$chr)
bim_res$n <- NULL

# get number of datasets
n_datasets <- length(bim_files)

# some features do not use all datasets in their meta-analysis
# and thus have missing columns which messes up concatenation
## prepare columns according to the number of datasets
data_cols <- paste0(rep(c("beta_", "sd_", "z_"),3), paste0("tissue_", sort(rep(0:(n_datasets-1),3) ) ) )

all_cols <- c("Variant", "Allele", data_cols, "fixed_beta", "fixed_sd", "fixed_z", "Random_Z" )

dummy_cols <- setNames(data.frame(matrix(ncol = length(all_cols), nrow = 1) ), all_cols) 

## functions
add_P <- function(data){
    data$Fixed_P <- 2*pnorm(q=abs(data$fixed_z), lower.tail=FALSE)
    data$Random_P <- 2*pnorm(q=abs(data$Random_Z), lower.tail=FALSE)
    data$Fixed_bonf <- p.adjust(data$Fixed_P, method = "bonferroni")
    data$Random_bonf <- p.adjust(data$Random_P, method = "bonferroni")
    data$Fixed_FDR <- p.adjust(data$Fixed_P, method = "BH" )
    data$Random_FDR <- p.adjust(data$Random_P, method = "BH" )
    return(data)
}

top_assoc <- list()

for(feature in features_loc){
    
    message(" * Processing feature: ", feature)
    
    d <- read_tsv(files_loc[ feature ], col_types = list(Allele = "c"))
    
    # ignore empty or malformed files
    if( nrow(d) == 0 ){
      next
    }
    
    if( !"fixed_z" %in% names(d) ){ next}
    # standardise columns
    d <- bind_rows(dummy_cols, d)
    d <- d[2:nrow(d), ]

    d$feature <- feature
    # add P values
    d <- add_P(d)
    # add in SNP info
    # standardise columns HERE
    d <- left_join(d, bim_res, by = "Variant") %>%
        arrange(pos) %>%
        select(feature, variant_id = Variant, chr, pos, ref, alt, everything() )
    # write out
    out_file <- paste0(prefix, chrom, "_", feature, "_all_nominal.tsv", ".gz")
    write_tsv(d, out_file, col_names = FALSE) 
    top <- arrange(d, Random_P ) %>% head(1)
    top_assoc[[feature]] <- top
}

# Check if any top associations were found
if (length(top_assoc) == 0) {
  message("No top associations found. Writing an empty top file.")
  write_tsv(data.frame(), top_file)
} else {
  # Combine top associations and write to file
  top_res <- bind_rows(top_assoc)
  write_tsv(top_res, top_file)
}

message(" * Collation complete for chromosome: ", chrom)