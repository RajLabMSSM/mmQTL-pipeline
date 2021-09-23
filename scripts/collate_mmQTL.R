## Collate MMQTL results

## for each feature in a chromosome
## read in result file
## match on SNP coordinates and sort 
## convert Random_Z to P-value
## Bonferroni adjust
## write top association to a file
## write sorted associations to a file - no header
## concatenate all the pieces into a single chromosome file

prefix <- "results/example/mmQTL/output/"
chrom <- "chr1"
pheno_meta <- "example/phenotype_metadata_72genes.tsv"

meta <- read_tsv(pheno_meta, col_names = c("chr", "start", "end", "feature") )

meta_loc <- meta[ meta$chr %in% chrom, ]

all_files <- list.files(prefix, pattern = "test_peak_1_statistical_signal", recursive = TRUE)

files_loc <- all_files[ dirname(all_files) %in% meta_loc$feature ] 

features_loc <- dirname(files_loc)
files_loc <- paste0( prefix, files_loc)

names(files_loc) <- features_loc

## get SNP coordinates from the coordinate BIM files for each dataset
bim_files <- list.files( "results/example/", pattern = paste0("*", chrom, ".bim" ), recursive = TRUE )
## read in and get distinct rows
bim_res <- map_df(bim_files, read_tsv, col_names = c("chr", "Variant","n", "pos", "ref", "alt") ) %>% distinct()

bim_res$chr <- paste0("chr", bim_res$chr)
bim_res$n <- NULL

## functions
add_P <- function(data){
    data$Fixed_P <- 2*pnorm(q=abs(data$fixed_z), lower.tail=FALSE)
    data$Random_P <- 2*pnorm(q=abs(data$Random_Z), lower.tail=FALSE)
    data$Fixed_q <- p.adjust(data$Fixed_P, method = "bonferroni")
    data$Random_q <- p.adjust(data$Random_P, method = "bonferroni")
    return(data)
}

top_assoc <- list()
for(feature in features_loc){
    print(feature)
    f <- files_loc[ feature ]
    d <- read_tsv(f)
    # add P values
    d <- add_P(d)
    d <- left_join(d, bim_res, by = "Variant") %>%
        arrange(pos) %>%
        select(chr, pos, Variant, ref, alt, everything() )
    # write out
    out_file <- paste0(prefix, chrom, "_", feature, "_all_nominal.tsv" )
    write_tsv(d, out_file, col_names = FALSE) 
    top <- arrange(d, Random_q ) %>% head(1)
    top_assoc[[feature]] <- top
    
}

top_res <- bind_rows(top_assoc)

top_file <- paste0(prefix, chrom, "_top_assoc.tsv" )
write_tsv(top_res, top_file)
