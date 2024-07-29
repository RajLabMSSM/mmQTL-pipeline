library(optparse)
library(tidyverse)

option_list <- list(make_option(c('--inputDir'), help = '', default = './'))

option.parser <- OptionParser(option_list = option_list)
opt <- parse_args(option.parser)

inDir <- opt$inputDir

files <- list.files(path = inDir, pattern = "*all_sites_pileup_annotation.tsv.gz", full.names = TRUE, recursive = TRUE)
data <- map(files, read_tsv)

all_sites <- map(data, "ESid2") %>% unlist()
all_sites <- unique(all_sites)
all_sites <- data.frame(all_sites)
all_meta <- all_sites %>%
  mutate(chr = str_split_fixed(all_sites, ":", 5)[,1],
         chrom = gsub("chr","",chr)) %>%
  filter(chrom %in% c(1:22)) %>%
  mutate(CHROM = as.numeric(chrom),
         end = str_split_fixed(all_sites, ":", 5)[,2],
         start = as.numeric(end)-1) %>%
  arrange(., CHROM, start) %>%
  rename("feature" = "all_sites") %>%
  filter(grepl("A:G|T:C",feature)) %>%
  select(chr, start, end, feature)

annotations <- function(x){
  y <- x %>% 
    rename("feature" = "ESid2") %>%
    filter(feature %in% all_meta$feature) %>%
    select(feature, Func.refGene, Gene.refGene, ExonicFunc.refGene, AAChange.refGene, rmsk, REDIportal_info, ensembl_id)
}

all_meta2 <- lapply(data,annotations)
all_meta2 <- do.call("rbind",all_meta2)
all_meta2 <- all_meta2 %>% distinct(feature, .keep_all = TRUE)

meta_file <- paste0(inDir, "pheno_meta.tsv")
meta_file2 <- paste0(inDir, "annotations_meta.tsv")
write.table(all_meta, file = meta_file, row.names = FALSE, quote = FALSE, sep = "\t")
write_tsv(all_meta2, file = meta_file2)
