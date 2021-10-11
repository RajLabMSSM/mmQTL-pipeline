# library(leafcutter, quietly=TRUE)
suppressMessages(library(argparser, quietly=TRUE))
suppressMessages(library(dplyr, quietly=TRUE))
suppressMessages(library(foreach, quietly=TRUE))

# leafcutter functions:

# now takes strand into account - match on strand!


#' Make a data.frame of meta data about the introns
#' @param introns Names of the introns
#' @return Data.frame with chr, start, end, cluster id
#' @export
get_intron_meta <- function(introns) {
  intron_meta <- do.call(rbind, strsplit(introns,":"))
  colnames(intron_meta) <- c("chr","start","end","clu")
  intron_meta <- as.data.frame(intron_meta, stringsAsFactors=FALSE)
  intron_meta$start <- as.numeric(intron_meta$start)
  intron_meta$end <- as.numeric(intron_meta$end)
  # add in strand info
  intron_meta$strand <- do.call(rbind, strsplit(intron_meta$clu, "_"))[,3]
  intron_meta
}

#' Work out which gene each cluster belongs to. Note the chromosome names used in the two inputs must match.
#' @param intron_meta Data frame describing the introns, usually from get_intron_meta
#' @param exons_table Table of exons, see e.g. /data/gencode19_exons.txt.gz
#' @param flip Whether to flip strand - for when junctions are reverse to reference
#' @return Data.frame with cluster ids and genes separated by commas
#' @import dplyr
#' @export
map_clusters_to_genes <- function(intron_meta, exons_table, flip = FALSE) {
  gene_df <- foreach (chr=sort(unique(exons_table$chr)), .combine=rbind) %dopar% {
    
    intron_chr <- intron_meta[ intron_meta$chr==chr, ]
    exons_chr <- exons_table[exons_table$chr==chr, ]
    
    if(flip == TRUE){
        intron_chr$strand <- ifelse(intron_chr$strand == "+", "-", "+")
    }
    print(head(exons_chr))
    print(head(intron_chr))
    
    exons_chr$temp <- exons_chr$start
    intron_chr$temp <- intron_chr$end
    three_prime_matches <- inner_join( intron_chr, exons_chr, by=c("strand", "temp") )

    exons_chr$temp <- exons_chr$end
    intron_chr$temp <- intron_chr$start
    five_prime_matches <- inner_join( intron_chr, exons_chr, by=c( "strand", "temp") )

    all_matches <- rbind(three_prime_matches, five_prime_matches)[ , c("clu", "gene_name")]

    all_matches <- all_matches[!duplicated(all_matches),]

    if (nrow(all_matches)==0) return(data.frame(clu = NA, gene_name = NA))
    all_matches$clu <- paste(chr,all_matches$clu,sep=':')
    all_matches
  }
  clu_df <- gene_df %>% group_by(clu) %>% summarize(genes=paste(gene_name, collapse = ","))
  class(clu_df) <- "data.frame"
  clu_df
}


p <- arg_parser("LeafCutter: map clusters to genes")
p <- add_argument(p, "intron_counts_file", help="Intron counts file from LeafCutter, typically <prefix>_perind.counts.gz")
p <- add_argument(p, "exon_file", help="File listing all unique exons in annotation. Must have columns: chr, start, end, strand, gene_id[, gene_name].")
p <- add_argument(p, "output_name", help="Output file name")
p <- add_argument(p, "--strand", short="-o", help="whether junctions are stranded or not - default is 0 - unstranded. 1 = forward stranded; 2 = reverse stranded", default="0")
p <- add_argument(p, "--output_dir", short="-o", help="Output directory", default=".")
argv <- parse_args(p)

cat("LeafCutter: mapping clusters to genes\n")
intron_counts <- read.table(argv$intron_counts_file, header=TRUE, check.names=FALSE, row.names=1)
intron_meta <- get_intron_meta(rownames(intron_counts))

exon_table <- read.table(argv$exon_file, header=TRUE, stringsAsFactors=FALSE)
# JH - my table only has gene_name
#stopifnot(is.element('gene_name', colnames(exon_table)))
stopifnot(is.element('gene_id', colnames(exon_table)))
exon_table[, 'gene_name'] <- exon_table[, 'gene_id']

m <- map_clusters_to_genes(intron_meta, exon_table)
message(paste0(" * ", nrow(m), " gene-cluster matches!"))
m2 <- map_clusters_to_genes(intron_meta, exon_table, flip = TRUE) 
message(paste0(" * with strand flipped: ", nrow(m2), " gene-cluster matches!"))

if(nrow(m) > nrow(m2) ){ 
    out_file <- m
}else{
    out_file <- m2
}
write.table(out_file, file.path(argv$output_dir, argv$output_name), sep = "\t", quote=FALSE, row.names=FALSE)
