## Prepare Leafcutter matrix for meta-analysis
## this way all the samples are clustered together first
## then split off into separate datasets later

## Jack Humphrey
## 2021

## arguments:
## junction list
## exons
## prefix 

library(optparse)
library(readr)

option_list <- list(
   make_option(c('--prefix'), help = 'stem of out file, including output directory', default = "results/example/example"),
   make_option(c('--junctions'), help = 'list of junction files', default = ""),
   make_option(c('--exons'), help = 'the exon metadata', default = "")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

#print(opt)

prefix <- opt$prefix
exons <- opt$exons
junctions  <- opt$junctions

outFile <- basename(prefix)
outFolder <- dirname(prefix)

if( !dir.exists(outFolder) ){ dir.create(outFolder) }

junc_files <- readLines(junctions)
stopifnot( all(file.exists(junc_files) ) )

## step 1: Leafcutter clustering

cmd <- paste0("python scripts/leafcutter/leafcutter_cluster_regtools.py -j ", junctions, " -o ", outFile, " -r ", outFolder, " -l 250000 --checkchrom=True" )
print(" * clustering junctions ")
print(cmd)
system(cmd)

count_matrix <- paste0( prefix, "_perind.counts.gz" )
junc_matrix <- paste0( prefix, "_perind_numers.counts.gz")
stopifnot( file.exists( count_matrix) )
stopifnot( file.exists( junc_matrix ) )

gene_map <- paste0(prefix, "_gene_map.tsv")

## step 2: map gene names to clusters
#Rscript scripts/leafcutter/map_clusters_to_genes.R --output_dir junctions/ junctions/test_perind.counts.gz input/isoseq_exon_pheno_meta.tsv test
cmd <- paste("ml R/4.0.3; Rscript  scripts/leafcutter/map_clusters_to_genes.R --output_dir . ", count_matrix, " ", exons, " ", gene_map )
print(cmd)
print(" * mapping genes to clusters" )
system(cmd)

stopifnot( file.exists( gene_map) )

#save.image("debug.RData")
## step 3: prepare junction count matrix and metadata
message(" * preparing matrix and metadata")

junc_df <- read.table( junc_matrix, header=TRUE, check.names = FALSE)
gene_df <- read.table( gene_map, header=TRUE, check.names = FALSE)
gene_df$clu <- gsub(".*:", "", gene_df$clu)

get_intron_meta <- function(introns){
    intron_meta = do.call(rbind, strsplit(introns, ":"))
    colnames(intron_meta) = c("chr", "start", "end", "clu")
    intron_meta = as.data.frame(intron_meta, stringsAsFactors = F)
    intron_meta$start = as.numeric(intron_meta$start)
    intron_meta$end = as.numeric(intron_meta$end)
    intron_meta$middle = 0.5 * (intron_meta$start + intron_meta$end)
    intron_meta
}

# build metadata
intron_meta <- get_intron_meta(row.names(junc_df) )
intron_meta$feature <- row.names(junc_df)
intron_meta$gene <- gene_df$genes[ match(intron_meta$clu, gene_df$clu)]
intron_meta$gene <- gsub(",", "\\+", intron_meta$gene)

intron_meta$group <- intron_meta$clu
intron_meta <- intron_meta[, c("chr", "start", "end", "feature", "group", "gene") ]

junc_df$feature <- intron_meta$feature
row.names(junc_df) <- NULL

junc_df <- junc_df[, c("feature", names(junc_df)[1:ncol(junc_df)-1] ) ]

## write out
write_tsv(junc_df, paste0(prefix, "_leafcutter_junction_matrix.tsv.gz") )
write_tsv(intron_meta, paste0(prefix, "_leafcutter_junction_metadata.tsv.gz") )

# take out trash
trash <- paste( "rm", paste0(prefix, c("_sortedlibs", "_pooled", "_refined"), collapse = " " ), paste0(outFolder, "/*sorted.gz" ) ) 
message(" * removing temp files")
system(trash)


