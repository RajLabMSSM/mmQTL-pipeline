## Make a phenotype matrix file from a GTF
## Jack Humphrey
## 2021
# arguments:
## GTF
## prefix - for output
## mode - gene or transcript or other - for SUPPA etc
in_gtf <- "/sc/arion/projects/ad-omics/data/references//hg38_reference/GENCODE/gencode.v38.primary_assembly/gencode.v38.primary_assembly.annotation.gtf"
prefix <- "test"
mode <- "gene"

shhh <- suppressPackageStartupMessages

library(optparse)

option_list <- list(
    make_option(c('--prefix'), help = 'stem of out file', default = "test"),
    make_option(c('--mode'), help = 'whether to get out gene or transcript', default = "gene"),
    make_option(c('--gtf'), help = 'the input GTF file', default = "")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

prefix <- opt$prefix
in_gtf <- opt$gtf
mode <- opt$mode

stopifnot( file.exists(in_gtf) )
stopifnot( mode %in% c("gene", "transcript" ) )


out_file <- paste0(prefix, "_", mode, "_pheno_matrix.tsv")

shhh(library(rtracklayer))

message( " * mode selected is ", mode )
message( " * reading in GTF file..." )

gtf <- import(in_gtf)

# create phenotype metadata from the GTF used
# if gene mode then output info for each gene
# if transcript then output info for each transcript
# TODO - how to parse SUPPA and/or txrevise?
if( mode == "gene" ){
   gtf_loc <- gtf[ gtf$type == "gene" ]
   features = gtf_loc$gene_id
}

if( mode == "transcript" ){
    gtf_loc <- gtf[ gtf$type == "transcript" ]
    features <- gtf_loc$transcript_id
}

out <- data.frame( 
    chr = seqnames(gtf_loc),
    start = start(gtf_loc),
    end = end(gtf_loc),
    feature = features,
    stringsAsFactors = FALSE
)

message( " * ", nrow(out), " features kept for metadata" )

if( mode == "transcript" ){
    out$group <- gtf_loc$gene_id
}

message( " * writing to ", out_file)

readr::write_tsv(out, out_file)

