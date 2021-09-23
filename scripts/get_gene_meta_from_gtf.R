## Make a phenotype matrix file from a GTF
## Jack Humphrey
## 2021

library(rtracklayer)

# arguments:
## GTF
## prefix - for output
## mode - gene or transcript or other - for SUPPA etc
in_gtf <- "/sc/arion/projects/ad-omics/data/references//hg38_reference/GENCODE/gencode.v38.primary_assembly/gencode.v38.primary_assembly.annotation.gtf"
prefix <- "test"
mode <- "gene"
out_file <- paste0(prefix, "_pheno_matrix.tsv")

gtf <- import(in_gtf)

save.image("debug.RData")

if( mode == "gene" ){
   gtf_loc <- gtf[ gtf$type == "gene" ]
   features = gtf$gene_id
}

if( mode == "transcript" ){
    gtf_loc <- gtf[ gtf$type == "transcript" ]
    features <- gtf$transcript_id
}

out <- data.frame( 
    chr = seqnames(gtf_loc),
    start = start(gtf_loc),
    end = end(gtf_loc),
    feature = features,
    stringsAsFactors = FALSE
)


if( mode == "transcript" ){
    out$group <- gtf_loc$gene_id
}

readr::write_tsv(out, out_file)

