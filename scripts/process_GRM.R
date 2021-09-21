# convert GRM 

# Roman Kosoy and Jack Humphrey
# 2021

library(optparse)

option_list <- list(
    make_option(c('--prefix'), help = 'stem of out file', default = "results/example/example")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

prefix <- opt$prefix


#WorkingFolder = "/sc/hydra/projects/roussp01b/MICROGLIA/genotyping/Combined_Phase1and2/preparing_caQTL_files/"
#GRM_BIN_file = paste0(WorkingFolder, "PLINK_for_eQTL_GRM")

ReadGRMBin=function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}

GRM <- ReadGRMBin(prefix)

## to create the matrix:

GRM_matrix <- matrix(NA, ncol = length(GRM$id$V2), nrow = length(GRM$id$V2), dimnames = list(GRM$id$V2, GRM$id$V2))
diag(GRM_matrix) <- GRM$diag

gdata::lowerTriangle(GRM_matrix, diag=FALSE, byrow=FALSE) <- GRM$off
gdata::upperTriangle(GRM_matrix, diag=FALSE, byrow=TRUE) <- GRM$off

# GRM_matrix[1:5,1:5]

out_file <- paste0(prefix, ".tsv")

write.table(GRM_matrix, out_file, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


#saveRDS(GRM_result, out_file)

#paste0(WorkingFolder, "PLINK_for_eQTL_GRM_nonBinary.RDs"))
#m(list = ls())

