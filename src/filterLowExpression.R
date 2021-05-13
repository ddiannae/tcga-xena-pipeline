###############################################################################
##  LOSS OF LONG RANGE CO-EXPRESSION IN CANCER
##  Analysis of RNA-Seq data.
###############################################################################
## Diana Garcia - diana.gco@gmail.com
## Date: April, 2021
## 
## Original code. Dr. Cristobal Fresno - cristobalfresno@gmail.com
## Date:  2016-12-12
###############################################################################
## Data Normalization
##  - Filter genes with median expression count under 10
###############################################################################
args <- commandArgs(trailingOnly = T)

if (length(args) < 2 ) {
  stop("Incorrect number of arguments", call.=FALSE)
} else {
  TISSUE = args[1]
  DATADIR = args[2]
}
DATADIR <- paste(DATADIR, TISSUE, sep="/")
RDATA <- paste(DATADIR, "rdata", sep="/")
##########################################################################
load(file=paste(RDATA, "raw_full.RData", sep="/"))
{ 
  ### We keep only genes with mean expression count > 10 
  exp_genes <- apply(full$M, 1, function(x) mean(x)>10)
  egtable <- table(exp_genes)
  cat("There are", egtable[[1]], "genes with mean expression count < 10", egtable[[2]], "with mean count > 10 \n")
  
  ##Filtering low expression genes
  mean10 <- list(M = full$M[exp_genes, rownames(full$targets)], annot=full$annot[exp_genes, ], 
                 targets=full$targets)
  rownames(mean10$annot) <- rownames(mean10$M)
  
  cat("Saving mean10_ProteinCoding.RData \n") 
  save(mean10, file=paste(RDATA, "mean10_proteinCoding.RData", sep="/"), compress="xz")
}
###########################################################################
