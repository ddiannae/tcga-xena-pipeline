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
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(dplyr)

##########################################################################
load(snakemake@input[[1]])
{ 
  ### We keep only genes with mean expression count > 10 
  exp_genes <- apply(full$M, 1, function(x) mean(x)>10)
  egtable <- table(exp_genes)
  cat("There are", egtable[[1]], "genes with mean expression count < 10", egtable[[2]], "with mean count > 10 \n")
  exp_genes <- names(exp_genes[exp_genes == TRUE])

  ##Filtering low expression genes
  mean10 <- list(M = full$M[exp_genes, full$targets$id],
                 annot = full$annot %>% filter(gene_id %in% exp_genes),
                 targets = full$targets)
                 
  mean10$annot <- mean10$annot %>% arrange(gene_id)
  mean10$M <- mean10$M[order(row.names(mean10$M)), mean10$targets$id]
  
  stopifnot(all(mean10$annot$gene_id == rownames(mean10$M)))
  cat(paste0("Genes in matrix and annotation match positions \n"))
  
  stopifnot(all(mean10$targets$id == colnames(mean10$M)))
  cat(paste0("Samples in matrix and annotation match positions \n"))
  
  row.names(mean10$annot) <- mean10$annot$gene_id
  
  cat("Saving mean10_protein_coding.RData \n") 
  save(mean10, file=snakemake@output[[1]], compress="xz")
}
###########################################################################
