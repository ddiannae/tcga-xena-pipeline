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

  non_zero_genes <- apply(full$M, 1, function(x) sum(x == 0) < ncol(full$M)/2)
  egtable <- table(non_zero_genes)
  cat("There are", egtable[[1]], "genes with 0 values in more than 50% samples", egtable[[2]], 
      " genes with 0 values in less than 50% samples \n")
  non_zero_genes <- names(non_zero_genes[non_zero_genes == TRUE])
  
  genes_pass <- intersect(exp_genes, non_zero_genes)
  
  ##Filtering low expression genes
  full <- list(M = round(full$M[genes_pass, full$targets$id]),
                 annot = full$annot %>% filter(gene_id %in% genes_pass),
                 targets = full$targets)
                 
  full$annot <- full$annot %>% arrange(gene_id)
  full$M <- full$M[order(row.names(full$M)), full$targets$id]
  
  stopifnot(all(full$annot$gene_id == rownames(full$M)))
  cat(paste0("Genes in matrix and annotation match positions \n"))
  
  stopifnot(all(full$targets$id == colnames(full$M)))
  cat(paste0("Samples in matrix and annotation match positions \n"))
  
  row.names(full$annot) <- full$annot$gene_id
  
  cat("Saving full_filtered.RData \n") 
  save(full, file=snakemake@output[[1]], compress="xz")
}
###########################################################################
