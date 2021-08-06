log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(NOISeq)
library(readr)
library(dplyr)

STEP1 = snakemake@params[["step1"]]
STEP2 = snakemake@params[["step2"]]
STEP3 = snakemake@params[["step3"]]
##########################################
## ARSyN to reduce batch effect
##########################################
{
  cat("Loading data\n")
  load(snakemake@input[[1]])
  
  mydata <- NOISeq::readData(
    data = full$M, 
    factors = full$target %>% select(group))
  
  cat("Performing ARSyN for batch correction")
  myARSyN <- ARSyNseq(mydata, norm = "n", logtransf = FALSE)
  
  ##Saving everything
  full_arsyn <- list(M = assayData(myARSyN)$exprs, annot = full$annot, 
                                targets = full$targets)
  
  stopifnot(nrow(full$M) == nrow(full_arsyn$annot))
  stopifnot(all(rownames(full_arsyn$M) == rownames(full_arsyn$annot)))
  
  full <- full_arsyn
  cat("Saving", paste(STEP1, STEP2, STEP3,  "norm_arsyn.RData", sep = "_"), "\n")
  save(full, file=snakemake@output[["arsyn_rdata"]], compress="xz")
  
  cat("Generating data matrices with arsyn for Aracne\n")
  ## Data matrices for Aracne
  #normal samples
  normal <- as_tibble(full$M[,full$targets$group == "normal"])
  normal <- bind_cols(gene=as.character(full$annot$gene_id), normal)
  
  #cancer samples
  cancer <- as_tibble(full$M[,full$targets$group == "cancer"])
  cancer <- bind_cols(gene=as.character(full$annot$gene_id), cancer)
  
  #gene_ids
  symbols <- full$annot %>% select(gene_id)
  
  cat("Saving arsyn data\n")
  
  cat("Saving", paste(STEP1, STEP2, STEP3,  "norm_arsyn_normal.tsv", sep = "_"), "\n")
  write_tsv(normal, file=snakemake@output[["normal_matrix"]])
  
  cat("Saving", paste(STEP1, STEP2, STEP3, "norm_arsyn_cancer.tsv", sep = "_"), "\n")
  write_tsv(cancer, file=snakemake@output[["cancer_matrix"]])
  
  cat("Saving", paste(STEP1, STEP2, STEP3,  "norm_arsyn_genelist.txt", sep = "_"), "\n")
  write_tsv(symbols, file=snakemake@output[["gene_list"]], col_names = F)
}
