###############################################################################
## Data Normalization
##    - Filter genes with median expression count under 10
##    - Try Normalization combinations:
##          Within lane. For lenght: RPKM, loess, full
##          Within lane. For GC content: loess, full
##          Between lanes. TMM, full
##    - Pick the best set of normalization methods based on:
##          R2 values of cubic spline regression models for GC and lenght bias
##          RNA composition bias. Number of samples that passed the test
###############################################################################
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

##Usefull Libraries
library(ggplot2)
library(reshape2)
library(readr)
library(dplyr)
library(EDASeq)
library(NOISeq)
library(DESeq2)

STEP1 = snakemake@params[["step1"]]
STEP2 = snakemake@params[["step2"]]
STEP3 = snakemake@params[["step3"]]

RDATA <-paste(snakemake@params[["tissue_dir"]], "rdata", sep="/")
w <- 1024
h <- 1024
p <- 24
##########################################################################
{##### USER SELECTED NORMALIZATION
  cat("Loading data\n")
  load(snakemake@input[[1]])
  cat("Performing normalization. Step1: ", STEP1, ", Step2: ", STEP2, " Step3: ", STEP3, "\n")
  
  first <- strsplit(STEP1, "-")[[1]]
  second <- strsplit(STEP2, "-")[[1]]
  
  if(first[1] == "length") {
    if(first[2] != "no") {
      norm_counts <- withinLaneNormalization(full$M, full$annot$length, which = first[2])  
    } else {
      norm_counts <- full$M
    }
    if(second[1] == "gc") {
      if(second[2] != "no") {
        norm_counts <- withinLaneNormalization(norm_counts, full$annot$gc, which = second[2])    
      }
    }
  } else if (first[1] == "gc") {
    if(first[2] != "no") {
      norm_counts <- withinLaneNormalization(full$M, full$annot$gc, which = first[2])
    } else {
      norm_counts <- full$M
    }
    if(second[1] == "length") {
      if(second[2] != "no") {
        norm_counts <- withinLaneNormalization(norm_counts, full$annot$length, which = second[2])  
      }
    }
  } else {
    norm_counts <- full$M
  }

  if(STEP3 == "no") {
    norm_counts <- full$M
  } else if (STEP3 == "tmm") {
    norm_counts <- tmm(norm_counts, long = 1000, lc = 0, k = 0)
  } else {
    norm_counts <- betweenLaneNormalization(norm_counts, which = STEP3, offset = FALSE)
  }
  cat("Normalization done. Step1: ", STEP1, ", Step2: ", STEP2, " Step3: ", STEP3, "\n")
  
  # Saving normalized data
  full$M <- norm_counts
  
  norm_data_cpm10 <- filtered.data(full$M, factor=full$targets$group, 
                                   norm=TRUE, cv.cutoff=100, cpm=10)
  
  filtered <- nrow(full$M) - nrow(norm_data_cpm10)
  cat("After normalization. There are", filtered, "genes with counts per million mean < 10", 
      nrow(norm_data_cpm10), "with counts per million mean > 10 \n")
  
  full <-list(M = norm_data_cpm10, 
                    annot = full$annot %>% filter(gene_id %in% rownames(norm_data_cpm10)),
                    targets = full$targets)
  
  stopifnot(nrow(full$M) == nrow(full$annot))
  stopifnot(rownames(full$M) == rownames(full$annot))
  
  cat("Saving", paste(STEP1, STEP2, STEP3, "norm_full.RData", sep = "_"), "\n")
  save(full, file=snakemake@output[["norm_rdata"]], compress="xz")
  
  cat("Generating data matrices for Aracne\n")
  ## Data matrices for Aracne

  #normal samples
  normal <- as_tibble(full$M[, full$targets$group == "normal"])
  normal <- bind_cols(gene=as.character(full$annot$gene_id), normal)
  
  #cancer samples
  cancer <- as_tibble(full$M[, full$targets$group == "cancer"])
  cancer <- bind_cols(gene=as.character(full$annot$gene_id), cancer)
  
  #gene_ids
  symbols <- full$annot %>% select(gene_id)
  
  cat("Saving data\n")
  
  cat("Saving", paste(STEP1, STEP2, STEP3, "norm_normal.tsv", sep = "_"), "\n")
  write_tsv(normal, file=snakemake@output[["normal_matrix"]])
  
  cat("Saving", paste(STEP1, STEP2, STEP3, "norm_cancer.tsv", sep = "_"), "\n")
  write_tsv(cancer, file=snakemake@output[["cancer_matrix"]])
  
  cat("Saving", paste(STEP1, STEP2, STEP3, "norm_genelist.txt", sep = "_"), "\n")
  write_tsv(symbols, file=snakemake@output[["gene_list"]], col_names = F)
}##########################################
#########################################