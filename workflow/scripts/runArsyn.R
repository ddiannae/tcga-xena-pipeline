log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(NOISeq)
library(readr)
library(dplyr)

STEP1 = snakemake@params[["step1"]]
STEP2 = snakemake@params[["step2"]]
STEP3 = snakemake@params[["step3"]]
PLOTSDIR <-paste(snakemake@params[["tissue_dir"]], "plots","arsyn", sep="/")
dir.create(PLOTSDIR, recursive = TRUE)

w <- 1024
h <- 1024
p <- 24
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
  pca.dat <- dat(myARSyN, type = "PCA", logtransf = F)
  pca.results <- pca.dat@dat$result
  
  ## Variance explained by each component
  png(file=paste(PLOTSDIR, paste(STEP1, STEP2, STEP3,  "corrected_cpm10_arsyn_pca_variance.png", sep = "_"), sep="/"),
      width = w, height = h, pointsize = p)
  barplot(pca.results$var.exp[,1], xlab = "PC", ylab = "Explained variance", ylim = c(0,0.4))
  dev.off()
  cat("PCA variance arsyn plot generated.\n")
  
  ## Loading plot
  png(file=paste(PLOTSDIR, paste(STEP1, STEP2, STEP3,  "corrected_cpm10_arsyn_pca_loading.png", sep = "_"), sep="/"), 
      width = w, height = h, pointsize = p)
  plot(pca.results$loadings[,1:2], col = 1, pch = 20, cex = 0.5,
       xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
       ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
       main = "PCA loadings",
       xlim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1),
       ylim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1))  
  dev.off()
  cat("PCA loading arsyn plot generated.\n")
  
  ## Score plot
  mycol <- as.character(full$targets$group)
  mycol[mycol == "normal"] <- "black"
  mycol[mycol == "cancer"] <- "red2"
  
  png(file=paste(PLOTSDIR, paste(STEP1, STEP2, STEP3,  "corrected_cpm10_arsyn_pca_score.png", sep = "_"), sep="/"), 
      width = w*2, height = h, pointsize = p)
  par(mfrow = c(1,2))
  
  # PC1 & PC2
  rango <- diff(range(pca.results$scores[,1:2]))
  plot(pca.results$scores[,1:2], col = "white",
       xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
       ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
       main = "PCA scores",
       xlim = range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1),
       ylim = range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1))
  points(pca.results$scores[,1], pca.results$scores[,2], col = mycol, cex = 1.5)  
  legend("topright", c("normal", "cancer"), col = c("black", "red2"), ncol = 2, pch = 1)
  
  # PC1 & PC3
  rango2 = diff(range(pca.results$scores[,c(1,3)]))
  plot(pca.results$scores[,c(1,3)], col = "white",
       xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
       ylab = paste("PC 3 ", round(pca.results$var.exp[3,1]*100,0), "%", sep = ""),
       main = "PCA scores",
       xlim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1),
       ylim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1))
  points(pca.results$scores[,1], pca.results$scores[,3], col = mycol, cex = 1.5)
  legend("topright", c("normal", "cancer"), col = c("black", "red2"), ncol = 2, pch = 1)
  dev.off()
  cat("PCA scores arsyn plot generated.\n")

  cat('ARSyN data. Final dimensions: ', paste(dim(assayData(myARSyN)$exprs), collapse=", "), '.\n')
  
  ##Saving everything
  full_arsyn <- list(M = assayData(myARSyN)$exprs, annot = full$annot, 
                                targets = full$targets)
  
  stopifnot(nrow(full$M) == nrow(full_arsyn$annot))
  stopifnot(all(row.names(full_arsyn$M) == row.names(full_arsyn$annot)))
  
  full <- full_arsyn
  cat("Saving", paste(STEP1, STEP2, STEP3,  "norm_cpm10_arsyn.RData", sep = "_"), "\n")
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
  
  cat("Saving", paste(STEP1, STEP2, STEP3,  "norm_cpm10_arsyn_normal.tsv", sep = "_"), "\n")
  write_tsv(normal, file=snakemake@output[["normal_matrix"]])
  
  cat("Saving", paste(STEP1, STEP2, STEP3, "norm_cpm10_arsyn_cancer.tsv", sep = "_"), "\n")
  write_tsv(cancer, file=snakemake@output[["cancer_matrix"]])
  
  cat("Saving", paste(STEP1, STEP2, STEP3,  "norm_cpm10_genelist.txt", sep = "_"), "\n")
  write_tsv(symbols, file=snakemake@output[["gene_list"]], col_names = F)
}
