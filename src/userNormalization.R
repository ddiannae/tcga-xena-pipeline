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
##Usefull Libraries
library(ggplot2)
library(reshape2)
library(readr)
library(dplyr)
library(NOISeq)
library(EDASeq)
library(DESeq2)

args <- commandArgs(trailingOnly = T)

if (length(args) < 4 ) {
  stop("Incorrect number of arguments", call.=FALSE)
} else {
  TISSUE = args[1]
  DATADIR = args[2]
  STEP1 = args[3]
  STEP2 = args[4]
  STEP3 = args[5]
}
DATADIR <- paste(DATADIR, TISSUE, sep="/")
PLOTSDIR <-paste(DATADIR, "plots", sep="/")
RDATA <-paste(DATADIR, "rdata", sep="/")
PLOTSNORMDIR <- paste(PLOTSDIR, "normalization", sep = "/")
w <- 1024
h <- 1024
p <- 24
##########################################################################
{##### USER SELECTED NORMALIZATION
  cat("Loading data\n")
  load(file=paste(RDATA, "mean10_proteinCoding.RData", sep="/"))
  cat("Performing normalization. Step1: ", STEP1, ", Step2: ", STEP2, " Step3: ", STEP3, "\n")
  
  first <- strsplit(STEP1, "-")[[1]]
  second <- strsplit(STEP2, "-")[[1]]
  
  if(first[1] == "length") {
    norm_counts <- withinLaneNormalization(mean10$M, mean10$annot$length, which = first[2])
    if(second[1] == "gc") {
      norm_counts <- withinLaneNormalization(norm_counts, mean10$annot$gc, which = second[2])  
    }
  } else if (first[1] == "gc") {
    norm_counts <- withinLaneNormalization(mean10$M, mean10$annot$gc, which = first[2])
    if(second[1] == "length") {
      norm_counts <- withinLaneNormalization(norm_counts, mean10$annot$length, which = second[2])  
    }
  } else {
    norm_counts <- mean10$M
  }

  if (STEP3 == "tmm") {
    norm_counts <- tmm(norm_counts, long = 1000, lc = 0, k = 0)
  } else {
    norm_counts <- betweenLaneNormalization(norm_counts, which = STEP3, offset = FALSE)
  }
  cat("Normalization done. Step1: ", STEP1, ", Step2: ", STEP2, " Step3: ", STEP3, "\n")
  
}##########################################
## EXPLORATORY ANALYSIS (NOISeq package)
##########################################
{  
  mydata <- NOISeq::readData(
    data = norm_counts, 
    length = mean10$annot %>% select(gene_id, length), 
    biotype = mean10$annot %>% select(gene_id, gene_type), 
    chromosome = mean10$annot %>% select(chr, start, end), 
    factors = mean10$target %>% select(group),
    gc = mean10$annot %>% select(gene_id, gc))
  
  ### Length bias 
  mylengthbias <- dat(mydata, factor="group", norm = TRUE, type = "lengthbias")
  png(paste(PLOTSDIR, paste(STEP1, STEP2, STEP3, "corrected_length_bias.png", sep = "_"), sep = "/"),
      width=w, height=h, pointsize=p)
  explo.plot(mylengthbias, samples = NULL, toplot = "global")
  dev.off()
  cat("Length bias plot generated\n")
  
  ## GC Bias
  mygcbias <- dat(mydata, factor = "group", norm = TRUE, type = "GCbias")
  png(paste(PLOTSDIR, paste(STEP1, STEP2, STEP3, "corrected_gc_bias.png", sep = "_"), sep = "/"), 
      width=w, height=h, pointsize=p)
  explo.plot(mygcbias, samples = NULL, toplot = "global")
  dev.off()
  cat("GC bias plot generated\n")
  
  #RNA Composition
  myrnacomp <- dat(mydata, norm = TRUE, type = "cd")
  dtable <- table(myrnacomp@dat$DiagnosticTest[,  "Diagnostic Test"])
  if (is.na(dtable["PASSED"])) dtable <- data.frame(PASSED = 0)
  cat("Passed:", dtable["PASSED"], "samples. Proportion: ", 
      dtable["PASSED"]/( dtable["PASSED"] +  dtable["FAILED"]), "\n")
  
  png(paste(PLOTSDIR,  paste(STEP1, STEP2, STEP3,  "corrected_rna_composition.png", sep = "_"), sep = "/"), 
      width=w, height=h, pointsize=p)
  explo.plot(myrnacomp, samples = 1:12)
  dev.off()
  cat("RNA composition plot generated\n")
  
  mycountsbio <- dat(mydata, factor = NULL, type = "countsbio")
  png(paste(PLOTSDIR, paste(STEP1, STEP2, STEP3,  "corrected", "protein_coding_boxplot.png", sep = "_"), sep="/"), 
      width=w*2, height=h, pointsize=p)
  explo.plot(mycountsbio, toplot = "protein_coding", samples = NULL, plottype = "boxplot")
  dev.off()
  cat("Counts distribution plot for protein coding and all samples generated\n")
  
  png(paste(PLOTSDIR, paste(STEP1, STEP2, STEP3,  "corrected", "protein_coding_barplot.png", sep = "_"), sep="/"), 
      width=w*2, height=h, pointsize=p)
  explo.plot(mycountsbio, toplot = "protein_coding", samples = NULL, plottype = "barplot")
  dev.off()
  cat("Counts distribution barplot for protein coding biotype and all samples generated\n")
  
  mycountsbio <- dat(mydata, factor = "group", type = "countsbio")
  png(paste(PLOTSDIR, paste(STEP1, STEP2, STEP3,  "corrected", "protein_coding_boxplot_group.png", sep = "_"), sep="/"), 
      width=w*2, height=h, pointsize=p)
  explo.plot(mycountsbio, toplot = "protein_coding", samples = NULL, plottype = "boxplot")
  dev.off()
  cat("Counts distribution plot for protein coding grouped generated\n")
  
  # Saving normalized data
  norm_data <- mean10
  norm_data$M <- norm_counts
  
  cat("Saving", paste(STEP1, STEP2, STEP3, "norm.RData", sep = "_"), "\n")
  save(norm_data, file=paste(RDATA, paste(STEP1, STEP2, STEP3,  "norm.RData", sep = "_"), sep="/"), compress="xz")
}##########################################
## LOW EXPRESSION FILTERING
##########################################
{
  cat("Counts per million < 10 filter\n")
  ##Density plot
  pl <- ggplot(data = melt(log(norm_data$M+1)), aes(x=value, group=Var2, colour=Var2)) + 
    geom_density(show.legend = F)
  
  png(file = paste(PLOTSDIR, paste(STEP1, STEP2, STEP3,  "corrected_densitylog.png", sep = "_"), sep="/"), 
      width = w, height = h, pointsize = p)
  print(pl)
  dev.off()
  cat("Density plot for all samples generated\n")
  
  norm_data_cpm10 <- filtered.data(norm_data$M, factor=norm_data$targets$group, 
                             norm=TRUE, cv.cutoff=100, cpm=10)

  filtered <- nrow(norm_data$M) - nrow(norm_data_cpm10)
  cat("After normalization. There are", filtered, "genes with counts per million mean < 10", 
      nrow(norm_data_cpm10), "with counts per million mean > 10 \n")
  
  norm_data_cpm10 <-list(M = norm_data_cpm10, 
                    annot = norm_data$annot[row.names(norm_data$annot) %in% row.names(norm_data_cpm10),],
                    targets = norm_data$targets)
  
  stopifnot(nrow(norm_data_cpm10$M) == nrow(norm_data_cpm10$annot))
  stopifnot(row.names(norm_data_cpm10$M) == row.names(norm_data_cpm10$annot))
  
  pl <- ggplot(data=melt(log(norm_data_cpm10$M+1)), aes(x=value, group=Var2, colour=Var2)) +
    geom_density(show.legend = F)
  png(file=paste(PLOTSDIR, paste(STEP1, STEP2, STEP3, "corrected_cpm10_densitylog.png", sep = "_"), sep="/"),
      width = w, height = h, pointsize = p)
  print(pl)
  dev.off()
  cat("Density plot for all samples after cpm10 filter generated\n")
  
  cat("Saving", paste(STEP1, STEP2, STEP3, "norm_cpm10.RData", sep = "_"), "\n")
  save(norm_data_cpm10, file=paste(RDATA, paste(STEP1, STEP2, STEP3, "norm_cpm10.RData", sep = "_"), sep="/"), compress="xz")
  
  cat("Generating data matrices for Aracne\n")
  ## Data matrices for Aracne
  ## ALL = healthy | cancer
  M <- as_tibble(norm_data_cpm10$M)
  M <- bind_cols(gene=as.character(norm_data_cpm10$annot$gene_id), M)
  
  #normal samples
  normal <- as_tibble(norm_data_cpm10$M[,norm_data_cpm10$targets$group == "normal"])
  normal <- bind_cols(gene=as.character(norm_data_cpm10$annot$gene_id), normal)
  
  #cancer samples
  cancer <- as_tibble(norm_data_cpm10$M[,norm_data_cpm10$targets$group == "cancer"])
  cancer <- bind_cols(gene=as.character(norm_data_cpm10$annot$gene_id), cancer)
  
  #gene_ids
  symbols <- as.character(norm_data_cpm10$annot$gene_id)
  
  cat("Saving data\n")
  
  cat("Saving", paste(STEP1, STEP2, STEP3,"norm_cpm10_all.tsv", sep = "_"), "\n")
  write.table(M, file=paste(RDATA, paste(STEP1, STEP2, STEP3,"norm_cpm10_all.tsv", sep = "_"), sep="/"), 
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("Saving", paste(STEP1, STEP2, STEP3, "norm_cpm10_normal.tsv", sep = "_"), "\n")
  write.table(normal, file = paste(RDATA, paste(STEP1, STEP2, STEP3, "norm_cpm10_normal.tsv", sep = "_"), sep="/"), 
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("Saving", paste(STEP1, STEP2, STEP3, "norm_cpm10_cancer.tsv", sep = "_"), "\n")
  write.table(cancer, file = paste(RDATA, paste(STEP1, STEP2, STEP3, "norm_cpm10_cancer.tsv", sep = "_"), sep="/"), 
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("Saving", paste(STEP1, STEP2, STEP3, "norm_cpm10_genelist.txt", sep = "_"), "\n")
  write.table(symbols, file = paste(RDATA, paste(STEP1, STEP2, STEP3, "norm_cpm10_genelist.txt", sep = "_"), sep="/"), 
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}##########################################
## PCA
##########################################
{
### PCA with cpm (filter but no ARSyn)
  pca.results <- pca.results <- PCA.GENES(t(log2(1 + norm_data_cpm10$M)))
  
  ## Variance explained by each component
  png(file=paste(PLOTSDIR, paste(STEP1, STEP2, STEP3, "corrected_cpm10_pca_variance.png", sep = "_"), sep="/"), 
      width = w, height = h, pointsize = p)
  barplot(pca.results$var.exp[,1], xlab = "PC", ylab = "Explained variance")
  dev.off()
  cat("PCA variance norm plot generated.\n")
  
  ## Loading plot
  png(file=paste(PLOTSDIR, paste(STEP1, STEP2, STEP3, "corrected_cpm10_pca_loading.png", sep = "_"), sep="/"), 
      width = w, height = h, pointsize = p)
  plot(pca.results$loadings[,1:2], col = 1, pch = 20, cex = 0.5,
       xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
       ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
       main = "PCA loadings",
       xlim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1),
       ylim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1))  
  dev.off()
  cat("PCA loading norm plot generated.\n")
  
  ## Score plot
  mycol <- as.character(norm_data_cpm10$targets$group)
  mycol[mycol == "normal"] <- "black"
  mycol[mycol == "cancer"] <- "red2"
  
  png(file=paste(PLOTSDIR,  paste(STEP1, STEP2, STEP3, "corrected_cpm10_pca_score.png", sep = "_"), sep="/"), 
      width = 2*w, height = h, pointsize = p)
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
  cat("PCA scores norm plot generated.\n")
}##########################################
## ARSyN to reduce batch effect
##########################################
{
  mydata <- NOISeq::readData(
    data = norm_data_cpm10$M, 
    factors = mean10$target %>% select(group))
  
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
  mycol <- as.character(norm_data_cpm10$targets$group)
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
  
  pl <- ggplot(data=melt(log(assayData(myARSyN)$exprs+1)), aes(x=value, group=Var2, colour=Var2)) + 
    geom_density(show.legend = F)
  
  png(file=paste(PLOTSDIR, paste(STEP1, STEP2, STEP3,  "corrected_cpm10_arsyn_densitylog.png", sep = "_"), sep="/"),
      width = w, height = h, pointsize = p)
  print(pl)
  dev.off()
  
  mydata <- NOISeq::readData(
    data =  assayData(myARSyN)$exprs, 
    length = norm_data_cpm10$annot %>% select(gene_id, length), 
    biotype = norm_data_cpm10$annot %>% select(gene_id, gene_type), 
    chromosome = norm_data_cpm10$annot %>% select(chr, start, end), 
    factors = norm_data_cpm10$target %>% select(group),
    gc = norm_data_cpm10$annot %>% select(gene_id, gc))

  mycountsbio <- dat(mydata, factor = "group", type = "countsbio")
  png(paste(PLOTSDIR, paste(STEP1, STEP2, STEP3,  "corrected_arsyn", "protein_coding_boxplot_group.png", sep = "_"), 
            sep="/"), width=w*2, height=h, pointsize=p)
    explo.plot(mycountsbio, toplot = "protein_coding", samples = NULL, plottype = "boxplot")
  dev.off()

  cat('ARSyN data. Final dimensions: ', paste(dim(assayData(myARSyN)$exprs), collapse=", "), '.\n')
  
  ##Saving everything
  norm_data_cpm10_arsyn <- list(M = assayData(myARSyN)$exprs, annot = norm_data_cpm10$annot, 
                                targets = norm_data_cpm10$targets)

  stopifnot(nrow(norm_data_cpm10_arsyn$M) == nrow(norm_data_cpm10_arsyn$annot))
  stopifnot(all(row.names(norm_data_cpm10_arsyn$M) == row.names(norm_data_cpm10_arsyn$annot)))
  
  cat("Saving", paste(STEP1, STEP2, STEP3,  "norm_cpm10_arsyn.RData", sep = "_"), "\n")
  save(norm_data_cpm10_arsyn, file=paste(RDATA, paste(STEP1, STEP2, STEP3,  "norm_cpm10_arsyn.RData", sep = "_"), 
                                         sep="/"), compress="xz")
  
  cat("Generating data matrices with arsyn for Aracne\n")
  ## Data matrices for Aracne
  ## ALL = healthy | cancer
  M <- as_tibble(norm_data_cpm10_arsyn$M)
  M <- bind_cols(gene=as.character(norm_data_cpm10_arsyn$annot$gene_id), M)
  
  #normal samples
  normal <- as_tibble(norm_data_cpm10_arsyn$M[,norm_data_cpm10_arsyn$targets$group == "normal"])
  normal <- bind_cols(gene=as.character(norm_data_cpm10_arsyn$annot$gene_id), normal)
  
  #cancer samples
  cancer <- as_tibble(norm_data_cpm10_arsyn$M[,norm_data_cpm10_arsyn$targets$group == "cancer"])
  cancer <- bind_cols(gene=as.character(norm_data_cpm10_arsyn$annot$gene_id), cancer)
  
  #gene_ids
  symbols <-as.character(norm_data_cpm10_arsyn$annot$gene_id)
  
  cat("Saving arsyn data\n")
  
  cat("Saving", paste(STEP1, STEP2, STEP3, "norm_cpm10_arsyn_all.tsv", sep = "_"), "\n")
  write.table(M, file=paste(RDATA, paste(STEP1, STEP2, STEP3,  "norm_cpm10_arsyn_all.tsv", sep = "_"), sep="/"), 
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("Saving", paste(STEP1, STEP2, STEP3,  "norm_cpm10_arsyn_normal.tsv", sep = "_"), "\n")
  write.table(normal, file = paste(RDATA, paste(STEP1, STEP2, STEP3, "norm_cpm10_arsyn_normal.tsv", sep = "_"), sep="/"), 
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("Saving", paste(STEP1, STEP2, STEP3, "norm_cpm10_arsyn_cancer.tsv", sep = "_"), "\n")
  write.table(cancer, file = paste(RDATA, paste(STEP1, STEP2, STEP3,  "norm_cpm10_arsyn_cancer.tsv", sep = "_"), sep="/"), 
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("Saving", paste(STEP1, STEP2, STEP3,  "norm_cpm10_genelist.txt", sep = "_"), "\n")
  write.table(symbols, file = paste(RDATA, paste(STEP1, STEP2, STEP3,  "norm_cpm10_genelist.txt", sep = "_"), sep="/"), 
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
 

