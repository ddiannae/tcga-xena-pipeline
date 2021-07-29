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
## Quality Control 
##    -EXPLORATORY ANALYSIS (NOISeq package)
##    -Reading data into NOISeq package -> mydata
##    -Plots
##        -Biodetection plot
##        -Count distribution per biotype
##        -Count distribution per sample
##        -Count distribution per Experimental factors
##    -Bias
##        -Length bias detection
##        -GC bias
##        -RNA composition
##    -PCA
##############################################################################
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(dplyr)
library(NOISeq)
library(ggplot2)

PLOTSDIR <-paste(snakemake@params[["tissue_dir"]], "plots", snakemake@params[["plots_type"]], sep="/")
dir.create(PLOTSDIR, recursive = TRUE)
w <- 1024
h <- 1024
p <- 24

load(snakemake@input[[1]])
##########################################
## EXPLORATORY ANALYSIS (NOISeq package)
##########################################
{
  ## Reading data into NOISeq package -> mydata
  mydata <- NOISeq::readData(
    data = full$M, 
    length = full$annot %>% select(gene_id, length) %>% as.data.frame(), 
    biotype = full$annot %>% select(gene_id, gene_type) %>% as.data.frame(), 
    chromosome = full$annot %>% select(chr, start, end) %>% as.data.frame(), 
    factors = full$targets %>% select(group) %>% as.data.frame(),
    gc = full$annot %>% select(gene_id, gc) %>% as.data.frame())

}
##########################################
## Plots
##########################################
{
  # Biodetection plot. Per group.
  mybiodetection <- dat(mydata, type="biodetection", factor="group", k=0)
  png(filename=paste(PLOTSDIR, "biodetection.Rd_%03d.png", sep="/"), 
      width=w, height=h, pointsize=p)
  explo.plot(mybiodetection)
  dev.off()
  cat("Biodetection plots generated\n")
  
  ## Count distribution per biotype. Using count per million, only for one sample
  mycountsbio <- dat(mydata, factor = NULL, type = "countsbio")
  png(filename=paste(PLOTSDIR, "countsbio.png", sep="/"), width=w, height=h, pointsize=p)
  explo.plot(mycountsbio, toplot = 1, samples = 1, plottype = "boxplot")
  dev.off()
  cat("Counts distribution plot per biotype and one sample generated\n")
  #What about expression level?

  ## Count distribution per sample
  mycountsbio <- dat(mydata, factor = NULL, type = "countsbio")
  png(paste(PLOTSDIR, "protein_coding_boxplot.png", sep="/"), width=w*2, height=h, pointsize=p)
  explo.plot(mycountsbio, toplot = "protein_coding",
      samples = NULL, plottype = "boxplot")
  dev.off()
  cat("Counts distribution plot for protein coding and all samples generated\n")

  png(paste(PLOTSDIR, "protein_coding_barplot.png", sep="/"), width=w*2, height=h, pointsize=p)
  explo.plot(mycountsbio, toplot = "protein_coding",
      samples = NULL, plottype = "barplot")
  dev.off()
  cat("Counts distribution barplot for protein coding biotype and all samples generated\n")

  mycountsbio <- dat(mydata, factor = "group", type = "countsbio")
  ## Count distribution per Experimental factors
  png(paste(PLOTSDIR, "protein_coding_boxplot_group.png", sep="/"), width=w, height=h, pointsize=p)
  explo.plot(mycountsbio, toplot = "protein_coding",
      samples = NULL, plottype = "boxplot")
  dev.off()
  cat("Counts distribution boxplot for protein coding biotype and group generated\n")

  png(paste(PLOTSDIR, "protein_coding_barplot_group.png", sep="/"),
      width=w, height=h, pointsize=p)
  explo.plot(mycountsbio, toplot = "protein_coding",
      samples = NULL, plottype = "barplot")
  dev.off()
  cat("Counts distribution barplot for protein coding biotype and group generated\n")
  # How much sensitivity we loose? 

}##########################################
## Bias
##########################################
{
  ## Length bias detection
  mylengthbias <- dat(mydata, factor="group", type="lengthbias")
  png(paste(PLOTSDIR, "length_bias.png", sep="/"), width=w, height=h, pointsize=p)
  explo.plot(mylengthbias, samples = NULL, toplot = "global")
  dev.off()
  cat("Lenght bias plot generated\n")
  # Do we see a clear pattern?

  ## GC bias
  mygcbias <- dat(mydata, factor = "group", type="GCbias")
  png(paste(PLOTSDIR, "gc_bias.png", sep="/"), width=w, height=h, pointsize=p)
  explo.plot(mygcbias, samples = NULL, toplot = "global")
  dev.off()
  cat("GC bias plot generated\n")
  # Do we see a clear pattern?

  ## RNA composition
  mycomp <- dat(mydata, type="cd")
  png(paste(PLOTSDIR, "rna_composition.png", sep="/"), width=w, height=h, pointsize=p)
  explo.plot(mycomp, samples=1:12)
  dev.off()
  cat("RNA composition plot generated\n")
  # Are samples comparable?
}
#############################
## PCA Analysis with NOISeq
##########################################
{
  pca.dat <- dat(mydata, type = "PCA", logtransf = F)
  pca.results <- pca.dat@dat$result
  
  ## Variance explained by each component
  png(file=paste(PLOTSDIR, "pca_variance.png", sep="/"),
      width = w, height = h, pointsize = p)
  barplot(pca.results$var.exp[,1], xlab = "PC", ylab = "Explained variance")
  dev.off()
  cat("PCA variance plot generated.\n")
  
  ## Loading plot
  png(file=paste(PLOTSDIR, "pca_loading.png", sep="/"), 
      width = w, height = h, pointsize = p)
  plot(pca.results$loadings[,1:2], col = 1, pch = 20, cex = 0.5,
       xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
       ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
       main = "PCA loadings",
       xlim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1),
       ylim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1))  
  dev.off()
  cat("PCA loading plot generated.\n")
  
  ## Score plot
  mycol <- as.character(full$targets$group)
  mycol[mycol == "normal"] <- "black"
  mycol[mycol == "cancer"] <- "red2"
  
  png(file=paste(PLOTSDIR, "pca_score.png", sep="/"), 
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
  cat("PCA scores plot generated.\n")
}
