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
###############################################################################
library("BiocParallel")
library("parallel")
library("NOISeq")
library("EDASeq")
# register(SnowParam(workers=detectCores()-1, progress=TRUE))#Windows
register(MulticoreParam(workers=detectCores()-1, progress=TRUE))#Linux
options(width=80)
###############################################################################
## Filter genes with low expression
cat("#################\n")
cat("Step 3: Normalization\n")
cat("#################\n")
DATADIR <- '/pipeline/data/'
RDATA <- paste(DATADIR, "rdata", sep = "")
PLOTSDIR <-paste(DATADIR, "plots", sep = "")
PLOTSNORMDIR <- paste(PLOTSDIR, "normalization", sep = "/")
dir.create(PLOTSNORMDIR)
w <- 1024
h <- 1024
p <- 24

##########################################################################
if(FALSE) {### NORMALIZATION METHODS TESTING
  ###########################################################################
  load(file=paste(RDATA, "RawFull.RData", sep="/"))
  { ### We keep only genes with mean expression count > 10 
  exp.genes <- apply(full$M, 1, function(x) mean(x)>10)
  egtable <- table(exp.genes)
  # FALSE  TRUE 
  #  2234 17215
  cat("There are", egtable[[1]], "genes with mean expression count < 10", egtable[[2]], "with mean count > 10 \n")
  
  ##Filtering low expression genes
  mean10 <- list(M=full$M[exp.genes, ], Annot=full$Annot[exp.genes, ], Targets=full$Targets)
  rownames(mean10$Annot) <- rownames(mean10$M)
  cat("Saving Mean10.RData \n") 
  save(mean10, file=paste(RDATA, "Mean10.RData", sep="/"), compress="xz")
  }
  cat("Testing normalization methods\n.")
  mydataM10EDA <- EDASeq::newSeqExpressionSet(
    counts=mean10$M,
    featureData=mean10$Annot,
    phenoData=data.frame(
      conditions=mean10$Targets$Group,
      row.names=colnames(mean10$M)))
  
  lenght.norm <- c("RPKM", "loess", "full")
  gc.norm <- c("loess", "full")
  between.nom <- c("tmm", "full")
  normalization.results <- data.frame()
  
  ## This function gets the relevant statistics for the regression methods for GC and Length bias
  getRegressionStatistics <- function(regressionmodel) {
    name <- names(regressionmodel)
    print(name)
    rsquared <- summary(regressionmodel[[1]])$r.squared
    print(rsquared)
    fstatistic <- summary(regressionmodel[[1]])$fstatistic
    pvalue <- signif(pf(q = fstatistic[1], df1 = fstatistic[2], df2 = fstatistic[3], lower.tail = FALSE), 2)
    return(list("name" = name, "r2" = rsquared, "p" = pvalue))
  }
  
  ## This function will retrieve the NOISeq results to test the normalization combination
  getNOISeqResults <- function(step1, step2, step3, n.counts, m10.data) {
    ### Check the NOISEq results 
    mydata <- NOISeq::readData(
      data = n.counts, 
      length = m10.data$Annot[, c("EnsemblID", "Length")], 
      biotype = m10.data$Annot[, c("EnsemblID", "Type")], 
      chromosome = m10.data$Annot[, c("Chr", "Start", "End")], 
      factors = m10.data$Targets[, "Group",drop=FALSE], 
      gc = m10.data$Annot[, c("EnsemblID", "GC")])
    nsamples <- dim(m10.data$Targets)[1]
    
    ### Length bias 
    mylengthbias <- dat(mydata, factor="Group", norm = TRUE, type="lengthbias")
    l.stats.1 <- getRegressionStatistics(mylengthbias@dat$RegressionModels[1])
    l.stats.2 <- getRegressionStatistics(mylengthbias@dat$RegressionModels[2])
    
    png(paste(PLOTSNORMDIR, paste(step1, step2, step3, "Lenghtbias.png", sep = "_"), sep="/"), width=w, height=h, pointsize=p)
    explo.plot(mylengthbias, samples = NULL, toplot = "global")
    dev.off()
    
    ## GC Bias
    mygcbias <- dat(mydata, factor = "Group", norm = TRUE, type ="GCbias")
    gc.stats.1 <- getRegressionStatistics(mygcbias@dat$RegressionModels[1])
    gc.stats.2 <- getRegressionStatistics(mygcbias@dat$RegressionModels[2])
    
    png(paste(PLOTSNORMDIR, paste(step1, step2, step3, "GCbias.png", sep = "_"), sep="/"), width=w, height=h, pointsize=p)
    explo.plot(mygcbias, samples = NULL, toplot = "global")
    dev.off()
    
    dtable <- data.frame(PASSED = NA)
    #RNA Composition
    rnacomp <-  tryCatch({
      myrnacomp <- dat(mydata, norm = TRUE, type="cd")
      dtable <- table(myrnacomp@dat$DiagnosticTest[,  "Diagnostic Test"])
      if (is.na(dtable["PASSED"])) dtable <- data.frame(PASSED = 0)
      png(paste(PLOTSNORMDIR,  paste(step1, step2, step3, "RNAComposition.png", sep = "_"), sep="/"), width=w, height=h, pointsize=p)
      explo.plot(myrnacomp, samples = 1:12)
      dev.off()
    }, error = function(cond) {
      cat("Error running rnacomposition test.\n")
      return(NA)
    })
    
    norm.set.results <- data.frame(step1, step2, step3, 
                                   l.stats.1$r2, l.stats.1$p, l.stats.2$r2, l.stats.2$p,
                                   gc.stats.1$r2, gc.stats.1$p, gc.stats.2$r2, gc.stats.2$p, 
                                   dtable["PASSED"], dtable["PASSED"]/nsamples)
    colnames(norm.set.results) <- c("Step1", "Step2", "Step3", 
                                    paste("Lenght", l.stats.1$name, "R2", sep = "."), paste("Lenght", l.stats.1$name, "p-value", sep = "."),  
                                    paste("Lenght", l.stats.2$name, "R2", sep = "."), paste("Lenght", l.stats.2$name, "p-value", sep = "."), 
                                    paste("GC", gc.stats.1$name, "R2", sep = "."), paste("GC", gc.stats.1$name, "p-value", sep = "."), 
                                    paste("GC", gc.stats.2$name, "R2", sep = "."), paste("GC", gc.stats.2$name, "p-value", sep = "."),
                                    "RNA.PassedSamples", "RNA.PassedProportion")
    return(norm.set.results)
  } 
  
  ## We try with length normalization first
  for (ln in lenght.norm) {
    if (ln == "RPKM") {
      ln.data <- rpkm(as.matrix(counts(mydataM10EDA)), long=mean10$Annot$Length)
    } else {
      ln.data <- withinLaneNormalization(counts(mydataM10EDA), mean10$Annot$Length, which = ln)
    }
    for (gcn in gc.norm) {
      gcn.data <- withinLaneNormalization(ln.data, mean10$Annot$GC, which = gcn)
      for (bn in between.nom) {
        if (bn == "tmm") {
          between.data <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
        } else {
          between.data <- betweenLaneNormalization(gcn.data, which = bn, offset = FALSE)
        }
        cat("Testing with length normalization: ", ln, ", GC normalization: ", gcn, " and between lanes normalization: ", bn, "\n")
        norm.noiseq.results <- getNOISeqResults(paste("Length", ln, sep = "."), paste("GC", gcn, sep = "."), paste("Between", bn, sep =  "."),
                                                between.data, mean10)
        normalization.results <- rbind(normalization.results, norm.noiseq.results)                      
      }
    }
  }
  
  ## Now we try with GC normalization first
  for (gcn in gc.norm) {
    gcn.data <- withinLaneNormalization(counts(mydataM10EDA), mean10$Annot$GC, which = gcn)
    for (ln in lenght.norm) {
      if (ln == "RPKM") {
        ln.data <- rpkm(gcn.data, long=mean10$Annot$Length)
      } else {
        ln.data <- withinLaneNormalization(gcn.data, mean10$Annot$Length, which = ln)
      }
      for (bn in between.nom) {
        if (bn == "tmm") {
          between.data <- tmm(ln.data, long = 1000, lc = 0, k = 0)
        } else {
          between.data <- betweenLaneNormalization(ln.data, which = bn, offset = FALSE)
        }
        cat("Testing with GC normalization: ", gcn, ",  length normalization: ", ln, " and between lane normalization: ", bn, "\n")
        norm.noiseq.results <- getNOISeqResults(paste("GC", gcn, sep = "."), paste("Length", ln, sep = "."), paste("Between", bn, sep =  "."), 
                                                between.data, mean10)
        normalization.results <- rbind(normalization.results, norm.noiseq.results)                      
      }
    }
  }
  cat("End of normalization texting\n")
  cat("Saving NormalizationResults.tsv\n") 
  write.table(normalization.results, file=paste(RDATA, "NormalizationResults.tsv", sep="/"), 
              quote = F, sep = "\t", row.names = F)
} else {
  {##### USER SELECTED NORMALIZATION
  load(file=paste(RDATA, "Mean10.RData", sep="/"))
  ln.data <- withinLaneNormalization(mean10$M, mean10$Annot$Length, which = "full")
  gcn.data <- withinLaneNormalization(ln.data, mean10$Annot$GC, which = "full")
  norm.counts <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
  step1 <- "Length.full"
  step2 <- "GC.full"
  step3 <- "Between.tmm"
  ##### After normalization is done
  ##### Step1, step2, step3
  }##########################################
  ## EXPLORATORY ANALYSIS (NOISeq package)
  ##########################################
  {  
  mydata <- NOISeq::readData(
    data = norm.counts, 
    length = mean10$Annot[, c("EnsemblID", "Length")], 
    biotype = mean10$Annot[, c("EnsemblID", "Type")], 
    chromosome = mean10$Annot[, c("Chr", "Start", "End")], 
    factors = mean10$Targets[, "Group",drop=FALSE], 
    gc = mean10$Annot[, c("EnsemblID", "GC")])
  
  ### Length bias 
  mylengthbias <- dat(mydata, factor="Group", norm = TRUE, type="lengthbias")
  png(paste(PLOTSDIR, paste(step1, step2, step3, "corrected", "Lenghtbias.png", sep = "_"), sep="/"), width=w, height=h, pointsize=p)
  explo.plot(mylengthbias, samples = NULL, toplot = "global")
  dev.off()
  cat("Lenght bias plot generated\n")
  
  ## GC Bias
  mygcbias <- dat(mydata, factor = "Group", norm = TRUE, type ="GCbias")
  png(paste(PLOTSDIR, paste(step1, step2, step3, "corrected", "GCbias.png", sep = "_"), sep="/"), width=w, height=h, pointsize=p)
  explo.plot(mygcbias, samples = NULL, toplot = "global")
  dev.off()
  cat("GC bias plot generated\n")
  
  #RNA Composition
  myrnacomp <- dat(mydata, norm = TRUE, type="cd")
  png(paste(PLOTSDIR,  paste(step1, step2, step3, "corrected", "RNAComposition.png", sep = "_"), sep="/"), width=w, height=h, pointsize=p)
  explo.plot(myrnacomp, samples = 1:12)
  dev.off()
  cat("RNA composition plot generated\n")
  
  mycountsbio <- dat(mydata, factor = NULL, type = "countsbio")
  png(paste(PLOTSDIR, paste(step1, step2, step3, "corrected", "protein_coding_boxplot.png", sep = "_"), sep="/"), width=w*2, height=h, pointsize=p)
  explo.plot(mycountsbio, toplot = "protein_coding", samples = NULL, plottype = "boxplot")
  dev.off()
  cat("Counts distribution plot for protein coding and all samples generated\n")
  
  png(paste(PLOTSDIR, paste(step1, step2, step3, "corrected", "protein_coding_barplot.png", sep = "_"), sep="/"), width=w*2, height=h, pointsize=p)
  explo.plot(mycountsbio, toplot = "protein_coding", samples = NULL, plottype = "barplot")
  dev.off()
  cat("Counts distribution barplot for protein coding biotype and all samples generated\n")
  
  # Saving normalized data
  norm.data <- mean10
  norm.data$M <- norm.counts
  
  cat("Saving", paste(step1, step2, step3, "Norm.RData", sep = "_"), "\n")
  save(norm.data, file=paste(RDATA, paste(step1, step2, step3, "Norm.RData", sep = "_"), sep="/"), compress="xz")
  }##########################################
  ## LOW EXPRESSION FILTERING
  ##########################################
  {
  library("ggplot2")
  library("reshape2")
  cat("Counts per million < 10 filter\n")
  ##Density plot
  pl<-ggplot(data=melt(log(norm.data$M+1)), aes(x=value, group=Var2, colour=Var2))+geom_density(show.legend = F)
  pdf(file=paste(PLOTSDIR, paste(step1, step2, step3, "corrected_densitylog.pdf", sep = "_"), sep="/"))
  pl
  dev.off()
  cat("Density plot for all samples generated\n")
  
  norm.data.cpm10 <- filtered.data(norm.data$M, factor=norm.data$Targets$Group, 
                             norm=TRUE, cv.cutoff=100, cpm=10)

  filtered <- nrow(norm.data$M)-nrow(norm.data.cpm10)
  cat("After normalization. There are", filtered, "genes with counts per million mean < 10", 
      nrow(norm.data.cpm10), "with counts per million mean > 10 \n")
  
  norm.data.cpm10 <-list(M = norm.data.cpm10, 
                    Annot = norm.data$Annot[row.names(norm.data$Annot) %in% row.names(norm.data.cpm10),],
                    Targets = norm.data$Targets)
  
  stopifnot(nrow(norm.data.cpm10$M) == nrow(norm.data.cpm10$Annot))
  stopifnot(row.names(norm.data.cpm10$M) == row.names(norm.data.cpm10$Annot))
  pl<-ggplot(data=melt(log(norm.data.cpm10$M+1)), aes(x=value, group=Var2, colour=Var2))+geom_density(show.legend = F)
  pdf(file=paste(PLOTSDIR, paste(step1, step2, step3, "corrected_cpm10_densitylog.pdf", sep = "_"), sep="/"))
  pl
  dev.off()
  cat("Density plot for all samples after cpm10 filter generated\n")
  
  cat("Saving", paste(step1, step2, step3, "Norm_cpm10.RData", sep = "_"), "\n")
  save(norm.data.cpm10, file=paste(RDATA, paste(step1, step2, step3, "Norm_cpm10.RData", sep = "_"), sep="/"), compress="xz")
  }##########################################
  ## PCA
  ##########################################
  {
  ### PCA with cpm (filter but no ARSyn)
  pca.results <- pca.results <- PCA.GENES(t(log2(1 + norm.data.cpm10$M)))
  
  ## Variance explained by each component
  pdf(file=paste(PLOTSDIR, paste(step1, step2, step3, "corrected_cpm10_PCAVariance.pdf", sep = "_"), sep="/"), 
      width = 4*2, height = 4*2)
  barplot(pca.results$var.exp[,1], xlab = "PC", ylab = "Explained variance")
  dev.off()
  cat("PCA variance norm plot generated.\n")
  
  ## Loading plot
  pdf(file=paste(PLOTSDIR, paste(step1, step2, step3, "corrected_cpm10_PCALoading.pdf", sep = "_"), sep="/"), 
      width = 4*2, height = 4*2)
  plot(pca.results$loadings[,1:2], col = 1, pch = 20, cex = 0.5,
       xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
       ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
       main = "PCA loadings",
       xlim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1),
       ylim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1))  
  dev.off()
  cat("PCA loading norm plot generated.\n")
  
  ## Score plot
  mycol <- as.character(norm.data.cpm10$Targets$Group)
  mycol[mycol == 'N'] <- "black"
  mycol[mycol == 'T'] <- "red2"
  
  pdf(file=paste(PLOTSDIR,  paste(step1, step2, step3, "corrected_cpm10_PCAScore.pdf", sep = "_"), sep="/"), 
      width = 5*2, height = 5)
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
  legend("topright", c("N", "T"), col = c("black", "red2"), ncol = 2, pch = 1)
  
  # PC1 & PC3
  rango2 = diff(range(pca.results$scores[,c(1,3)]))
  plot(pca.results$scores[,c(1,3)], col = "white",
       xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
       ylab = paste("PC 3 ", round(pca.results$var.exp[3,1]*100,0), "%", sep = ""),
       main = "PCA scores",
       xlim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1),
       ylim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1))
  points(pca.results$scores[,1], pca.results$scores[,3], col = mycol, cex = 1.5)
  legend("topright", c("N", "T"), col = c("black", "red2"), ncol = 2, pch = 1)
  dev.off()
  cat("PCA scores norm plot generated.\n")
  }##########################################
  ## ARSyN to reduce batch effect
  ##########################################
  {
  mydata <- NOISeq::readData(
    data = norm.data.cpm10$M, 
    factors = norm.data.cpm10$Targets[, "Group",drop = FALSE])
  
  cat("Performing ARSyN for batch correction")
  myARSyN <- ARSyNseq(mydata, norm = "n", logtransf = FALSE)
  pca.dat <- dat(myARSyN, type = "PCA", logtransf = F)
  pca.results <- pca.dat@dat$result
  
  ## Variance explained by each component
  pdf(file=paste(PLOTSDIR, paste(step1, step2, step3, "corrected_cpm10_arsyn_PCAVariance.pdf", sep = "_"), sep="/"),
      width = 4*2, height = 4*2)
  barplot(pca.results$var.exp[,1], xlab = "PC", ylab = "Explained variance", ylim = c(0,0.4))
  dev.off()
  cat("PCA variance arsyn plot generated.\n")
  
  ## Loading plot
  pdf(file=paste(PLOTSDIR, paste(step1, step2, step3, "corrected_cpm10_arsyn_PCALoading.pdf", sep = "_"), sep="/"), 
      width = 4*2, height = 4*2)
  plot(pca.results$loadings[,1:2], col = 1, pch = 20, cex = 0.5,
       xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
       ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
       main = "PCA loadings",
       xlim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1),
       ylim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1))  
  dev.off()
  cat("PCA loading arsyn plot generated.\n")
  
  ## Score plot
  mycol <- as.character(norm.data.cpm10$Targets$Group)
  mycol[mycol == 'N'] <- "black"
  mycol[mycol == 'T'] <- "red2"
  
  pdf(file=paste(PLOTSDIR, paste(step1, step2, step3, "corrected_cpm10_arsyn_PCAScoreARSyN.pdf", sep = "_"), sep="/"), 
      width = 5*2, height = 5)
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
  legend("topright", c("N", "T"), col = c("black", "red2"), ncol = 2, pch = 1)
  
  # PC1 & PC3
  rango2 = diff(range(pca.results$scores[,c(1,3)]))
  plot(pca.results$scores[,c(1,3)], col = "white",
       xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
       ylab = paste("PC 3 ", round(pca.results$var.exp[3,1]*100,0), "%", sep = ""),
       main = "PCA scores",
       xlim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1),
       ylim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1))
  points(pca.results$scores[,1], pca.results$scores[,3], col = mycol, cex = 1.5)
  legend("topright", c("N", "T"), col = c("black", "red2"), ncol = 2, pch = 1)
  dev.off()
  cat("PCA scores arsyn plot generated.\n")
  
  pl<-ggplot(data=melt(log(assayData(myARSyN)$exprs+1)), aes(x=value, group=Var2, colour=Var2))+geom_density(show.legend = F)
  pdf(file=paste(PLOTSDIR, paste(step1, step2, step3, "corrected_cpm10_arsyn_densitylog.pdf", sep = "_"), sep="/"))
  pl
  dev.off()
  
  cat('ARSyN data. Final dimensions: ', paste(dim(assayData(myARSyN)$exprs), collapse=", "), '.\n')
  
  ##Saving everything
  norm.data_cpm10_arsyn <- list(M = assayData(myARSyN)$exprs, Annot = norm.data.cpm10$Annot, 
                                Targets = norm.data.cpm10$Targets)

  stopifnot(nrow(norm.data_cpm10_arsyn$M) == nrow(norm.data_cpm10_arsyn$Annot))
  stopifnot(all(row.names(norm.data_cpm10_arsyn$M) == row.names(norm.data_cpm10_arsyn$Annot)))
  
  cat("Saving", paste(step1, step2, step3, "Norm_cpm10_arsyn.RData", sep = "_"), "\n")
  save(norm.data_cpm10_arsyn, file=paste(RDATA, paste(step1, step2, step3, "Norm_cpm10_arsyn.RData", sep = "_"), sep="/"), compress="xz")
  
  cat("Generating data matrices with arsyn for Aracne\n")
  ## Data matrices for Aracne
  ## ALL = healthy | cancer
  M <- as.data.frame(norm.data_cpm10_arsyn$M)
  M <- cbind(gene=as.character(norm.data_cpm10_arsyn$Annot$EnsemblID), M)
  
  #normal samples
  normal <- as.data.frame(norm.data_cpm10_arsyn$M[,norm.data_cpm10_arsyn$Targets$Group == "N"])
  normal <- cbind(gene=as.character(norm.data_cpm10_arsyn$Annot$EnsemblID), normal)
  
  #tumor samples
  tumor <- as.data.frame(norm.data_cpm10_arsyn$M[,norm.data_cpm10_arsyn$Targets$Group == "T"])
  tumor <- cbind(gene=as.character(norm.data_cpm10_arsyn$Annot$EnsemblID), tumor)
  
  #EnsemblIDs
  symbols <-as.character(norm.data_cpm10_arsyn$Annot$EnsemblID)
  
  cat("Saving arsyn data\n")
  
  cat("Saving", paste(step1, step2, step3, "Norm_cpm10_arsyn_all.tsv", sep = "_"), "\n")
  write.table(M, file=paste(RDATA, paste(step1, step2, step3, "Norm_cpm10_arsyn_all.tsv", sep = "_"), sep="/"), 
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("Saving", paste(step1, step2, step3, "Norm_cpm10_arsyn_normal.tsv", sep = "_"), "\n")
  write.table(normal, file = paste(RDATA, paste(step1, step2, step3, "Norm_cpm10_arsyn_normal.tsv", sep = "_"), sep="/"), 
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("Saving", paste(step1, step2, step3, "Norm_cpm10_arsyn_tumor.tsv", sep = "_"), "\n")
  write.table(tumor, file = paste(RDATA, paste(step1, step2, step3, "Norm_cpm10_arsyn_tumor.tsv", sep = "_"), sep="/"), 
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("Saving", paste(step1, step2, step3, "Norm_cpm10_genelist.txt", sep = "_"), "\n")
  write.table(symbols, file = paste(RDATA, paste(step1, step2, step3, "Norm_cpm10_genelist.txt", sep = "_"), sep="/"), 
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  cat("Saving", paste(step1, step2, step3, "Norm_cpm10_arsyn_all.RData", sep = "_"), "\n")
  save(norm.data_cpm10_arsyn, file=paste(RDATA, paste(step1, step2, step3, "Norm_cpm10_arsyn.RData", sep = "_"), sep="/"), 
       compress="xz")
  }
}
#############################
##########################################################################
## GREAT JOB!!! YOU MADE IT TILL THE END!!!!
###########################################################################
cat("#################\n")
cat("End of Step 3: NORMALIZATION\n")
cat("#################\n")

