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

load(file=paste(RDATA, "RawFull.RData", sep="/"))
exp.genes <- apply(full$M, 1, function(x) mean(x)>10)
egtable <- table(exp.genes)
# FALSE  TRUE 
#  2234 17215
cat("There are", egtable[[1]], "genes with mean expression count < 10", egtable[[2]], "with mean count > 10 \n")

##Filtering low expression genes
mean10 <- list(M=full$M[exp.genes, ], Annot=full$Annot[exp.genes, ], Targets=full$Targets)

save(mean10, file=paste(RDATA, "Mean10.RData", sep="/"), compress="xz")
cat("Saving RawFull.RData \n")
##########################################################################
#Normalization
cat("Testing normalization methods\n.")
rownames(mean10$Annot) <- rownames(mean10$M)
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

## This function will retrieve the NOISeq results to test the normalization combination
getNOISeqResults <- function(step1, step2, step3, norm.data, m10.data) {
  ### Check the NOISEq results 
  mydata <- NOISeq::readData(
    data = norm.data, 
    length = m10.data$Annot[, c("EnsemblID", "Length")], 
    biotype = m10.data$Annot[, c("EnsemblID", "Type")], 
    chromosome = m10.data$Annot[, c("Chr", "Start", "End")], 
    factors = m10.data$Targets[, "Group",drop=FALSE], 
    gc = m10.data$Annot[, c("EnsemblID", "GC")])
  nsamples <- dim(m10.data$Targets)[1]
  ### Length bias 
  mylengthbias <- dat(mydata, factor="Group", norm=T, type="lengthbias")
  regmodels <- mylengthbias@dat$RegressionModels
  name1 <- names(regmodels[1])
  rsquared.1 <- summary(regmodels[[1]])$r.squared
  fstatistic.1 <- summary(regmodels[[1]])$fstatistic
  pvalue.1 <- signif(pf(q = fstatistic.1[1], df1 = fstatistic.1[2], df2 = fstatistic.1[3], lower.tail = FALSE), 2)
  name2 <- names(regmodels[2])
  rsquared.2 <- summary(regmodels[[2]])$r.squared
  fstatistic.2 <- summary(regmodels[[2]])$fstatistic
  pvalue.2 <- signif(pf(q = fstatistic.2[1], df1 = fstatistic.2[2], df2 = fstatistic.2[3], lower.tail = FALSE), 2)
  
  png(paste(PLOTSNORMDIR, paste("Lengthbias_", step1, "_", step2, "_", step3, ".png", sep = ""), sep="/"), width=w, height=h, pointsize=p)
  explo.plot(mylengthbias, samples = NULL, toplot = "global")
  dev.off()
  
  ## GC Bias
  mygcbias <- dat(mydata, factor = "Group", norm = T, type ="GCbias")
  regmodels.gc <- mylengthbias@dat$RegressionModels
  name1.gc <- names(regmodels.gc[1])
  rsquared.1.gc <- summary(regmodels.gc[[1]])$r.squared
  fstatistic.1.gc <- summary(regmodels.gc[[1]])$fstatistic
  pvalue.1.gc <- signif(pf(q = fstatistic.1.gc[1], df1 = fstatistic.1.gc[2], df2 = fstatistic.1.gc[3], lower.tail = FALSE), 2)
  name2.gc <- names(regmodels.gc[2])
  rsquared.2.gc <- summary(regmodels.gc[[2]])$r.squared
  fstatistic.2.gc <- summary(regmodels.gc[[2]])$fstatistic
  pvalue.2.gc <- signif(pf(q = fstatistic.2.gc[1], df1 = fstatistic.2.gc[2], df2 = fstatistic.2.gc[3], lower.tail = FALSE), 2)
  
  png(paste(PLOTSNORMDIR, paste("GCbias_", step1, "_", step2, "_", step3, ".png", sep = ""), sep="/"), width=w, height=h, pointsize=p)
  explo.plot(mygcbias, samples = NULL, toplot = "global")
  dev.off()
  
  dtable <- data.frame(PASSED = NA)
  #RNA Composition
  rnacomp <-  tryCatch({
    myrnacomp <- dat(mydata, norm=F, type="cd")
    dtable <- table(myrnacomp@dat$DiagnosticTest[,  "Diagnostic Test"])
    png(paste(PLOTSNORMDIR,  paste("RNAComposition_", step1, "_", step2, "_", step3, ".png", sep = "") , sep="/"), width=w, height=h, pointsize=p)
    explo.plot(myrnacomp, samples = 1:12)
    dev.off()
  }, error = function(cond) {
    cat("Error running rnacomposition test.\n")
    return(NA)
  })
  
  norm.set.results <- data.frame(step1, step2, step3, 
                                 rsquared.1, pvalue.1, rsquared.2, pvalue.2,
                                 rsquared.1.gc, pvalue.1.gc, rsquared.2.gc, pvalue.2.gc, dtable["PASSED"]/nsamples)
  colnames(norm.set.results) <- c("Step1", "Step2", "Step3", 
                                  paste("Lenght", name1, "R2", sep = "."), paste("Lenght", name1, "p-value", sep = "."),  
                                  paste("Lenght", name2, "R2", sep = "."), paste("Lenght", name2, "p-value", sep = "."), 
                                  paste("GC", name1.gc, "R2", sep = "."), paste("GC", name1.gc, "p-value", sep = "."), 
                                  paste("GC", name2.gc, "R2", sep = "."), paste("GC", name2.gc, "p-value", sep = "."), "RNA.PassedSamples")
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
      cat("Testing with length normalization: ", ln, ", GC normalization: ", gcn, " and between lanes normalization: ", bn, "\n")
      if (bn == "tmm") {
        normalized.data <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
      } else {
        normalized.data <- betweenLaneNormalization(gcn.data, which = bn, offset = TRUE)
      }
  
      norm.noiseq.results <- getNOISeqResults(paste("Length", ln, sep = "."), paste("GC", gcn, sep = "."), paste("Between", bn, sep =  "."),
                                              normalized.data, mean10)
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
      cat("Testing with GC normalization: ", gcn, ",  length normalization: ", ln, " and between lane normalization: ", bn, "\n")
      if (bn == "tmm") {
        normalized.data <- tmm(ln.data, long = 1000, lc = 0, k = 0)
      } else {
        normalized.data <- betweenLaneNormalization(ln.data, which = bn, offset = TRUE)
      }
      
      norm.noiseq.results <- getNOISeqResults(paste("GC", ln, sep = "."), paste("Length", gcn, sep = "."), paste("Between", bn, sep =  "."), 
                                              normalized.data, mean10)
      normalization.results <- rbind(normalization.results, norm.noiseq.results)                      
    }
  }
}

write.table(normalization.results, file=paste(RDATA, "NormalizationResults.tsv", sep="/"), 
            quote = F, sep = "\t", row.names = F)
# ##withinLaneNormalization   #"loess","median","upper"
# ##betweenLaneNormalization  #"loess","median","upper"
# 
# #data2 <- withinLaneNormalization(as.matrix(counts(data)),feature$gc, which="full")
# #data2 <- withinLaneNormalization(as.matrix(counts(data)),feature$length, 
# #                                 which="full")
# #data2 <- withinLaneNormalization(data2,feature$length, 
# #                                 which="full")
# 
# ## base<-"RPKM_FULL-GC_TMM" En 0
# # data2 <- rpkm(as.matrix(counts(mydataM10EDA)), long=mean10$Annot$Length)
# # data2 <- withinLaneNormalization(data2,mean10$Annot$GC, which="full")
# # data3 <- tmm(data2, long = 1000, lc = 0, k = 0)
# 
# ## base<-"FULL-Length_FULL-GC_TMM" #En 1
# # data2 <- withinLaneNormalization(counts(mydataM10EDA),mean10$Annot$Length, which="full")
# # data2 <- withinLaneNormalization(data2,mean10$Annot$GC, which="full")
# # data3 <- tmm(data2, long = 1000, lc = 0, k = 0)
# 
# base<-"FULLGC_FULLLength_TMM" #En 2
# data2 <- withinLaneNormalization(counts(mydataM10EDA),mean10$Annot$GC, which="full")
# data2 <- withinLaneNormalization(data2,mean10$Annot$Length, which="full")
# data3 <- tmm(data2, long = 1000, lc = 0, k = 0)
# 
# ##base<-"FULL-GC_Loess-Length_TMM" #En 0
# # data2 <- withinLaneNormalization(counts(mydataM10EDA),mean10$Annot$GC, which="full")
# # data2 <- withinLaneNormalization(data2,mean10$Annot$Length, which="loess")
# # data3 <- tmm(data2, long = 1000, lc = 1, k = 0)
# 
# ## base<-"FULL-GC_FULL-Length_TMMlc1" #En 2
# # data2 <- withinLaneNormalization(counts(mydataM10EDA),mean10$Annot$GC, which="full")
# # data2 <- withinLaneNormalization(data2,mean10$Annot$Length, which="full")
# # data3 <- tmm(data2, long = 1000, lc = 1, k = 0)
# 
# ## base<-"FULL-GC_FULL-Length_BetweenFULL" #En 1
# # data2 <- withinLaneNormalization(counts(mydataM10EDA),mean10$Annot$GC, which="full")
# # data2 <- withinLaneNormalization(data2,mean10$Annot$Length, which="full")
# # data3 <- betweenLaneNormalization(data2, which="full",offset=TRUE)
# 
# ## base<-"FULL-GC_TMM_long-lc1" #En 3
# # data2 <- withinLaneNormalization(counts(mydataM10EDA),mean10$Annot$GC, which="full")
# # data3 <- tmm(data2, long = mean10$Annot$Length, lc = 1, k = 0)
# 
# ##base<-"FULL-GC_FULL-Length_TMM_V2" #En 2
# # data2 <- withinLaneNormalization(mydataM10EDA,"GC", which="full", num.bins=10, offset=FALSE, round=TRUE)
# # data2 <- withinLaneNormalization(data2,"Length", which="full")
# # data3 <- tmm(counts(data2), long = 1000, lc = 0)
# 
# 
# 
# ## Length bias detection
# 
# ## Length bias detection
# mylengthbias <- dat(mydata, factor="Group", norm=FALSE, type="lengthbias")
# png(paste(PLOTSDIR, "lengthbias_corrected.png", sep="/"), width=w, height=h, pointsize=p)
# explo.plot(mylengthbias, samples=1:2)
# dev.off()
# #Do we see a clear pattern?
# 
# ##GC bias
# mygcbiasRaw <- NOISeq::dat(mydata, factor = "Group", norm=FALSE, type="GCbias")
# png(paste(PLOTSDIR, "GCbias_corrected.png", sep="/"), width=w, height=h, pointsize=p)
# explo.plot(mygcbiasRaw)
# dev.off()
# 
# mycomp <- dat(mydata, norm=FALSE, type="cd")
# png(paste(PLOTSDIR, "RNAComposition_corrected.png", sep="/"), width=w, height=h, pointsize=p)
# explo.plot(mycomp, samples=1:12)
# dev.off()
# 
# table(mycomp@dat$DiagnosticTest[,  "Diagnostic Test"])
# 
# ##
# ##Metodo				Length 	GC	RNA FAILED/PASSED
# ##"RPKM_FULL-GC_TMM"			Si	No	880/0 Ondulaciones laterales
# ##"FULL-Length_FULL-GC_TMM"		No++	No	872/8
# ##"FULL-GC_FULL-Length_TMM"		No++	No	867/13  Este??  Este es el mejorcito en distribuciones
# ##"FULL-GC_Loess-Length_TMM" 		Si	Si	866/14
# ##"FULL-GC_FULL-Length_TMMlc1" 		Si	Si	868/12
# ##"FULL-GC_FULL-Length_BetweenFULL"	No--	No-- 	872/8   Este??   
# ##"FULL-GC_TMM_long-lc1"		Si	Si	874/6
# 
# ##QCreport(mydata, factor = "Group", norm=TRUE)
# ##Da error el reporte final
# 
# mycountsbio <- dat(mydata, factor = NULL, type = "countsbio", norm=FALSE)
# ## Count distribution per sample
# png(file=paste(PLOTSDIR, "Distribution_corrected.png", sep="/"), width=2000, height=2000)
# explo.plot(mycountsbio, toplot = "protein_coding", samples = NULL, plottype = "boxplot")
# dev.off()
# png(file=paste(PLOTSDIR, "DistributionBar_corrected.png", sep="/"), width=2000, height=2000)
# explo.plot(mycountsbio, toplot = "protein_coding", samples = NULL, plottype = "barplot")
# dev.off()
# 
# ##En 4 el FULLGC_FULLLength_TMM
# png(file=paste(PLOTSDIR, "Boxplot_corrected.png", sep="/"), width=2000, height=2000)
# boxplot(log(data3+1))
# dev.off()
# 
# ##Limpiando
# FULLGC_FULLLength_TMM<-mean10
# FULLGC_FULLLength_TMM$M<-data3
# 
# save(FULLGC_FULLLength_TMM, file=paste(RDATA, "FULLGC_FULLLength_TMM.RData", sep="/"), compress="xz")
# ##########################################################################
# ##Multidimesional PCA noise analysis
# ##########################################################################
# library("NOISeq")
# library("EDASeq")
# load(file=paste(RDATA, "FULLGC_FULLLength_TMM.RData", sep="/"))
# 
# #### 3) NOISE FILTERING
# library("ggplot2")
# library("reshape2")
# 
# pl<-ggplot(data=melt(log(FULLGC_FULLLength_TMM$M+1)),
#           aes(x=value, group=Var2, colour=Var2))+geom_density()
# pdf(file=paste(PLOTSDIR, "DensityRAWFullLog.pdf", sep="/"))
# pl
# dev.off()
# 
# myfilterRaw<-filtered.data(FULLGC_FULLLength_TMM$M, factor=FULLGC_FULLLength_TMM$Targets$Group, 
#                            norm=TRUE, cv.cutoff=100, cpm=10)
# # [1] "Filtering out low count features..."
# # [1] "15378 features are to be kept for differential expression analysis with filtering method 1" cpm=1
# # [1] "15281 features are to be kept for differential expression analysis with filtering method 1" cpm=10 Corrige y no es tan agresivo
# # [1] "15134 features are to be kept for differential expression analysis with filtering method 1" cpm=20
# # [1] "15002 features are to be kept for differential expression analysis with filtering method 1" cpm=30
# # [1] "14767 features are to be kept for differential expression analysis with filtering method 1" cpm=50
# 
# dim(FULLGC_FULLLength_TMM$M)
# # [1] 17215   881
# dim(myfilterRaw)
# # [1] 15281   881
# nrow(FULLGC_FULLLength_TMM$M)-nrow(myfilterRaw)
# # [1] 1934 ##Filtrados por baja expression cpm=10
# # [1] 2081 ##Filtrados por baja expression cpm=20
# # [1] 2213 ##Filtrados por baja expression cpm=30
# # [1] 2448 ##Filtrados por baja expression cpm=50
# 
# ##Propagando a la anotacion
# myfilterRaw<-list(M=myfilterRaw, 
#                   Annot=FULLGC_FULLLength_TMM$Annot[row.names(FULLGC_FULLLength_TMM$Annot)%in%row.names(myfilterRaw),],
#                   Targets=FULLGC_FULLLength_TMM$Targets)
# stopifnot(nrow(myfilterRaw$M)==nrow(myfilterRaw$Annot))
# stopifnot(row.names(myfilterRaw$M)==row.names(myfilterRaw$Annot))
# 
# 
# pl<-ggplot(data=melt(log(myfilterRaw$M+1)),aes(x=value, group=Var2, colour=Var2))+geom_density()
# pdf(file=paste(PLOTSDIR, "DensityFilter10.pdf", sep="/"))
# pl
# dev.off()
# 
# #### 2) PCA EXPLORATION
# setwd("/pipeline/ARSyN")
# source("sourceARSyN.R")
# setwd("/pipeline")
# 
# #pca.results = PCA.GENES(t(log2(1+data3)))
# pca.results <- PCA.GENES(t(log2(1+myfilterRaw$M)))
# traditional.pca<-prcomp(t(log2(1+myfilterRaw$M)))
# summary(traditional.pca)$importance[,1:10]
# #                             PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8      PC9     PC10
# # Standard deviation     46.35244 39.58127 26.07990 22.85858 20.19392 20.02734 17.23167 16.28167 14.59049 13.17168
# # Proportion of Variance  0.14091  0.10275  0.04461  0.03427  0.02674  0.02631  0.01947  0.01739  0.01396  0.01138
# # Cumulative Proportion   0.14091  0.24366  0.28826  0.32253  0.34928  0.37558  0.39506  0.41244  0.42640  0.43778
# 
# ## Variance explained by each component
# pdf(file=paste(PLOTSDIR, "PCAVarianceRAW.pdf", sep="/"), width = 4*2, height = 4*2)
# barplot(pca.results$var.exp[,1], xlab = "PC", ylab = "explained variance", ylim = c(0,0.4))
# dev.off()
# 
# ## Loading plot
# pdf(file=paste(PLOTSDIR, "PCALoadingRaw.pdf", sep="/"), width = 4*2, height = 4*2)
# plot(pca.results$loadings[,1:2], col = 1, pch = 20, cex = 0.5,
#      xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
#      ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
#      main = "PCA loadings",
#      xlim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1),
#      ylim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1))  
# dev.off()
# 
# ## Score plot
# mycol <- as.character(myfilterRaw$Targets$Group)
# mycol[mycol == 'N'] <- "black"
# mycol[mycol == 'T'] <- "red2"
# 
# pdf(file=paste(PLOTSDIR, "PCAScoreRaw.pdf", sep="/"), width = 5*2, height = 5)
# par(mfrow = c(1,2))
# 
# # PC1 & PC2
# rango <- diff(range(pca.results$scores[,1:2]))
# plot(pca.results$scores[,1:2], col = "white",
#      xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
#      ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
#      main = "PCA scores",
#      xlim = range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1),
#      ylim = range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1))
# points(pca.results$scores[,1], pca.results$scores[,2], col = mycol, cex = 1.5)  
# legend("topright", unique(as.character(myfilterRaw$Targets$Group)), col = rep(unique(mycol),2), ncol = 2)
# 
# # PC1 & PC3
# rango2 = diff(range(pca.results$scores[,c(1,3)]))
# plot(pca.results$scores[,c(1,3)], col = "white",
#      xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
#      ylab = paste("PC 3 ", round(pca.results$var.exp[3,1]*100,0), "%", sep = ""),
#      main = "PCA scores",
#      xlim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1),
#      ylim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1))
# points(pca.results$scores[,1], pca.results$scores[,3], col = mycol, cex = 1.5)
# legend("topright", unique(as.character(myfilterRaw$Targets$Group)), col = rep(unique(mycol),2), ncol = 2)
# 
# dev.off()
# 
# ##Probamos si ARSyN reduce el ruido########################################################
# myARSyN <- ARSyN(data=log2(myfilterRaw$M + 1), Covariates=t(as.matrix(myfilterRaw$Targets$Group)))
# 
# pca.results <- PCA.GENES(t(myARSyN))
# traditional.pca<-prcomp(t(myARSyN))
# summary(traditional.pca)$importance[,1:10]
# #                             PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8     PC9     PC10
# # Standard deviation     29.93246 5.611833 5.587665 5.541001 5.512378 5.460539 5.407716 5.380455 5.34742 5.324276
# # Proportion of Variance  0.14271 0.005020 0.004970 0.004890 0.004840 0.004750 0.004660 0.004610 0.00455 0.004520
# # Cumulative Proportion   0.14271 0.147720 0.152700 0.157590 0.162430 0.167170 0.171830 0.176440 0.18100 0.185510
# 
# ## Variance explained by each component
# pdf(file=paste(PLOTSDIR, "PCAVarianceARSyN.pdf", sep="/"), width = 4*2, height = 4*2)
# barplot(pca.results$var.exp[,1], xlab = "PC", ylab = "explained variance", ylim = c(0,0.4))
# dev.off()
# 
# ## Loading plot
# pdf(file=paste(PLOTSDIR, "PCALoadingARSyN.pdf", sep="/"), width = 4*2, height = 4*2)
# plot(pca.results$loadings[,1:2], col = 1, pch = 20, cex = 0.5,
#      xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
#      ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
#      main = "PCA loadings",
#      xlim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1),
#      ylim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1))  
# dev.off()
# 
# ## Score plot
# mycol <- as.character(myfilterRaw$Targets$Group)
# mycol[mycol == 'N'] <- "black"
# mycol[mycol == 'T'] <- "red2"
# 
# pdf(file=paste(PLOTSDIR, "PCAScoreARSyN.pdf", sep="/"), width = 5*2, height = 5)
# par(mfrow = c(1,2))
# 
# # PC1 & PC2
# rango <- diff(range(pca.results$scores[,1:2]))
# plot(pca.results$scores[,1:2], col = "white",
#      xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
#      ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
#      main = "PCA scores",
#      xlim = range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1),
#      ylim = range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1))
# points(pca.results$scores[,1], pca.results$scores[,2], col = mycol, cex = 1.5)  
# legend("topright", unique(as.character(myfilterRaw$Targets$Group)), col = rep(unique(mycol),2), ncol = 2)
# 
# # PC1 & PC3
# rango2 = diff(range(pca.results$scores[,c(1,3)]))
# plot(pca.results$scores[,c(1,3)], col = "white",
#      xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
#      ylab = paste("PC 3 ", round(pca.results$var.exp[3,1]*100,0), "%", sep = ""),
#      main = "PCA scores",
#      xlim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1),
#      ylim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1))
# points(pca.results$scores[,1], pca.results$scores[,3], col = mycol, cex = 1.5)
# legend("topright", unique(as.character(myfilterRaw$Targets$Group)), col = rep(unique(mycol),2), ncol = 2)
# 
# dev.off()
# 
# 
# ##Cómo quedaron las densidades y demás con esta modificación
# pdf(file=paste(PLOTSDIR, "BoxplotARSyN.pdf", sep="/"), width=2000, height=2000)
# boxplot(myARSyN)
# dev.off()
# 
# pl<-ggplot(data=melt(myARSyN), aes(x=value, group=Var2, colour=Var2))+geom_density()
# pdf(file=paste(PLOTSDIR, "DensityARSyN.pdf", sep="/"))
# pl
# dev.off()
# 
# dim(myARSyN)
# # [1] 15281   881
# 
# ##Limpiando
# FULLGC_FULLLength_TMM_CPM10_ARSYN<-list(M=myARSyN, Annot=myfilterRaw$Annot, Targets=myfilterRaw$Targets)
# dim(FULLGC_FULLLength_TMM_CPM10_ARSYN$M)
# # [1] 15281   881
# stopifnot(nrow(FULLGC_FULLLength_TMM_CPM10_ARSYN$M)==nrow(FULLGC_FULLLength_TMM_CPM10_ARSYN$Annot))
# stopifnot(all(row.names(FULLGC_FULLLength_TMM_CPM10_ARSYN$M)==row.names(FULLGC_FULLLength_TMM_CPM10_ARSYN$Annot)))
# 
# save(FULLGC_FULLLength_TMM_CPM10_ARSYN, file=paste(RDATA, "FULLGC_FULLLength_TMM_CPM10_ARSYN.RData", sep="/"), compress="xz")
# ##Matrices de datos para las redes
# #Todos= sanos | enfermos
# M<-as.data.frame(FULLGC_FULLLength_TMM_CPM10_ARSYN$M)
# M<-cbind(gene=as.character(FULLGC_FULLLength_TMM_CPM10_ARSYN$Annot$EnsemblId), M)
# #sanos
# sanos<-as.data.frame(FULLGC_FULLLength_TMM_CPM10_ARSYN$M[,FULLGC_FULLLength_TMM_CPM10_ARSYN$Targets$Group=="N"])
# sanos<-cbind(gene=as.character(FULLGC_FULLLength_TMM_CPM10_ARSYN$Annot$EnsemblId), sanos)
# #enfermos
# enfermos<-as.data.frame(FULLGC_FULLLength_TMM_CPM10_ARSYN$M[,FULLGC_FULLLength_TMM_CPM10_ARSYN$Targets$Group=="T"])
# enfermos<-cbind(gene=as.character(FULLGC_FULLLength_TMM_CPM10_ARSYN$Annot$EnsemblId), enfermos)
# #Symbolos
# symbol<-as.character(FULLGC_FULLLength_TMM_CPM10_ARSYN$Annot$EnsemblId)
# ##Guardando
# write.table(M, file = paste(RDATA, "M.tab", sep="/"), sep="\t", quote=FALSE, row.names=FALSE)
# write.table(sanos, file = paste(RDATA, "Sanos.txt", sep="/"), sep="\t", quote=FALSE, row.names=FALSE)
# write.table(enfermos, file = paste(RDATA, "Enfermos.txt", sep="/"), sep="\t", quote=FALSE, row.names=FALSE)
# write.table(symbol, file = paste(RDATA, "ListaGenes.txt", sep="/"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
# save(FULLGC_FULLLength_TMM_CPM10_ARSYN, file=paste(RDATA, "FULLGC_FULLLength_TMM_CPM10_ARSYN.RData", sep="/"), compress="xz")

