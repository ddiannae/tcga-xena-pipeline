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
rownames(mean10$Annot) <- rownames(mean10$M)
save(mean10, file=paste(RDATA, "Mean10.RData", sep="/"), compress="xz")
cat("Saving RawFull.RData \n")
##########################################################################
#Normalization
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
  mylengthbias <- dat(mydata, factor="Group", norm = TRUE, type="lengthbias")
  regmodels <- mylengthbias@dat$RegressionModels
  name1 <- names(regmodels[1])
  rsquared.1 <- summary(regmodels[[1]])$r.squared
  fstatistic.1 <- summary(regmodels[[1]])$fstatistic
  pvalue.1 <- signif(pf(q = fstatistic.1[1], df1 = fstatistic.1[2], df2 = fstatistic.1[3], lower.tail = FALSE), 2)
  name2 <- names(regmodels[2])
  rsquared.2 <- summary(regmodels[[2]])$r.squared
  fstatistic.2 <- summary(regmodels[[2]])$fstatistic
  pvalue.2 <- signif(pf(q = fstatistic.2[1], df1 = fstatistic.2[2], df2 = fstatistic.2[3], lower.tail = FALSE), 2)
  
  png(paste(PLOTSNORMDIR, paste(step1, step2, step3, "Lenghtbias.png", sep = "_"), sep="/"), width=w, height=h, pointsize=p)
  explo.plot(mylengthbias, samples = NULL, toplot = "global")
  dev.off()
  
  ## GC Bias
  mygcbias <- dat(mydata, factor = "Group", norm = TRUE, type ="GCbias")
  regmodels.gc <- mylengthbias@dat$RegressionModels
  name1.gc <- names(regmodels.gc[1])
  rsquared.1.gc <- summary(regmodels.gc[[1]])$r.squared
  fstatistic.1.gc <- summary(regmodels.gc[[1]])$fstatistic
  pvalue.1.gc <- signif(pf(q = fstatistic.1.gc[1], df1 = fstatistic.1.gc[2], df2 = fstatistic.1.gc[3], lower.tail = FALSE), 2)
  name2.gc <- names(regmodels.gc[2])
  rsquared.2.gc <- summary(regmodels.gc[[2]])$r.squared
  fstatistic.2.gc <- summary(regmodels.gc[[2]])$fstatistic
  pvalue.2.gc <- signif(pf(q = fstatistic.2.gc[1], df1 = fstatistic.2.gc[2], df2 = fstatistic.2.gc[3], lower.tail = FALSE), 2)
  
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
                                 rsquared.1, pvalue.1, rsquared.2, pvalue.2,
                                 rsquared.1.gc, pvalue.1.gc, rsquared.2.gc, pvalue.2.gc, 
                                 dtable["PASSED"], dtable["PASSED"]/nsamples)
  colnames(norm.set.results) <- c("Step1", "Step2", "Step3", 
                                  paste("Lenght", name1, "R2", sep = "."), paste("Lenght", name1, "p-value", sep = "."),  
                                  paste("Lenght", name2, "R2", sep = "."), paste("Lenght", name2, "p-value", sep = "."), 
                                  paste("GC", name1.gc, "R2", sep = "."), paste("GC", name1.gc, "p-value", sep = "."), 
                                  paste("GC", name2.gc, "R2", sep = "."), paste("GC", name2.gc, "p-value", sep = "."),
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
      cat("Testing with length normalization: ", ln, ", GC normalization: ", gcn, " and between lanes normalization: ", bn, "\n")
      if (bn == "tmm") {
        normalized.data <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
      } else {
        normalized.data <- betweenLaneNormalization(gcn.data, which = bn, offset = FALSE)
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
        normalized.data <- betweenLaneNormalization(ln.data, which = bn, offset = FALSE)
      }
      
      norm.noiseq.results <- getNOISeqResults(paste("GC", gcn, sep = "."), paste("Length", ln, sep = "."), paste("Between", bn, sep =  "."), 
                                              normalized.data, mean10)
      normalization.results <- rbind(normalization.results, norm.noiseq.results)                      
    }
  }
}

write.table(normalization.results, file=paste(RDATA, "NormalizationResults.tsv", sep="/"), 
            quote = F, sep = "\t", row.names = F)
 
