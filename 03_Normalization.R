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

## This function gets the 
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
      cat("Testing with length normalization: ", ln, ", GC normalization: ", gcn, " and between lanes normalization: ", bn, "\n")
      if (bn == "tmm") {
        between.data <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
      } else {
        between.data <- betweenLaneNormalization(gcn.data, which = bn, offset = FALSE)
      }
  
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
      cat("Testing with GC normalization: ", gcn, ",  length normalization: ", ln, " and between lane normalization: ", bn, "\n")
      if (bn == "tmm") {
        between.data <- tmm(ln.data, long = 1000, lc = 0, k = 0)
      } else {
        between.data <- betweenLaneNormalization(ln.data, which = bn, offset = FALSE)
      }
      
      norm.noiseq.results <- getNOISeqResults(paste("GC", gcn, sep = "."), paste("Length", ln, sep = "."), paste("Between", bn, sep =  "."), 
                                              between.data, mean10)
      normalization.results <- rbind(normalization.results, norm.noiseq.results)                      
    }
  }
}

write.table(normalization.results, file=paste(RDATA, "NormalizationResults.tsv", sep="/"), 
            quote = F, sep = "\t", row.names = F)
