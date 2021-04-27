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
##  - Try Normalization combinations:
##      Within lane. For lenght: RPKM, loess, full
##      Within lane. For GC content: loess, full
##      Between lanes. TMM, full
##  - Integrate results in one pdf
###############################################################################
library(ggplot2)
library(reshape2)
library(grid)
library(png)
library(gridExtra)
library(readr)
library(dplyr)
library(NOISeq)
library(EDASeq)
library(DESeq2)
library(cqn)

###############################################################################
args <- commandArgs(trailingOnly = T)

if (length(args) < 2 ) {
  stop("Incorrect number of arguments", call.=FALSE)
} else {
  TISSUE = args[1]
  DATADIR = args[2]
}
DATADIR <- paste(DATADIR, TISSUE, sep="/")
RDATA <- paste(DATADIR, "rdata", sep="/")
PLOTSDIR <-paste(DATADIR, "plots", sep="/")
PLOTSNORMDIR <- paste(PLOTSDIR, "normalization", sep = "/")
dir.create(PLOTSNORMDIR)
w <- 1024
h <- 1024
p <- 24

##########################################################################
load(file=paste(RDATA, "raw_full.RData", sep="/"))
{ 
  ### We keep only genes with mean expression count > 10 
  exp_genes <- apply(full$M, 1, function(x) mean(x)>10)
  egtable <- table(exp_genes)
  cat("There are", egtable[[1]], "genes with mean expression count < 10", egtable[[2]], "with mean count > 10 \n")
  
  ##Filtering low expression genes
  mean10 <- list(M = full$M[exp_genes, ], annot=full$annot[exp_genes, ], targets=full$targets)
  rownames(mean10$annot) <- rownames(mean10$M)
  
  cat("Filtering protein coding features \n")
  protein_coding <- rownames(mean10$annot[mean10$annot$gene_type == "protein_coding", ])
  cat("There are ", length(protein_coding), " features with type: protein_coding \n")
  mean10 <- list(M = mean10$M[protein_coding, ], annot=mean10$annot[protein_coding, ], targets=mean10$targets)
  
  cat("Saving mean10_ProteinCoding.RData \n") 
  save(mean10, file=paste(RDATA, "mean10_proteinCoding.RData", sep="/"), compress="xz")
}
 ###########################################################################
{ 
  cat("Testing normalization methods\n.")
  mydataM10EDA <- EDASeq::newSeqExpressionSet(
    counts=mean10$M,
    featureData=mean10$annot,
    phenoData=data.frame(
      conditions=mean10$targets$group,
      row.names=colnames(mean10$M)))
  
  lenght_norm <- c("full")
  gc_norm <- c( "full" )
  between_nom <- c("full")
  
  #lenght_norm <- c("full", "loess", "median", "upper")
  #gc_norm <- c( "full", "loess", "median", "upper")
  #between_nom <- c("full", "median", "tmm", "upper")
  #normalization_results <- data.frame()
  
  ## This function gets the relevant statistics for the regression methods for GC and Length bias
  getRegressionStatistics <- function(regressionmodel) {
    name <- names(regressionmodel)
    rsquared <- summary(regressionmodel[[1]])$r.squared
    print(rsquared)
    fstatistic <- summary(regressionmodel[[1]])$fstatistic
    pvalue <- signif(pf(q = fstatistic[1], df1 = fstatistic[2], 
                        df2 = fstatistic[3], lower.tail = FALSE), 2)
    return(list("name" = name, "r2" = rsquared, "p" = pvalue))
  }
  
  ## This function will retrieve the NOISeq results to test the normalization combination
  getNOISeqResults <- function(step1, step2, step3, n_counts, m10_data) {
    ### Check the NOISEq results 
    mydata <- NOISeq::readData(
      data = n_counts, 
      length = m10_data$annot %>% select(gene_id, width), 
      biotype = m10_data$annot %>% select(gene_id, gene_type), 
      chromosome = m10_data$annot %>% select(chr, start, end), 
      factors = m10_data$target %>% select(group),
      gc = m10_data$annot %>% select(gene_id, gc))

    nsamples <- dim(m10_data$targets)[1]
    
    ### Length bias 
    mylengthbias <- dat(mydata, factor="group", norm = TRUE, type="lengthbias")
    l_stats_1 <- getRegressionStatistics(mylengthbias@dat$RegressionModels[1])
    l_stats_2 <- getRegressionStatistics(mylengthbias@dat$RegressionModels[2])

    ## GC Bias
    mygcbias <- dat(mydata, factor = "group", norm = TRUE, type ="GCbias")
    gc_stats_1 <- getRegressionStatistics(mygcbias@dat$RegressionModels[1])
    gc_stats_2 <- getRegressionStatistics(mygcbias@dat$RegressionModels[2])
   
    #RNA Composition
    myrnacomp <- dat(mydata, norm = TRUE, type="cd")
    dtable <- table(myrnacomp@dat$DiagnosticTest[,  "Diagnostic Test"])
    if (is.na(dtable["PASSED"])) dtable <- data.frame(PASSED = 0)
    
    pngPlots <- c(paste(PLOTSNORMDIR, "lenght_bias.png", sep="/"), 
                  paste(PLOTSNORMDIR, "gc_bias.png", sep="/"), 
                  paste(PLOTSNORMDIR, "rna_composition.png", sep="/"))
    
    png(pngPlots[1], width=w/2, height=h/2)
    explo.plot(mylengthbias, samples = NULL, toplot = "global")
    dev.off()
    
    png(pngPlots[2], width=w/2, height=h/2)
    explo.plot(mygcbias, samples = NULL, toplot = "global")
    dev.off()
    
    png(pngPlots[3],width=w/2, height=h/2)
    explo.plot(myrnacomp, samples = 1:12)
    dev.off()
    
    thePlots <- lapply (pngPlots, function(pngFile) {
      rasterGrob(readPNG(pngFile, native = FALSE), interpolate = FALSE)
    })
    
    plotname <- paste(step1, step2, step3, sep = "_")
    png(paste(PLOTSNORMDIR, paste(plotname, "png", sep ="."), sep="/"),  height=h/2, width=w*(3/2))
    par(oma = c(0, 0, 1.5, 0))
    plot.new()
    do.call(grid.arrange, c(thePlots,  ncol = 3, nrow=1))
    mtext(plotname, outer = TRUE, cex = 1.5)
    dev.off()
    
    unlink(pngPlots)
    
    norm_set_results <- data.frame(step1, step2, step3, 
                                   l_stats_1$r2, l_stats_1$p, l_stats_2$r2, l_stats_2$p,
                                   gc_stats_1$r2, gc_stats_1$p, gc_stats_2$r2, gc_stats_2$p, 
                                   dtable["PASSED"], dtable["PASSED"]/nsamples)
    colnames(norm_set_results) <- c("step1", "step2", "step3", 
                                    paste("lenght", l_stats_1$name, "R2", sep = "_"), paste("lenght", l_stats_1$name, "p-value", sep = "_"),  
                                    paste("lenght", l_stats_2$name, "R2", sep = "_"), paste("lenght", l_stats_2$name, "p-value", sep = "_"), 
                                    paste("gc", gc_stats_1$name, "R2", sep = "_"), paste("gc", gc_stats_1$name, "p-value", sep = "_"), 
                                    paste("gc", gc_stats_2$name, "R2", sep = "_"), paste("gc", gc_stats_2$name, "p-value", sep = "_"),
                                    "rna_passed_samples", "rna_passed_proportion")
    return(norm_set_results)
  } 
  df_normalizations <- expand.grid(gcn = gc_norm, ln = lenght_norm, bn = between_nom)
  ## We try with GC normalization first
  for (gcn in gc_norm) {
    gcn_data <- withinLaneNormalization(counts(mydataM10EDA), mean10$annot$gc, which = gcn)
    for (ln in lenght_norm) {
      ln_data <- withinLaneNormalization(gcn_data, mean10$annot$width, which = ln)
      for (bn in between_nom) {
        if (bn == "tmm") {
          between_data <- tmm(ln_data, long = 1000, lc = 0, k = 0)
        } else {
          between_data <- betweenLaneNormalization(ln_data, which = bn, offset = FALSE)
        }
        cat("Testing with GC normalization: ", gcn, ", length normalization: ", ln, " and between lane normalization: ", bn, "\n")
        norm_noiseq_results <- getNOISeqResults(paste("gc", gcn, sep = "_"), paste("length", ln, sep = "_"), paste("between", bn, sep =  "_"), 
                                                between_data, mean10)
        normalization_results <- bind_rows(normalization_results, norm_noiseq_results)                      
      }
    }
  }
  
  ## We try with length normalization now
  for (ln in lenght_norm) {
    ln_data <- withinLaneNormalization(counts(mydataM10EDA), mean10$annot$width, which = ln)
    for (gcn in gc_norm) {
      gcn_data <- withinLaneNormalization(ln_data, mean10$annot$gc, which = gcn)
      for (bn in between_nom) {
        if (bn == "tmm") {
          between_data <- tmm(gcn_data, long = 1000, lc = 0, k = 0)
        } else {
          between_data <- betweenLaneNormalization(gcn_data, which = bn, offset = FALSE)
        }
        cat("Testing with length normalization: ", ln, ", GC normalization: ", gcn, " and between lanes normalization: ", bn, "\n")
        norm_noiseq_results <- getNOISeqResults(paste("length", ln, sep = "_"), paste("GC", gcn, sep = "_"), paste("between", bn, sep =  "_"),
                                                between_data, mean10)
        normalization_results <- bind_rows(normalization_results, norm_noiseq_results)                      
      }
    }
  }

  ## Finally, we test with cqn
  ##Length, GC content and size correction
  cat("Testing with cqn normalization\n")
  y_DESeq <- DESeqDataSetFromMatrix(countData=mean10$M, 
                                  colData=mean10$targets, design= ~group)
  y_DESeq <- estimateSizeFactors(y_DESeq)
  
  cqn_mean10<- cqn(mean10$M, lengths=mean10$annot$width,
                   x = mean10$annot$gc, sizeFactors=sizeFactors(y_DESeq), verbose=TRUE)
  normalized_cqn <- cqn_mean10$y + cqn_mean10$offset
  norm_noiseq_results <- getNOISeqResults("gc_cqn", "length_cqn", "between_cqn", normalized_cqn, mean10)
  
  normalization_results <- bind_rows(normalization_results, norm_noiseq_results)
  pngPlots <- list.files(path = PLOTSNORMDIR, pattern = "*.png", full.names = TRUE)
  
  thePlots <- lapply (pngPlots, function(pngFile) {
    rasterGrob(readPNG(pngFile, native = FALSE), interpolate = FALSE)
  })
  
  pdf(paste(PLOTSDIR, "normalization_plots.pdf", sep="/"), title="Normalization Plots")
  print(marrangeGrob(thePlots, nrow = 3, ncol = 1, top = NULL))
  dev.off()
  
  cat("End of normalization texting\n")
  cat("Saving NormalizationResults.tsv\n") 
  write_tsv(normalization_results, paste(RDATA, "normalization_results.tsv", sep="/"))
}