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
##  - Try Normalization combinations:
##      Within lane. For length: RPKM, loess, full
##      Within lane. For GC content: loess, full
##      Between lanes. TMM, full
##  - Generate plots
###############################################################################
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(ggplot2)
library(reshape2)
library(grid)
library(png)
library(gridExtra)
library(parallel)
library(readr)
library(dplyr)
library(NOISeq)
library(EDASeq)
library(DESeq2)

###############################################################################
MCCORES <- as.numeric(snakemake@params[["mccores"]])
PLOTSDIR <-paste(snakemake@params[["tissue_dir"]], "plots", sep="/")
PLOTSNORMDIR <- paste(PLOTSDIR, "normalization", sep = "/")
dir.create(PLOTSNORMDIR)
w <- 1024
h <- 1024
p <- 24

{
  load(snakemake@input[[1]])
  
  cat("Testing normalization methods\n.")
  mydataM10EDA <- EDASeq::newSeqExpressionSet(
    counts = mean10$M,
    featureData = mean10$annot %>% as.data.frame(),
    phenoData = data.frame(
      conditions = mean10$targets$group,
      row.names = mean10$targets$id))
  
  length_norm <- c("no", "full", "loess", "median", "upper")
  gc_norm <- c("no", "full", "loess", "median", "upper")
  between_norm <-c("full", "median", "tmm", "upper")
  
  ## This function gets the relevant statistics for the regression methods for GC and Length bias
  getRegressionStatistics <- function(regressionmodel) {
    name <- names(regressionmodel)
    rsquared <- summary(regressionmodel[[1]])$r.squared
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
      length = m10_data$annot %>% select(gene_id, length) %>% as.data.frame(), 
      biotype = m10_data$annot %>% select(gene_id, gene_type) %>% as.data.frame(), 
      chromosome = m10_data$annot %>% select(chr, start, end) %>% as.data.frame(),  
      factors = m10_data$targets %>% select(group) %>% as.data.frame(), 
      gc = m10_data$annot %>% select(gene_id, gc)%>% as.data.frame())

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
  
    cat("Step1:", step1, ", step2:", step2, " step3: ", step3, " calculations done\n")
    plotname <- paste(step1, step2, step3, sep = "_")  
    pngPlots <- c(paste0(PLOTSNORMDIR, "/", plotname, "_length_bias.png"), 
                  paste0(PLOTSNORMDIR, "/", plotname, "_gc_bias.png"), 
                  paste0(PLOTSNORMDIR, "/", plotname, "_rna_composition.png"))
    
    png(pngPlots[1], width=w/2, height=h/2)
    explo.plot(mylengthbias, samples = NULL, toplot = "global")
    dev.off()
    
    png(pngPlots[2], width=w/2, height=h/2)
    explo.plot(mygcbias, samples = NULL, toplot = "global")
    dev.off()
    
    png(pngPlots[3],width=w/2, height=h/2)
    explo.plot(myrnacomp, samples = 1:12)
    dev.off()
    
    cat("Step1:", step1, ", step2:", step2, " step3: ", step3, " plots saved\n")
    
    norm_set_results <- data.frame(step1, step2, step3, 
                                   l_stats_1$r2, l_stats_1$p, l_stats_2$r2, l_stats_2$p,
                                   gc_stats_1$r2, gc_stats_1$p, gc_stats_2$r2, gc_stats_2$p, 
                                   dtable["PASSED"], dtable["PASSED"]/nsamples)
    
    colnames(norm_set_results) <- c("step1", "step2", "step3", 
                                    paste("length", l_stats_1$name, "R2", sep = "_"), paste("length", l_stats_1$name, "p-value", sep = "_"),  
                                    paste("length", l_stats_2$name, "R2", sep = "_"), paste("length", l_stats_2$name, "p-value", sep = "_"), 
                                    paste("gc", gc_stats_1$name, "R2", sep = "_"), paste("gc", gc_stats_1$name, "p-value", sep = "_"), 
                                    paste("gc", gc_stats_2$name, "R2", sep = "_"), paste("gc", gc_stats_2$name, "p-value", sep = "_"),
                                    "rna_passed_samples", "rna_passed_proportion")
    return(norm_set_results)
  } 
  
  df_normalizations <- expand.grid(gcn = gc_norm, ln = length_norm, 
                                  bn = between_norm, stringsAsFactors = F)
  
  cat("Testing all ", nrow(df_normalizations), "normalization combinations\n")
  
  ## Normalizations with GC step first
  gc_norms <- mclapply(X = 1:nrow(df_normalizations), 
                        mc.cores = MCCORES, 
                        FUN = function(i){
                          
    gcn <- df_normalizations[i, "gcn"]
    ln <- df_normalizations[i, "ln"]
    bn <- df_normalizations[i, "bn"]
    
    if(gcn == "no") {
      gcn_data <- counts(mydataM10EDA)
    } else {
      gcn_data <- withinLaneNormalization(counts(mydataM10EDA), mean10$annot$gc, which = gcn)  
    }
    
    if(ln == "no") {
      ln_data <- gcn_data
    } else {
      ln_data <- withinLaneNormalization(gcn_data, mean10$annot$length, which = ln)
    }
        
    if (bn == "tmm") {
       between_data <- tmm(ln_data, long = 1000, lc = 0, k = 0)
    } else {
       between_data <- betweenLaneNormalization(ln_data, which = bn, offset = FALSE)
    }
    cat("Testing with GC normalization: ", gcn, ", length normalization: ", ln, " and between lane normalization: ", bn, "\n")
    norm_noiseq_results <- getNOISeqResults(paste("gc", gcn, sep = "_"), paste("length", ln, sep = "_"), paste("between", bn, sep =  "_"), 
                                                       between_data, mean10)
    return(norm_noiseq_results)
  })
  
  gc_norms <- bind_rows(gc_norms)
  
  ## Normalizations with Length step first
  ln_norms <- mclapply(X = 1:nrow(df_normalizations), 
                        mc.cores = MCCORES,
                        FUN = function(i){
                          
          gcn <- df_normalizations[i, "gcn"]
          ln <- df_normalizations[i, "ln"]
          bn <- df_normalizations[i, "bn"]
          
          if(ln == "no") {
            ln_data <- counts(mydataM10EDA)
          } else {
            ln_data <- withinLaneNormalization(counts(mydataM10EDA), mean10$annot$length, which = ln)
          }
          
          if(gcn == "no") {
            gcn_data <- ln_data
          } else {
            gcn_data <- withinLaneNormalization(ln_data, mean10$annot$gc, which = gcn)
          }
          
          if (bn == "tmm") {
            between_data <- tmm(ln_data, long = 1000, lc = 0, k = 0)
          } else {
            between_data <- betweenLaneNormalization(ln_data, which = bn, offset = FALSE)
          }
          cat("Testing with GC normalization: ", gcn, ", length normalization: ", ln, " and between lane normalization: ", bn, "\n")
          norm_noiseq_results <- getNOISeqResults(paste("length", ln, sep = "_"), paste("gc", gcn, sep = "_"), paste("between", bn, sep =  "_"), 
                                                  between_data, mean10)
          return(norm_noiseq_results)                                         
  })
  
  ln_norms <- bind_rows(ln_norms)
  
  normalization_results <- bind_rows(gc_norms, ln_norms)
 
  cat("End of normalization testing\n")
  cat("Saving normalization_results.tsv\n") 
  write_tsv(normalization_results, snakemake@output[[1]])
}
