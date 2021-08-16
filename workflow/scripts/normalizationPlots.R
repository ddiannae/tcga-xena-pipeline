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
##  Integrate normalization plots. No done with the normalizationTest script 
##  because it causes problems with parallelization
##  - Normalization combinations:
##      Within lane. For length: RPKM, loess, full
##      Within lane. For GC content: loess, full
##      Between lanes. TMM, full
##  - Integrate plots per normalization type in one png
##  - Integrate all results in one pdf
###############################################################################
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(ggplot2)
library(reshape2)
library(grid)
library(png)
library(gridExtra)
library(readr)
library(dplyr)

###############################################################################
args <- commandArgs(trailingOnly = T)

PLOTSDIR <-paste(snakemake@params[["tissue_dir"]], "plots", sep="/")
PLOTSNORMDIR <- paste(PLOTSDIR, "normalization", sep = "/")
w <- 1024
h <- 1024
p <- 24

{ 
  cat("Integrating plots \n.")
  
  length_norm <- c("no","full", "loess", "median", "upper")
  gc_norm <- c("no", "full", "loess", "median", "upper")
  between_norm <- c("no", "full", "median", "tmm", "upper")
  
  ## This function will retrieve plots for one normalization set and will create 
  ## one single plot
  savePlots <- function(step1, step2, step3) {
   
    plotname <- paste(step1, step2, step3, sep = "_")  
    pngPlots <- c(paste0(PLOTSNORMDIR, "/", plotname, "_length_bias.png"), 
                  paste0(PLOTSNORMDIR, "/", plotname, "_gc_bias.png"), 
                  paste0(PLOTSNORMDIR, "/", plotname, "_rna_composition.png"))
    
    thePlots <- lapply (pngPlots, function(pngFile) {
      rasterGrob(readPNG(pngFile, native = FALSE), interpolate = FALSE)
    })
    
    plotpath <- paste(PLOTSNORMDIR, paste(plotname, "png", sep ="."), sep="/")
    png(plotpath,  height=h/2, width=w*(3/2))
    par(oma = c(0, 0, 1.5, 0))
    plot.new()
    do.call(grid.arrange, c(thePlots,  ncol = 3, nrow=1))
    mtext(plotname, outer = TRUE, cex = 1.5)
    dev.off()

    unlink(pngPlots)
    return(plotpath)
  } 
  
  df_normalizations <- expand.grid(gcn = gc_norm, ln = length_norm, 
                                   bn = between_norm, stringsAsFactors = F)
  
  cat("Getting plots for all ", nrow(df_normalizations), "normalization combinations.\n")

  plots <- lapply(1:nrow(df_normalizations),  function(i){
                         gcn <- df_normalizations[i, "gcn"]
                         ln <- df_normalizations[i, "ln"]
                         bn <- df_normalizations[i, "bn"]
                         gcp <- NA
                         tryCatch(expr = {
                           gcp <- savePlots(paste("gc", gcn, sep = "_"), 
                                     paste("length", ln, sep = "_"), 
                                     paste("between", bn, sep =  "_"))
                         }, error = function(cond){
                          
                           cat(cond$message, "\n")
                         })
                         lp <- NA
                         tryCatch(expr = {
                           lp <- savePlots(paste("length", ln, sep = "_"), 
                                     paste("gc", gcn, sep = "_"), 
                                     paste("between", bn, sep =  "_"))
                         }, error = function(cond){
                           cat(cond$message, "\n")
                         })
              return(list(gcp, lp))         
          })
  
  plots <- sort(unlist(plots))
  
  thePlots <- lapply (plots, function(pngFile) {
    rasterGrob(readPNG(pngFile, native = FALSE), interpolate = FALSE)
  })
  
  cat("Saving normalization_plots.pdf\n") 
  pdf(paste(PLOTSDIR, "normalization_plots.pdf", sep="/"), title="Normalization Plots")
  print(marrangeGrob(thePlots, nrow = 3, ncol = 1, top = NULL))
  dev.off()
}
