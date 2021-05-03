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
##      Within lane. For lenght: RPKM, loess, full
##      Within lane. For GC content: loess, full
##      Between lanes. TMM, full
##  - Integrate plots per normalization type in one png
##  - Integrate all results in one pdf
###############################################################################
library(ggplot2)
library(reshape2)
library(grid)
library(png)
library(gridExtra)
library(parallel)
library(readr)
library(dplyr)

###############################################################################
args <- commandArgs(trailingOnly = T)

if (length(args) < 2 ) {
  stop("Incorrect number of arguments", call.=FALSE)
} else {
  TISSUE = args[1]
  DATADIR = args[2]
  MCCORES = args[3]
}
DATADIR <- paste(DATADIR, TISSUE, sep="/")
PLOTSDIR <-paste(DATADIR, "plots", sep="/")
PLOTSNORMDIR <- paste(PLOTSDIR, "normalization", sep = "/")
w <- 1024
h <- 1024
p <- 24

{ 
  cat("Integrating plots \n.")
  
  lenght_norm <- c("full", "loess", "median", "upper")
  gc_norm <- c( "full", "loess", "median", "upper")
  between_nom <- c("full", "median", "tmm", "upper")
  
  ## This function will retrieve plots for one normalization set and will create 
  ## one single plot
  savePlots <- function(step1, step2, step3) {
   
    plotname <- paste(step1, step2, step3, sep = "_")  
    pngPlots <- c(paste0(PLOTSNORMDIR, "/", plotname, "_lenght_bias.png"), 
                  paste0(PLOTSNORMDIR, "/", plotname, "_gc_bias.png"), 
                  paste0(PLOTSNORMDIR, "/", plotname, "_rna_composition.png"))
    
    thePlots <- lapply (pngPlots, function(pngFile) {
      rasterGrob(readPNG(pngFile, native = FALSE), interpolate = FALSE)
    })
    
    png(paste(PLOTSNORMDIR, paste(plotname, "png", sep ="."), sep="/"),  height=h/2, width=w*(3/2))
    par(oma = c(0, 0, 1.5, 0))
    plot.new()
    do.call(grid.arrange, c(thePlots,  ncol = 3, nrow=1))
    mtext(plotname, outer = TRUE, cex = 1.5)
    dev.off()

    unlink(pngPlots)
  } 
  
  df_normalizations <- expand.grid(gcn = gc_norm, ln = lenght_norm, 
                                   bn = between_nom, stringsAsFactors = F)
  
  cat("Getting plots for all ", nrow(df_normalizations), "normalization combinations.\n")

  plots <- mclapply(X = 1:nrow(df_normalizations), mc.cores = MCCORES, 
                       FUN = function(i){
                         gcn <- df_normalizations[i, "gcn"]
                         ln <- df_normalizations[i, "ln"]
                         bn <- df_normalizations[i, "bn"]
                         
                         savePlots(paste("gc", gcn, sep = "_"), 
                                   paste("length", ln, sep = "_"), 
                                   paste("between", bn, sep =  "_"))
                         
                         savePlots(paste("length", ln, sep = "_"), 
                                   paste("gc", gcn, sep = "_"), 
                                   paste("between", bn, sep =  "_"))
                       })
  
  ## Finally, we test with cqn
  cat("Getting plots for cqn normalization\n")
  savePlots("gc_cqn", "length_cqn", "between_cqn")
  
  pngPlots <- list.files(path = PLOTSNORMDIR, pattern = "*.png", full.names = TRUE)

  thePlots <- lapply (pngPlots, function(pngFile) {
    rasterGrob(readPNG(pngFile, native = FALSE), interpolate = FALSE)
  })
  
  cat("Saving normalization_plots.pdf\n") 
  pdf(paste(PLOTSDIR, "normalization_plots.pdf", sep="/"), title="Normalization Plots")
  print(marrangeGrob(thePlots, nrow = 3, ncol = 1, top = NULL))
  dev.off()
}
