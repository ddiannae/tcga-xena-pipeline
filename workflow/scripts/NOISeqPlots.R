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