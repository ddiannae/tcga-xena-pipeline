##############################################################################
## CHROMATIN AND GENE REGULATION: FROM GENE TO GENOME FOLDING
## Practical session: In silico analysis of RNA-Seq data.
###############################################################################
## Data preparation: 
##      -Quality Control & bias removal
## Author: 
##          Dr. Cristobal Fresno - cristobalfresno@gmail.com
## Date: 2016-12-12
###############################################################################
## Quality Control 
##      -Let's keep only the GC & length annotated genes
##      -EXPLORATORY ANALYSIS (NOISeq package)
##          -Reading data into NOISeq package -> mydata
##          -Plots
##              -Biodetection plot
##              -Count distribution per biotype
##              -Saturation plot
##              -Count distribution per sample
##              -Count distribution per Experimental factors
##          -Bias
##              -Length bias detection
##              -GC bias
##              -RNA composition
##          -Quality Control Report 
##     -Basal situation
##          GC content bias: Detected
##         Gene Length bias: Detected
##         RNA content bias: Detected
##     -CQN Correction
##          -Filter genes with mean < 10
##          -Length, GC content and size correction
##          -Normalized expression log2 scale
##          -Check the correction with NOISEQ
##                  -Length bias?
##                  -GC bias?
##                  -RNA composition?
##                  -Count distribution?
##                  -Cleaning up
##      -Multidimesional PCA noise analysis
##          -Distribution exploration
##              -Expression boxplot
##              -Removing low expression
##              -Annotation propagation
##          -PCA EXPLORATION
##              -Variance explained by each component
##              -Loading plot
##              -Score plot
###############################################################################
##Usefull Libraries
###############################################################################
library("BiocParallel")
library("parallel")
library("NOISeq")
library("cqn")
library("DESeq2")
library("ggplot2")
library("reshape2")
# register(SnowParam(workers=detectCores()-1, progress=TRUE))#Windows
register(MulticoreParam(workers=detectCores()-1, progress=TRUE))#Linux
options(width=80)
###############################################################################
##Quality Control
###############################################################################
DATADIR <- '/pipeline/data/'
args <- commandArgs(trailingOnly = TRUE)
DATADIR <- args[1]
RDATA <- paste(DATADIR, "rdata", sep="")
PLOTSDIR <-paste(DATADIR, "plots", sep="")
dir.create(PLOTSDIR)
w <- 1024
h <- 1024
p <- 24

load(file=paste(RDATA, "RawFull.RData", sep="/"))
{###Let's keep only the GC & length annotated genes
ids<-!is.na(full$Annot$GC) & !is.na(full$Annot$Length)

full$M<-full$M[ids,]
full$Annot<-full$Annot[ids,]
cat("Non GC and lenght annotated genes removed\n")

}##########################################
## EXPLORATORY ANALYSIS (NOISeq package)
##########################################
{## Reading data into NOISeq package -> mydata
row.names(full$M)<-full$Annot$EnsemblID
row.names(full$Annot)<-full$AnnotdEnsemblID
row.names(full$Targets)<-full$Targets$ID
full$Targets$Group<-factor(substr(full$Targets$ID, start=1, stop=1))
nsamples = table(full$Targets$Group)[1]
tsamples = table(full$Targets$Group)[2]

mydata <- NOISeq::readData(
  data = full$M, 
  length = full$Annot[, c("EnsemblID", "Length")], 
  biotype = full$Annot[, c("EnsemblID", "Type")], 
  chromosome = full$Annot[, c("Chr", "Start", "End")], 
  factors = full$Targets[, "Group",drop=FALSE], 
  gc = full$Annot[, c("EnsemblID", "GC")])

}##########################################
{## Plots
##########################################
# Biodetection plot
mybiodetection <- dat(mydata, type="biodetection", factor="Group", k=0)
png(filename=paste(PLOTSDIR, "biodetection.Rd_%03d.png", sep="/"),  width=w, height=h, pointsize=p)
explo.plot(mybiodetection)
dev.off()
cat("Biodetection plots generated\n")
#What do we need to see here?

## Count distribution per biotype
mycountsbio <- dat(mydata, factor = NULL, type = "countsbio")
png(filename=paste(PLOTSDIR, "countsbio.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mycountsbio, toplot = 1, samples = 1, plottype = "boxplot")
dev.off()
cat("Counts distribution plot per biotype generated\n")
#What about expression level?

## Count distribution per sample
mycountsbio <- dat(mydata, factor = NULL, type = "countsbio")
png(paste(PLOTSDIR, "protein_coding_boxplot.png", sep="/"), width=w*2, height=h, pointsize=p)
explo.plot(mycountsbio, toplot = "protein_coding",
    samples = NULL, plottype = "boxplot")
dev.off()
cat("Counts distribution plot per sample generated\n")

png(paste(PLOTSDIR, "protein_coding_barplot.png", sep="/"), width=w*2, height=h, pointsize=p)
explo.plot(mycountsbio, toplot = "protein_coding",
    samples = NULL, plottype = "barplot")
dev.off()
cat("Counts distribution barplot for protein coding biotype generated\n")

mycountsbio <- dat(mydata, factor = "Group", type = "countsbio")
## Count distribution per Experimental factors
png(paste(PLOTSDIR, "protein_coding_boxplot_group.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mycountsbio, toplot = "protein_coding",
    samples = NULL, plottype = "boxplot")
dev.off()
cat("Counts distribution boxplot for protein coding biotype and group generated\n")

png(paste(PLOTSDIR, "protein_coding_barplot_group.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mycountsbio, toplot = "protein_coding",
    samples = NULL, plottype = "barplot")
dev.off()
cat("Counts distribution barplot for protein coding biotype and group generated\n")
#How much sensitivity we loose? 

## Saturation plot
#mysaturation <- dat(mydata, k = 0, ndepth = 7, type = "saturation")
#png(paste(PLOTSDIR, "saturation.png", sep="/"), width=w, height=h, pointsize=p)
#explo.plot(mysaturation, toplot="protein_coding",
#           samples = c(1,3), yleftlim = NULL, yrightlim = NULL)
#dev.off()
#What about the depth of our samples?

}##########################################
{##Bias
##########################################
## Length bias detection
mylengthbias <- dat(mydata, factor="Group", norm=FALSE, type="lengthbias")
png(paste(PLOTSDIR, "lengthbias.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mylengthbias, samples=1:2)
dev.off()
cat("Lenght bias plot generated\n")
#Do we see a clear pattern?

##GC bias
mygcbiasRaw <- NOISeq::dat(mydata, factor = "Group", norm=FALSE, type="GCbias")
png(paste(PLOTSDIR, "GCbias.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mygcbiasRaw)
dev.off()
cat("GC bias plot generated\n")
#Do we see a clear pattern?

## RNA composition
mycomp <- dat(mydata, norm=FALSE, type="cd")
png(paste(PLOTSDIR, "RNAComposition.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mycomp, samples=1:12)
cat("RNA composition plot generated\n")
dev.off()
#Are the samples comparable?
##########################
## Quality Control Report
##########################
#A complete pdf report can be obtained using this function.
#QCreport(mydata, factor="Group", file=paste(PLOTSDIR, "QCReport.pdf", sep="/"))

###########################################################################
##Basal situation
##   GC content bias: Detected
##  Gene Length bias: Detected
##  RNA content bias: Detected
###########################################################################
}#############################
{## CQN Correction
# #############################
# ##Filter genes with mean < 10, i ., low expressionm genes
# expressedGenes<- apply(full$M, 1, function(x) mean(x)>10)
# 
# ##Filtering low expression genes
# mean10<-list(M=full$M[expressedGenes, ], Annot=full$Annot[expressedGenes, ],
#     Targets=full$Targets)
# 
# ##Length, GC content and size correction
# y_DESeq<-DESeqDataSetFromMatrix(countData=mean10$M, 
#     colData=mean10$Targets, design= ~Group)
# y_DESeq<-estimateSizeFactors(y_DESeq)
# 
# cqn.mean10<- cqn(mean10$M, lengths=mean10$Annot$Length,
#  x = mean10$Annot$GC, sizeFactors=sizeFactors(y_DESeq), verbose=TRUE)
# # RQ fit .....................................
# # SQN Using 'sigma' instead 'sig2' (= sigma^2) is preferred now
# 
# ##Normalized expression log2 scale
# normalized.cqn <- cqn.mean10$y + cqn.mean10$offset
# 
# ##Check the correction with NOISEQ
# mydataCQN <- NOISeq::readData(
#   data = normalized.cqn, 
#   length = mean10$Annot[, c("EnsemblID", "Length")], 
#   biotype = mean10$Annot[, c("EnsemblID", "Type")], 
#   chromosome = mean10$Annot[, c("Chr", "Start", "End")], 
#   factors = mean10$Targets[, "Group",drop=FALSE], 
#   gc = mean10$Annot[, c("EnsemblID", "GC")])
# 
# ## Length bias?
# mylengthbias <- NOISeq::dat(mydataCQN, factor = "Group", norm=TRUE,
#     type="lengthbias")
# png(paste(PLOTSDIR, "lengthbias_corrected.png", sep="/"),  width=w, height=h, pointsize=p)
# explo.plot(mylengthbias, samples = 1:2)
# dev.off()
# #Was it corrected?
# 
# ##GC bias?
# mygcbiasRaw <- NOISeq::dat(mydataCQN, factor = "Group", norm=TRUE, type="GCbias")
# png(paste(PLOTSDIR, "GCbias_corrected.png", sep="/"),  width=w, height=h, pointsize=p)
# explo.plot(mygcbiasRaw)
# dev.off()
# #Was it corrected?
# 
# ## RNA composition?
# mycomp <- dat(mydataCQN, norm=TRUE, type = "cd")
# png(paste(PLOTSDIR, "RNAComposition_corrected.png", sep="/"),  width=w, height=h, pointsize=p)
# explo.plot(mycomp, samples=1:12)
# dev.off()
# table(mycomp@dat$DiagnosticTest[,  "Diagnostic Test"])
# # FAILED PASSED 
# #     22     14 
# #Was it corrected? Which one is better?
# 
# ##Cleanning up
# #mean10 has the raw counts above mean > 10
# normalized.cqn<-list(M=normalized.cqn, Annot=mean10$Annot, 
#     Targets=mean10$Targets, cqn=cqn.mean10)
# 
# save(mean10, file=paste(RDATA, "mean10.RData", sep="/"), compress="xz")        
# save(normalized.cqn, file=paste(RDATA, "normalizedcqn.RData", sep="/"), compress="xz")    
    
}##########################################################################
{##Multidimesional PCA noise analysis
##########################################################################
## Distribution exploration
##########################################
# normalized.cqnMelt<-melt(normalized.cqn$M)
# names(normalized.cqnMelt) <- c("EntrezID", "Sample", "Expression")
# 
# ##Expression density
# p<-ggplot(data=normalized.cqnMelt,  
#     aes(x=Expression, group=Sample, colour=Sample))+geom_density()
# ggsave(paste(PLOTSDIR, "expression_densisty.png", sep="/"))
# 
# ##Expression boxplot
# p<-ggplot(data=normalized.cqnMelt,aes(y=Expression, x=Sample, group=Sample,    
#     colour=Sample))+geom_boxplot()
# ggsave(paste(PLOTSDIR, "expression_boxplot.png", sep="/"))
# 
# 
# ##########################################
# #### PCA EXPLORATION
# ##########################################
# pca.results<-prcomp(t(normalized.cqn$M))
# summary(pca.results)$importance[,1:10]
# 
# proportion<-round(summary(pca.results)$importance[2,]*100,0)
# ## Variance explained by each component
# png(paste(PLOTSDIR, "PCA.png", sep="/"), width=w, height=h, pointsize=p)
# screeplot(pca.results)
# dev.off()
# 
# ## Scatter plot on 1-3 PCA components
# Group <- normalized.cqn$Targets$Group
# p12<-ggplot(data=as.data.frame(pca.results$x), 
#         aes(x=PC1, y=PC2, colour=Group, shape=Group))+
#     geom_point(size=5)+
#     xlab(paste("PC 1 ", proportion[1], "%", sep = ""))+
#     ylab(paste("PC 2 ", proportion[2], "%", sep = ""))
# ggsave(paste(PLOTSDIR, "PCA_1-2.png", sep="/"))
# 
# p13<-ggplot(data=as.data.frame(pca.results$x), 
#         aes(x=PC1, y=PC3, colour=Group, shape=Group))+
#     geom_point(size=5)+
#     xlab(paste("PC 1 ", proportion[1], "%", sep = ""))+
#     ylab(paste("PC 3 ", proportion[3], "%", sep = ""))
# ggsave(paste(PLOTSDIR, "PCA_1-3.png", sep="/"))
#No need to filter samples, they can be accurately clustered
}##########################################################################
## GREAT JOB!!! YOU MADE IT TILL THE END!!!!
###########################################################################
