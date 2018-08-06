##############################################################################
## CHROMATIN AND GENE REGULATION: FROM GENE TO GENOME FOLDING
## Practical session: In silico analysis of RNA-Seq data.
###############################################################################
## Data preparation: 
##      -Quality Control
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
##      -PCA
###############################################################################
##Usefull Libraries
###############################################################################
library("BiocParallel")
library("parallel")
library("NOISeq")
# register(SnowParam(workers=detectCores()-1, progress=TRUE))#Windows
register(MulticoreParam(workers=detectCores()-1, progress=TRUE))#Linux
options(width=80)
###############################################################################
##Quality Control
###############################################################################
cat("#################\n")
cat("Step 2: QC\n")
cat("#################\n")
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
cat("Non GC and lenght annotated genes removed.\n")

row.names(full$M)<-full$Annot$EnsemblID
row.names(full$Annot)<-full$AnnotdEnsemblID
row.names(full$Targets)<-full$Targets$ID
full$Targets$Group<-factor(substr(full$Targets$ID, start=1, stop=1))

save(full, file=paste(RDATA, "RawFull.RData", sep="/"), compress="xz")
cat("Saving RawFull.RData \n")
}##########################################
## EXPLORATORY ANALYSIS (NOISeq package)
##########################################
{## Reading data into NOISeq package -> mydata
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
# Biodetection plot. Per group.
mybiodetection <- dat(mydata, type="biodetection", factor="Group", k=0)
png(filename=paste(PLOTSDIR, "biodetection.Rd_%03d.png", sep="/"),  width=w, height=h, pointsize=p)
explo.plot(mybiodetection)
dev.off()
cat("Biodetection plots generated\n")
#What do we need to see here?

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

## Saturation plot. 
## We're not getting the saturation plot because it takes too long and doesn't give
## us that useful information.
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
mylengthbias <- dat(mydata, factor="Group", type="lengthbias")
png(paste(PLOTSDIR, "lengthbias.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mylengthbias, samples=1:2)
dev.off()
cat("Lenght bias plot generated\n")
#Do we see a clear pattern?

##GC bias
mygcbiasRaw <- dat(mydata, factor = "Group", type="GCbias")
png(paste(PLOTSDIR, "GCbias.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mygcbiasRaw)
dev.off()
cat("GC bias plot generated\n")
#Do we see a clear pattern?

## RNA composition
mycomp <- dat(mydata, type="cd")
png(paste(PLOTSDIR, "RNAComposition.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mycomp, samples=1:12)
dev.off()
cat("RNA composition plot generated\n")
#Are the samples comparable?
##########################
## Quality Control Report
##########################
#A complete pdf report can be obtained using this function.
#QCreport(mydata, factor="Group", file=paste(PLOTSDIR, "QCReport.pdf", sep="/"))

}#############################
{## PCA Analysis with NOISeq
##########################################
pca.dat <- dat(mydata, type = "PCA", logtransf = F)
pca.results <- pca.dat@dat$result

## Variance explained by each component
pdf(file=paste(PLOTSDIR, "PCAVariance_raw.pdf", sep="/"), width = 4*2, height = 4*2)
barplot(pca.results$var.exp[,1], xlab = "PC", ylab = "Explained variance")
dev.off()
cat("PCA variance raw plot generated.\n")

## Loading plot
pdf(file=paste(PLOTSDIR, "PCALoading_raw.pdf", sep="/"), width = 4*2, height = 4*2)
plot(pca.results$loadings[,1:2], col = 1, pch = 20, cex = 0.5,
     xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
     ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
     main = "PCA loadings",
     xlim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1),
     ylim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1))  
dev.off()
cat("PCA loading raw plot generated.\n")

## Score plot
mycol <- as.character(full$Targets$Group)
mycol[mycol == 'N'] <- "black"
mycol[mycol == 'T'] <- "red2"

pdf(file=paste(PLOTSDIR, "PCAScore_raw.pdf", sep="/"), width = 5*2, height = 5)
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
legend("topright", c("N", "T"),fill = c("black", "red2"), ncol = 2, pch = 1)

# PC1 & PC3
rango2 = diff(range(pca.results$scores[,c(1,3)]))
plot(pca.results$scores[,c(1,3)], col = "white",
   xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
   ylab = paste("PC 3 ", round(pca.results$var.exp[3,1]*100,0), "%", sep = ""),
   main = "PCA scores",
   xlim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1),
   ylim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1))
points(pca.results$scores[,1], pca.results$scores[,3], col = mycol, cex = 1.5)
legend("topright", c("N", "T"),fill = c("black", "red2"), ncol = 2, pch = 1)
dev.off()
cat("PCA scores raw plot generated.\n")

}#############################
##########################################################################
## GREAT JOB!!! YOU MADE IT TILL THE END!!!!
###########################################################################
cat("#################\n")
cat("End of Step 2: QC\n")
cat("#################\n")
