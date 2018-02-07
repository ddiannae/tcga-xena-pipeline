###############################################################################
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
load(file="RawFull.RData")
{###Let's keep only the GC & length annotated genes
ids<-!is.na(full$Annot$GC) & !is.na(full$Annot$Length)
table(ids)
# FALSE  TRUE 
#  1083 19449

full$M<-full$M[ids,]
full$Annot<-full$Annot[ids,]

dim(full$M)
# [1] 56963     6
dim(full$Annot)
# [1] 56963    10

}##########################################
## EXPLORATORY ANALYSIS (NOISeq package)
##########################################
{## Reading data into NOISeq package -> mydata
row.names(full$M)<-full$Annot$EnsemblID
row.names(full$Annot)<-full$Annot$EnsemblID
row.names(full$Targets)<-full$Targets$ID
full$Targets$Group<-factor(substr(full$Targets$ID, start=1, stop=1))
table(full$Targets$Group)
#   N   T 
#  3  3 

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
## Biodetection plot
mybiodetection <- dat(mydata, type="biodetection", factor="Group", k=0)
par(mfrow = c(1,1))
explo.plot(mybiodetection)
dev.off()
#What do we need to see here?

## Count distribution per biotype
mycountsbio <- dat(mydata, factor = NULL, type = "countsbio")
explo.plot(mycountsbio, toplot = 1, samples = 1, plottype = "boxplot")
#What about expression level?

## Saturation plot
mysaturation <- dat(mydata, k = 0, ndepth = 7, type = "saturation")
explo.plot(mysaturation, toplot="protein_coding", 
    samples = c(1,3), yleftlim = NULL, yrightlim = NULL)
#What about the depth of our samples?    

## Count distribution per sample
mycountsbio <- dat(mydata, factor = NULL, type = "countsbio")
explo.plot(mycountsbio, toplot = "protein_coding", 
    samples = NULL, plottype = "boxplot")
explo.plot(mycountsbio, toplot = "protein_coding", 
    samples = NULL, plottype = "barplot")
mycountsbio <- dat(mydata, factor = "Group", type = "countsbio")
## Count distribution per Experimental factors
explo.plot(mycountsbio, toplot = "protein_coding", 
    samples = NULL, plottype = "boxplot")
explo.plot(mycountsbio, toplot = "protein_coding", 
    samples = NULL, plottype = "barplot")
#How much sensitivity we loose? 

}##########################################
{##Bias
##########################################
## Length bias detection
mylengthbias <- dat(mydata, factor="Group", norm=FALSE, type="lengthbias")
par(mfrow = c(1,1))
explo.plot(mylengthbias, samples=1:2)
#Do we see a clear pattern?

##GC bias
mygcbiasRaw <- NOISeq::dat(mydata, factor = "Group", norm=FALSE, type="GCbias")
par(mfrow = c(1,2))
explo.plot(mygcbiasRaw)
#Do we see a clear pattern?

## RNA composition
mycomp <- dat(mydata, norm=FALSE, type="cd")
explo.plot(mycomp, samples=1:12)
#Are the samples comparable?

##########################
## Quality Control Report
##########################
#A complete pdf report can be obtained using this function.
#QCreport(mydata, factor="Group")

###########################################################################
##Basal situation
##   GC content bias: Detected
##  Gene Length bias: Detected
##  RNA content bias: Detected
###########################################################################
table(mycomp@dat$DiagnosticTest[,  "Diagnostic Test"])
# FAILED PASSED 
#     30      6 
   
}#############################
{## CQN Correction
#############################
##Filter genes with mean < 10, i ., low expressionm genes
expressedGenes<- apply(full$M, 1, function(x) mean(x)>10)
table(expressedGenes)
# FALSE  TRUE 
#  2431 17018

##Filtering low expression genes
mean10<-list(M=full$M[expressedGenes, ], Annot=full$Annot[expressedGenes, ],
    Targets=full$Targets)

##Length, GC content and size correction
y_DESeq<-DESeqDataSetFromMatrix(countData=mean10$M, 
    colData=mean10$Targets, design= ~Group)
y_DESeq<-estimateSizeFactors(y_DESeq)

cqn.mean10<- cqn(mean10$M, lengths=mean10$Annot$Length,
 x = mean10$Annot$GC, sizeFactors=sizeFactors(y_DESeq), verbose=TRUE)
# RQ fit .....................................
# SQN Using 'sigma' instead 'sig2' (= sigma^2) is preferred now

##Normalized expression log2 scale
normalized.cqn <- cqn.mean10$y + cqn.mean10$offset

##Check the correction with NOISEQ
mydataCQN <- NOISeq::readData(
  data = normalized.cqn, 
  length = mean10$Annot[, c("EnsemblID", "Length")], 
  biotype = mean10$Annot[, c("EnsemblID", "Type")], 
  chromosome = mean10$Annot[, c("Chr", "Start", "End")], 
  factors = mean10$Targets[, "Group",drop=FALSE], 
  gc = mean10$Annot[, c("EnsemblID", "GC")])

## Length bias?
mylengthbias <- NOISeq::dat(mydataCQN, factor = "Group", norm=TRUE,
    type="lengthbias")
par(mfrow = c(1,2))
explo.plot(mylengthbias, samples = 1:2)
#Was it corrected?

##GC bias?
mygcbiasRaw <- NOISeq::dat(mydataCQN, factor = "Group", norm=TRUE, type="GCbias")
par(mfrow = c(1,2))
explo.plot(mygcbiasRaw)
#Was it corrected?

## RNA composition?
mycomp <- dat(mydataCQN, norm=TRUE, type = "cd")
explo.plot(mycomp, samples=1:12)
table(mycomp@dat$DiagnosticTest[,  "Diagnostic Test"])
# FAILED PASSED 
#     22     14 
#Was it corrected? Which one is better?

##Cleanning up
#mean10 has the raw counts above mean > 10
normalized.cqn<-list(M=normalized.cqn, Annot=mean10$Annot, 
    Targets=mean10$Targets, cqn=cqn.mean10)
save(mean10, file="mean10.RData", compress="xz")        
save(normalized.cqn, file="normalizedcqn.RData", compress="xz")    
    
}##########################################################################
{##Multidimesional PCA noise analysis
##########################################################################
## Distribution exploration
##########################################
normalized.cqnMelt<-melt(normalized.cqn$M)
names(normalized.cqnMelt) <- c("EntrezID", "Sample", "Expression")

##Expression density
p<-ggplot(data=normalized.cqnMelt,  
    aes(x=Expression, group=Sample, colour=Sample))+geom_density()
p

##Expression boxplot
p<-ggplot(data=normalized.cqnMelt,aes(y=Expression, x=Sample, group=Sample,    
    colour=Sample))+geom_boxplot()
p

##########################################
#### PCA EXPLORATION
##########################################
pca.results<-prcomp(t(normalized.cqn$M))
summary(pca.results)$importance[,1:10]
# PC1      PC2      PC3      PC4      PC5      PC6
# Standard deviation     63.96480 42.87586 36.63182 34.20711 29.15849 27.80530
# Proportion of Variance  0.21175  0.09514  0.06945  0.06056  0.04400  0.04001
# Cumulative Proportion   0.21175  0.30689  0.37634  0.43690  0.48090  0.52091
# PC7      PC8      PC9     PC10
# Standard deviation     25.88505 24.07073 23.44858 22.92690
# Proportion of Variance  0.03468  0.02999  0.02846  0.02720
# Cumulative Proportion   0.55559  0.58557  0.61403  0.64123

proportion<-round(summary(pca.results)$importance[2,]*100,0)
## Variance explained by each component
screeplot(pca.results)


## Scatter plot on 1-3 PCA components
Group <- normalized.cqn$Targets$Group
p12<-ggplot(data=as.data.frame(pca.results$x), 
        aes(x=PC1, y=PC2, colour=Group, shape=Group))+
    geom_point(size=5)+
    xlab(paste("PC 1 ", proportion[1], "%", sep = ""))+
    ylab(paste("PC 2 ", proportion[2], "%", sep = ""))
p12

p13<-ggplot(data=as.data.frame(pca.results$x), 
        aes(x=PC1, y=PC3, colour=Group, shape=Group))+
    geom_point(size=5)+
    xlab(paste("PC 1 ", proportion[1], "%", sep = ""))+
    ylab(paste("PC 3 ", proportion[3], "%", sep = ""))
p13
#No need to filter samples, they can be accurately clustered
}##########################################################################
## GREAT JOB!!! YOU MADE IT TILL THE END!!!!
###########################################################################