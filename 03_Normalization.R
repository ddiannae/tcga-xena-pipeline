load(file="RawFull.RData")
#############################
#####Normalizacion con EDASeq
#############################
# source("http://bioconductor.org/biocLite.R")
# biocLite("EDASeq")
library("EDASeq")

##Seguimos recomendación de EDASEQ a priori de la corrección
##Filtramos por media > 10
genesExpressos<- apply(full$M, 1, function(x) mean(x)>10)
table(genesExpressos)
# FALSE  TRUE 
#  2234 17215

##Filtering low expression genes
mean10<-list(M=full$M[genesExpressos, ], Annot=full$Annot[genesExpressos, ], Targets=full$Targets)

# feature <- data.frame(gc=biomart[,5], length=biomart[,8], rownames=biomart[,6], 
#                       biotype = biomart[,4], 
#                       Chr=biomart[,1],
#                       GeneStart=biomart[,2],
#                       GeneEnd=biomart[,3])

##Normalization
mydataM10EDA <- EDASeq::newSeqExpressionSet(
  counts=mean10$M,
  featureData=mean10$Annot,
  phenoData=data.frame(
    conditions=mean10$Targets$Group,
    row.names=colnames(mean10$M)))
mydataM10EDA 

##########################################################################
##save.image(file="Norm.RData", compress="xz")
load(file="Norm.RData")
library("NOISeq")
library("EDASeq")

##withinLaneNormalization   #"loess","median","upper"
##betweenLaneNormalization  #"loess","median","upper"

#data2 <- withinLaneNormalization(as.matrix(counts(data)),feature$gc, which="full")
#data2 <- withinLaneNormalization(as.matrix(counts(data)),feature$length, 
#                                 which="full")
#data2 <- withinLaneNormalization(data2,feature$length, 
#                                 which="full")

## base<-"RPKM_FULL-GC_TMM" En 0
# data2 <- rpkm(as.matrix(counts(mydataM10EDA)), long=mean10$Annot$Length)
# data2 <- withinLaneNormalization(data2,mean10$Annot$GC, which="full")
# data3 <- tmm(data2, long = 1000, lc = 0, k = 0)

## base<-"FULL-Length_FULL-GC_TMM" #En 1
# data2 <- withinLaneNormalization(counts(mydataM10EDA),mean10$Annot$Length, which="full")
# data2 <- withinLaneNormalization(data2,mean10$Annot$GC, which="full")
# data3 <- tmm(data2, long = 1000, lc = 0, k = 0)

base<-"FULLGC_FULLLength_TMM" #En 2
data2 <- withinLaneNormalization(counts(mydataM10EDA),mean10$Annot$GC, which="full")
data2 <- withinLaneNormalization(data2,mean10$Annot$Length, which="full")
data3 <- tmm(data2, long = 1000, lc = 0, k = 0)

##base<-"FULL-GC_Loess-Length_TMM" #En 0
# data2 <- withinLaneNormalization(counts(mydataM10EDA),mean10$Annot$GC, which="full")
# data2 <- withinLaneNormalization(data2,mean10$Annot$Length, which="loess")
# data3 <- tmm(data2, long = 1000, lc = 1, k = 0)

## base<-"FULL-GC_FULL-Length_TMMlc1" #En 2
# data2 <- withinLaneNormalization(counts(mydataM10EDA),mean10$Annot$GC, which="full")
# data2 <- withinLaneNormalization(data2,mean10$Annot$Length, which="full")
# data3 <- tmm(data2, long = 1000, lc = 1, k = 0)

## base<-"FULL-GC_FULL-Length_BetweenFULL" #En 1
# data2 <- withinLaneNormalization(counts(mydataM10EDA),mean10$Annot$GC, which="full")
# data2 <- withinLaneNormalization(data2,mean10$Annot$Length, which="full")
# data3 <- betweenLaneNormalization(data2, which="full",offset=TRUE)

## base<-"FULL-GC_TMM_long-lc1" #En 3
# data2 <- withinLaneNormalization(counts(mydataM10EDA),mean10$Annot$GC, which="full")
# data3 <- tmm(data2, long = mean10$Annot$Length, lc = 1, k = 0)

##base<-"FULL-GC_FULL-Length_TMM_V2" #En 2
# data2 <- withinLaneNormalization(mydataM10EDA,"GC", which="full", num.bins=10, offset=FALSE, round=TRUE)
# data2 <- withinLaneNormalization(data2,"Length", which="full")
# data3 <- tmm(counts(data2), long = 1000, lc = 0)

mydata <- NOISeq::readData(
  data = data3, 
  length = mean10$Annot[, c("EntrezID", "Length")], 
  biotype = mean10$Annot[, c("EntrezID", "Type")], 
  chromosome = mean10$Annot[, c("Chr", "Start", "End")], 
  factors = mean10$Targets[, "Group",drop=FALSE], 
  gc = mean10$Annot[, c("EntrezID", "GC")])

## Length bias detection
mylengthbias <- dat(mydata, factor = "Group", norm = FALSE, type = "lengthbias")
png(file=paste(base, "lengthbias.png", sep=""))
par(mfrow = c(1,2))
explo.plot(mylengthbias, samples = 1:2)
dev.off()
##GC bias
mygcbiasRaw <- NOISeq::dat(mydata, factor = "Group", norm = FALSE, type = "GCbias")
png(file=paste(base, "GCBias.png", sep=""))
par(mfrow = c(1,2))
explo.plot(mygcbiasRaw)
dev.off()
## RNA composition
mycomp <- dat(mydata, norm = FALSE, type = "cd")
png(file=paste(base, "RNA.png", sep=""))
explo.plot(mycomp, samples=1:12)
dev.off()
table(mycomp@dat$DiagnosticTest[,  "Diagnostic Test"])

##
##Metodo				Length 	GC	RNA FAILED/PASSED
##"RPKM_FULL-GC_TMM"			Si	No	880/0 Ondulaciones laterales
##"FULL-Length_FULL-GC_TMM"		No++	No	872/8
##"FULL-GC_FULL-Length_TMM"		No++	No	867/13  Este??  Este es el mejorcito en distribuciones
##"FULL-GC_Loess-Length_TMM" 		Si	Si	866/14
##"FULL-GC_FULL-Length_TMMlc1" 		Si	Si	868/12
##"FULL-GC_FULL-Length_BetweenFULL"	No--	No-- 	872/8   Este??   
##"FULL-GC_TMM_long-lc1"		Si	Si	874/6

##QCreport(mydata, factor = "Group", norm=TRUE)
##Da error el reporte final

mycountsbio <- dat(mydata, factor = NULL, type = "countsbio", norm=FALSE)
## Count distribution per sample
png(file=paste(base, "Distribution.png", sep=""), width=2000, height=2000)
explo.plot(mycountsbio, toplot = "protein_coding", samples = NULL, plottype = "boxplot")
dev.off()
png(file=paste(base, "DistributionBar.png", sep=""), width=2000, height=2000)
explo.plot(mycountsbio, toplot = "protein_coding", samples = NULL, plottype = "barplot")
dev.off()

##En 4 el FULLGC_FULLLength_TMM
png(file=paste(base, "Boxplot.png", sep=""), width=2000, height=2000)
boxplot(log(data3+1))
dev.off()

##Limpiando
FULLGC_FULLLength_TMM<-mean10
FULLGC_FULLLength_TMM$M<-data3

## save(FULLGC_FULLLength_TMM, file="FULLGC_FULLLength_TMM.RData", compress="xz")
##########################################################################
##Multidimesional PCA noise analysis
##########################################################################
library("NOISeq")
library("EDASeq")
load(file="FULLGC_FULLLength_TMM.RData")

#### 3) NOISE FILTERING
library("ggplot2")
library("reshape")

p<-ggplot(data=melt(log(FULLGC_FULLLength_TMM$M+1)),
          aes(x=value, group=X2, colour=X2))+geom_density()
pdf(file="DensityRAWFullLog.pdf")
p
dev.off()

myfilterRaw<-filtered.data(FULLGC_FULLLength_TMM$M, factor=FULLGC_FULLLength_TMM$Targets$Group, 
                           norm=TRUE, cv.cutoff=100, cpm=10)
# [1] "Filtering out low count features..."
# [1] "15378 features are to be kept for differential expression analysis with filtering method 1" cpm=1
# [1] "15281 features are to be kept for differential expression analysis with filtering method 1" cpm=10 Corrige y no es tan agresivo
# [1] "15134 features are to be kept for differential expression analysis with filtering method 1" cpm=20
# [1] "15002 features are to be kept for differential expression analysis with filtering method 1" cpm=30
# [1] "14767 features are to be kept for differential expression analysis with filtering method 1" cpm=50

dim(FULLGC_FULLLength_TMM$M)
# [1] 17215   881
dim(myfilterRaw)
# [1] 15281   881
nrow(FULLGC_FULLLength_TMM$M)-nrow(myfilterRaw)
# [1] 1934 ##Filtrados por baja expression cpm=10
# [1] 2081 ##Filtrados por baja expression cpm=20
# [1] 2213 ##Filtrados por baja expression cpm=30
# [1] 2448 ##Filtrados por baja expression cpm=50

##Propagando a la anotacion
myfilterRaw<-list(M=myfilterRaw, 
                  Annot=FULLGC_FULLLength_TMM$Annot[row.names(FULLGC_FULLLength_TMM$Annot)%in%row.names(myfilterRaw),],
                  Targets=FULLGC_FULLLength_TMM$Targets)
stopifnot(nrow(myfilterRaw$M)==nrow(myfilterRaw$Annot))
stopifnot(row.names(myfilterRaw$M)==row.names(myfilterRaw$Annot))


p<-ggplot(data=melt(log(myfilterRaw$M+1)),aes(x=value, group=X2, colour=X2))+geom_density()
pdf(file="DensityFilter10.pdf")
p
dev.off()

#### 2) PCA EXPLORATION
setwd("./ARSyN")
source("sourceARSyN.R")
setwd("../")

#pca.results = PCA.GENES(t(log2(1+data3)))
pca.results <- PCA.GENES(t(log2(1+myfilterRaw$M)))
traditional.pca<-prcomp(t(log2(1+myfilterRaw$M)))
summary(traditional.pca)$importance[,1:10]
#                             PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8      PC9     PC10
# Standard deviation     46.35244 39.58127 26.07990 22.85858 20.19392 20.02734 17.23167 16.28167 14.59049 13.17168
# Proportion of Variance  0.14091  0.10275  0.04461  0.03427  0.02674  0.02631  0.01947  0.01739  0.01396  0.01138
# Cumulative Proportion   0.14091  0.24366  0.28826  0.32253  0.34928  0.37558  0.39506  0.41244  0.42640  0.43778

## Variance explained by each component
pdf("PCAVarianceRAW.pdf", width = 4*2, height = 4*2)
barplot(pca.results$var.exp[,1], xlab = "PC", ylab = "explained variance", ylim = c(0,0.4))
dev.off()

## Loading plot
pdf("PCALoadingRaw.pdf", width = 4*2, height = 4*2)
plot(pca.results$loadings[,1:2], col = 1, pch = 20, cex = 0.5,
     xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
     ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
     main = "PCA loadings",
     xlim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1),
     ylim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1))  
dev.off()

## Score plot
mycol <- as.character(myfilterRaw$Targets$Group)
mycol[mycol == 'S'] <- "black"
mycol[mycol == 'E'] <- "red2"

pdf("PCAScoreRaw.pdf", width = 5*2, height = 5)
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
legend("topright", unique(as.character(myfilterRaw$Targets$Group)), col = rep(unique(mycol),2), ncol = 2)

# PC1 & PC3
rango2 = diff(range(pca.results$scores[,c(1,3)]))
plot(pca.results$scores[,c(1,3)], col = "white",
     xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
     ylab = paste("PC 3 ", round(pca.results$var.exp[3,1]*100,0), "%", sep = ""),
     main = "PCA scores",
     xlim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1),
     ylim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1))
points(pca.results$scores[,1], pca.results$scores[,3], col = mycol, cex = 1.5)
legend("topright", unique(as.character(myfilterRaw$Targets$Group)), col = rep(unique(mycol),2), ncol = 2)

dev.off()

##Probamos si ARSyN reduce el ruido########################################################
myARSyN <- ARSyN(data=log2(myfilterRaw$M + 1), Covariates=t(as.matrix(myfilterRaw$Targets$Group)))

pca.results <- PCA.GENES(t(myARSyN))
traditional.pca<-prcomp(t(myARSyN))
summary(traditional.pca)$importance[,1:10]
#                             PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8     PC9     PC10
# Standard deviation     29.93246 5.611833 5.587665 5.541001 5.512378 5.460539 5.407716 5.380455 5.34742 5.324276
# Proportion of Variance  0.14271 0.005020 0.004970 0.004890 0.004840 0.004750 0.004660 0.004610 0.00455 0.004520
# Cumulative Proportion   0.14271 0.147720 0.152700 0.157590 0.162430 0.167170 0.171830 0.176440 0.18100 0.185510

## Variance explained by each component
pdf("PCAVarianceARSyN.pdf", width = 4*2, height = 4*2)
barplot(pca.results$var.exp[,1], xlab = "PC", ylab = "explained variance", ylim = c(0,0.4))
dev.off()

## Loading plot
pdf("PCALoadingARSyN.pdf", width = 4*2, height = 4*2)
plot(pca.results$loadings[,1:2], col = 1, pch = 20, cex = 0.5,
     xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
     ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
     main = "PCA loadings",
     xlim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1),
     ylim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1))  
dev.off()

## Score plot
mycol <- as.character(myfilterRaw$Targets$Group)
mycol[mycol == 'S'] <- "black"
mycol[mycol == 'E'] <- "red2"

pdf("PCAScoreARSyN.pdf", width = 5*2, height = 5)
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
legend("topright", unique(as.character(myfilterRaw$Targets$Group)), col = rep(unique(mycol),2), ncol = 2)

# PC1 & PC3
rango2 = diff(range(pca.results$scores[,c(1,3)]))
plot(pca.results$scores[,c(1,3)], col = "white",
     xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
     ylab = paste("PC 3 ", round(pca.results$var.exp[3,1]*100,0), "%", sep = ""),
     main = "PCA scores",
     xlim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1),
     ylim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1))
points(pca.results$scores[,1], pca.results$scores[,3], col = mycol, cex = 1.5)
legend("topright", unique(as.character(myfilterRaw$Targets$Group)), col = rep(unique(mycol),2), ncol = 2)

dev.off()


##Cómo quedaron las densidades y demás con esta modificación
pdf(file="BoxplotARSyN.pdf", width=2000, height=2000)
boxplot(myARSyN)
dev.off()

p<-ggplot(data=melt(myARSyN), aes(x=value, group=X2, colour=X2))+geom_density()
pdf(file="DensityARSyN.pdf")
p
dev.off()

dim(myARSyN)
# [1] 15281   881

##Limpiando
FULLGC_FULLLength_TMM_CPM10_ARSYN<-list(M=myARSyN, Annot=myfilterRaw$Annot, Targets=myfilterRaw$Targets)
dim(FULLGC_FULLLength_TMM_CPM10_ARSYN$M)
# [1] 15281   881
stopifnot(nrow(FULLGC_FULLLength_TMM_CPM10_ARSYN$M)==nrow(FULLGC_FULLLength_TMM_CPM10_ARSYN$Annot))
stopifnot(all(row.names(FULLGC_FULLLength_TMM_CPM10_ARSYN$M)==row.names(FULLGC_FULLLength_TMM_CPM10_ARSYN$Annot)))

##Matrices de datos para las redes
#Todos= sanos | enfermos
M<-as.data.frame(FULLGC_FULLLength_TMM_CPM10_ARSYN$M)
M<-cbind(gene=as.character(FULLGC_FULLLength_TMM_CPM10_ARSYN$Annot$Symbol.y), M)
#sanos
sanos<-as.data.frame(FULLGC_FULLLength_TMM_CPM10_ARSYN$M[,FULLGC_FULLLength_TMM_CPM10_ARSYN$Targets$Group=="S"])
sanos<-cbind(gene=as.character(FULLGC_FULLLength_TMM_CPM10_ARSYN$Annot$Symbol.y), sanos)
#enfermos
enfermos<-as.data.frame(FULLGC_FULLLength_TMM_CPM10_ARSYN$M[,FULLGC_FULLLength_TMM_CPM10_ARSYN$Targets$Group=="E"])
enfermos<-cbind(gene=as.character(FULLGC_FULLLength_TMM_CPM10_ARSYN$Annot$Symbol.y), enfermos)
#Symbolos
symbol<-as.character(FULLGC_FULLLength_TMM_CPM10_ARSYN$Annot$Symbol.y)
##Guardando
write.table(M, file = "M.tab", sep="\t", quote=FALSE, row.names=FALSE)
write.table(sanos, file = "Sanos.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(enfermos, file = "Enfermos.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(symbol, file = "ListaGenes.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
save(FULLGC_FULLLength_TMM_CPM10_ARSYN, file="FULLGC_FULLLength_TMM_CPM10_ARSYN.RData", compress="xz")

