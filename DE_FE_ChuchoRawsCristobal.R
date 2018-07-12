##########################################################################################################
##############        Differential Expression and Functional Enrichment for RNA-Seq         ##############
##############                          EXERCISE 1                                          ##############
##########################################################################################################


# By Sonia Tarazona
# September, 2014




## This is a guided exercise. You only have to run the code and observe the results.



#### DATA

## Simulated data from A.fumigatus experiment
# 2 experimental conditions, 3 BIOLOGICAL replicates per condition
# 5% DEG
 source("http://bioconductor.org/biocLite.R")
# biocLite("edgeR")

setwd("~/RawCountsNotron")
#myRaw = read.delim("RawChico.txt", header = TRUE, as.is = TRUE, row.names = 1)
myRaw = read.delim("RawFull.txt", header = TRUE, as.is = TRUE, row.names = 1)

head(myRaw)
 dim(myRaw)

# myRaw = myRaw[, grep(".D", colnames(myRaw), fixed = TRUE)]
myRaw = myRaw[,order(colnames(myRaw))]
# tail(myRaw)

## Biological annotation (from Biomart)
biomart = read.delim("biomartHH.txt", header = TRUE, as.is = TRUE, row.names =1)
 tail(biomart)






##*****************************************************##



#### 1) EXPLORATORY ANALYSIS (NOISeq package)

 biocLite("NOISeq")
library(NOISeq)

## Reading data in NOISeq package

mybiotypes = biomart[,c(6,4)]
 head(mybiotypes)

mychroms = biomart[,1:3]
 head(mychroms)
rownames(mychroms) = biomart[,6]
colnames(mychroms) = c("Chr", "GeneStart", "GeneEnd")

myGC = biomart[,c(6,5)]
 head(myGC)
##En esta dices cuantos son los grupos, en each
myfactors = data.frame("group" = substr(colnames(myRaw), start = 1, stop = 1)
#                        , "day" = substr(colnames(myRaw), start = 2, stop = 2)
)
# myfactors
# myfactors = data.frame(myfactors, "cond" = apply(myfactors, 1, paste, collapse = "_"))
# myfactors = data.frame("group" = as.factor(rep(c("E","S"), each=9)))

# myfactors
mylength = biomart[,c(6, 8)]
 head(mylength)

mydata = NOISeq::readData(data = myRaw, length = mylength, biotype = mybiotypes, chromosome = mychroms, factors = myfactors, gc = myGC)




# head(myfactors, n =12)
#############################################
#############################################


##########################################
####TODO ESTO VA EN EL QCreport ##########
##########################################
## Biodetection plot
# mybiodetection <- dat(mydata, type = "biodetection", factor = "group", k = 0)
# par(mfrow = c(1,2))
# explo.plot(mybiodetection)
# ## Count distribution per biotype
# mycountsbio = dat(mydata, factor = NULL, type = "countsbio")
# explo.plot(mycountsbio, toplot = 1, samples = 1, plottype = "boxplot")
# ## Saturation plot
# mysaturation = dat(mydata, k = 0, ndepth = 7, type = "saturation")
# explo.plot(mysaturation, toplot = "protein_coding", samples = c(1,4), yleftlim = NULL, yrightlim = NULL)
# ## Count distribution per sample
# explo.plot(mycountsbio, toplot = "global", samples = NULL, plottype = "boxplot")
# explo.plot(mycountsbio, toplot = "global", samples = NULL, plottype = "barplot")
# ## Length bias detection
# mylengthbias = dat(mydata, factor = "group", norm = FALSE, type = "lengthbias")
# par(mfrow = c(1,2))
# explo.plot(mylengthbias, samples = 1:2)
# ## RNA composition
# mycomp = dat(mydata, norm = FALSE, type = "cd")
# explo.plot(mycomp)
########################################
########################################

## Quality Control Report
QCreport(mydata, factor = "group")




##*****************************************************##


#####Normalizacion con EDASeq
#############################
# biocLite("ggplot2")
 biocLite("reshape")
library("ggplot2")
library("reshape")
ggplot(data=melt(log(myRaw+1)),aes(x=value, group=variable, colour=variable))+geom_density()
ggplot(data=melt(log(uqua.counts+1)),aes(x=value, group=X2, colour=X2))+geom_density()


# source("http://bioconductor.org/biocLite.R")
# biocLite("EDASeq")
# help.start()

library(EDASeq)

dim(myGC)
head(mylength)

#feature <- data.frame(gc=myGC[,2],length=mylength[,2], rownames=biomart[,6], 
#                      biotype = mybiotypes, chromosome = mychroms)

feature <- data.frame(gc=biomart[,5], length=biomart[,8], rownames=biomart[,6], 
                      biotype = biomart[,4], 
                      Chr=biomart[,1],
                      GeneStart=biomart[,2],
                      GeneEnd=biomart[,3])

#filtrar por media mayor a un valor

#Hacemos que tengan la misma cantidad de filas
myRaw<-myRaw[row.names(myRaw)%in%feature$rownames, ]
feature<-feature[feature$rownames %in%row.names(myRaw), ]

#hacemos que tengan el mismo orden las filas
myRaw<-myRaw[order(row.names(myRaw)),]
feature<-feature[order(feature$rownames),]
row.names(feature)<-feature$rownames
all(row.names(myRaw)==feature$rownames)

common<- apply(myRaw,1,function(x) mean(x)>10)

table(common)
#FALSE  TRUE 
#2326 15332

##Filtering low expression genes
feature<-feature[common,]
myRaw<-as.matrix(myRaw[common,])

##Raw data QC
mydataRaw = NOISeq::readData(data=myRaw, 
                          length = feature[,c("rownames","length")], 
                          biotype = feature$mybiotypes, 
                          chromosome = feature[,c("Chr", "GeneStart", "GeneEnd")], 
                          factors = phenoData(data)@data, 
                          gc = feature[,c("rownames","gc")])

## Length bias detection
mylengthbiasRaw = NOISeq::dat(mydataRaw, factor = "conditions", norm = TRUE, type = "lengthbias")
par(mfrow = c(1,2))
explo.plot(mylengthbiasRaw)
##GCBais
mygcbiasRaw = NOISeq::dat(mydataRaw, factor = "conditions", norm = TRUE, type = "GCbias")
par(mfrow = c(1,2))
explo.plot(mygcbiasRaw)
## RNA composition
mycompRaw = NOISeq::dat(mydataRaw, norm = TRUE, type = "cd", refColumn =5)
explo.plot(mycompRaw, samples=sample(1:ncol(data2), 12))
table(mycompRaw@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#382     44

##Length, gc and RNA content bias detected

##Normalization
data <- EDASeq::newSeqExpressionSet(counts=myRaw,
                            featureData=feature,
                            phenoData=data.frame(
                              conditions=c(rep("E",327),rep("S",100)),
                              row.names=colnames(myRaw)))
data

#"loess","median","upper"
#data2 <- withinLaneNormalization(as.matrix(counts(data)),feature$gc, which="full")
#data2 <- withinLaneNormalization(as.matrix(counts(data)),feature$length, 
#                                 which="full")
#data2 <- withinLaneNormalization(data2,feature$length, 
#                                 which="full")

#RPKM + FC gc
data2 <- rpkm(as.matrix(counts(data)), long=feature$length)
data2 <- withinLaneNormalization(data2,feature$gc, which="full")

#dataOffset <- withinLaneNormalization(data,"length", which="full")
#dataOffset <- betweenLaneNormalization(data2,
#                                       which="full",offset=TRUE)

data3 <- tmm(data2, long = 1000, lc = 1, k = 0)
#data2 <- tmm(data2, long=feature$length)

mydata = NOISeq::readData(data=data3, 
                          length = feature[,c("rownames","length")], 
                          biotype = feature$mybiotypes, 
                          chromosome = feature[,c("Chr", "GeneStart", "GeneEnd")], 
                          factors = phenoData(data)@data, 
                          gc = feature[,c("rownames","gc")])


## Length bias detection
mylengthbias = NOISeq::dat(mydata, factor = "conditions", norm = TRUE, type = "lengthbias")
par(mfrow = c(1,2))
explo.plot(mylengthbias)
##GCBias
mygcbias = NOISeq::dat(mydata, factor = "conditions", norm = TRUE, type = "GCbias")
par(mfrow = c(1,2))
explo.plot(mygcbias)
## RNA composition
mycomp = NOISeq::dat(mydata, norm = TRUE, type = "cd", refColumn =5)
explo.plot(mycomp, samples=sample(1:ncol(data2), 12))

table(mycomp@dat$DiagnosticTest[,  "Diagnostic Test"])
#PASSED 
#426 
#write.table(data3, file = "NormMatrixRNAseq.txt")

load(file="DataNormalized.RData")
sanos<-data3[,328:ncol(data3)]
sanos<-cbind(gene=row.names(sanos),sanos)
enfermos<-data3[,1:327]
enfermos<-cbind(gene=row.names(enfermos),enfermos)
#write.table(sanos, file = "Sanos.tab", sep="\t", quote=FALSE, row.names=FALSE)
#write.table(enfermos, file = "Enfermos.tab", sep="\t", quote=FALSE, row.names=FALSE)

newEnfermos<-read.table(file="EnfermosNormalized.txt", sep="\t", header=TRUE)
head(newEnfermos)[,1:5]
head(enfermos)[,1:5]
enfermos<-as.data.frame(enfermos)
enfermos<-data.frame(enfermos, stringsAsFactors = FALSE)

dim(enfermos)==dim(newEnfermos)
all(enfermos[,1] == newEnfermos[,1])
all(enfermos[,-1] == newEnfermos[,-1])
enfermos[,-1][1:5,1:5] - newEnfermos[,-1][1:5,1:5]
enfermos[,-1]<-apply(enfermos[,-1], 2, function(x){as.numeric(as.character(x))})
head(enfermos)[,1:5]

rm(data3)
#save.image(file="CompleteDataNormalized.RData", compress="xz")
#save(list=c("myRaw", "data2", "data3", "feature"),file="DataNormalized.RData", compress="xz")
load(file="DataNormalized.RData")

library("ggplot2")
library("reshape")
boxplot(data3[,1:10])
boxplot(log(data3+1))

ggplot(data=melt(log(myRaw+1)),aes(x=value, group=X2, colour=X2))+geom_density()
ggplot(data=melt(log(data2+1)),aes(x=value, group=X2, colour=X2))+geom_density()
ggplot(data=melt(log(data2[,sample(1:ncol(data2), 12)]+1)),aes(x=value, group=X2, colour=X2))+geom_density()

ggplot(data=melt(log(data3[,sample(1:ncol(data2), 12)]+1)),aes(x=value, group=Var2, colour=Var2))+geom_density()

ggplot(data=melt(log(newEnfermos[ ,sample(2:ncol(newEnfermos), 12)]+1)),aes(x=value, group=variable, colour=variable))+geom_density()

ggplot(data=melt(log(newEnfermos[ ,sample(2:ncol(newEnfermos), 12)]+1)),aes(y=value, x=variable, group=variable, colour=variable))+geom_boxplot()

ggplot(data=melt((newEnfermos[ ,sample(2:ncol(newEnfermos), 12)])),aes(y=value, x=variable, group=variable, colour=variable))+geom_boxplot()

ggplot(data=melt(log(newEnfermos[ ,sample(2:ncol(newEnfermos), 12)]+1)),aes(y=value, x=variable, group=variable, colour=variable))+geom_boxplot()


ggplot(data=melt(log(uqua.counts+1)),aes(x=value, group=X2, colour=X2))+geom_density()

##Probando un gen
enfermos<-data3[,1:327]
id<-which(row.names(enfermos)%in%c("PYY", "OPRL1"))


plot(t(enfermos[id,]))
#aux<-prcomp(log((data2)+1))
#biplot(aux)
#plot(aux$rotation[,1], aux$rotation[,2], col=c(rep("red",3), rep("black",3)))


##Multidimesional PCA noise analysis

#### 3) NOISE FILTERING
library("NOISeq")
myfilterRaw<-filtered.data(data3, factor="conditions", norm=TRUE, cv.cutoff=100, cpm=10)
nrow(data3)
#[1] 15573
#[1] "Filtering out low count features..."
#[1] "9268 features are to be kept for differential expression analysis with filtering method 1" cpm=1
#[1] "6432 features are to be kept for differential expression analysis with filtering method 1" cpm=10
#15573-9268
#6305 filtrados

ggplot(data=melt(log(myfilterRaw[ ,sample(1:ncol(myfilterRaw), 12)]+1)),aes(y=value, x=X2, group=X2, colour=X2))+geom_boxplot()

ggplot(data=melt(log(myfilterRaw[ ,sample(1:ncol(myfilterRaw), 12)]+1)),aes(x=value, group=Var2, colour=Var2))+geom_density()


#### 2) PCA EXPLORATION

source("sourceARSyN.R")

#pca.results = PCA.GENES(t(log2(1+data3)))
pca.results = PCA.GENES(t(log2(1+myfilterRaw)))
traditional.pca<-prcomp(t(log2(1+myfilterRaw)))
summary(traditional.pca)

## Variance explained by each component
# pdf("Ex2_PCAvarExp14Oct.pdf", width = 4*2, height = 4*2)
barplot(pca.results$var.exp[,1], xlab = "PC", ylab = "explained variance", ylim = c(0,0.4))
# dev.off()

## Loading plot
# pdf("Ex2_LoadingPlot14Oct.pdf", width = 4*2, height = 4*2)
plot(pca.results$loadings[,1:2], col = 1, pch = 20, cex = 0.5,
     xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
     ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
     main = "PCA loadings",
     xlim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1),
     ylim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1))  
# dev.off()


## Score plot

# shapes for the plot
# mypch = rep(c(16,17), each = 6)

myfactors = data.frame("type" = substr(colnames(data2), start = 1, stop = 1),
                       "day" = substr(colnames(data2), start = 3, stop = 5)
)
myfactors$type
myfactors = data.frame(myfactors, "cond" = apply(myfactors, 1, paste, collapse = "_"))
# colors for the plot
mycol = as.character(myfactors$type)
mycol[mycol == 'S'] = "black"
mycol[mycol == 'E'] = "red2"


# pdf("Ex2_PCA.pdf", width = 5*2, height = 5)
par(mfrow = c(1,2))

# PC1 & PC2
rango = diff(range(pca.results$scores[,1:2]))
plot(pca.results$scores[,1:2], col = "white",
     xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
     ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
     main = "PCA scores",
     xlim = range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1),
     ylim = range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1))
points(pca.results$scores[,1], pca.results$scores[,2], col = mycol, cex = 1.5)  
legend("topright", unique(as.character(myfactors$type)),), col = rep(unique(mycol),2), ncol = 2)


# PC1 & PC3
rango2 = diff(range(pca.results$scores[,c(1,3)]))
plot(pca.results$scores[,c(1,3)], col = "white",
     xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
     ylab = paste("PC 3 ", round(pca.results$var.exp[3,1]*100,0), "%", sep = ""),
     main = "PCA scores",
     xlim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1),
     ylim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1))
points(pca.results$scores[,1], pca.results$scores[,3], col = mycol, cex = 1.5)
legend("topright", unique(as.character(myfactors$type)), col = rep(unique(mycol),2), ncol = 2)


# dev.off()
##Probamos si ArSym reduce el ruido########################################################
myARSyN <- ARSyN(data=log2(myfilterRaw + 1), Covariates=t(as.matrix(myfactors$type)))

#pca.results = PCA.GENES(t(log2(1+data3)))
pca.results = PCA.GENES(t(myARSyN))
traditional.pca<-prcomp(t(myARSyN))
summary(traditional.pca)

## Variance explained by each component
# pdf("Ex2_PCAvarExp14Oct.pdf", width = 4*2, height = 4*2)
barplot(pca.results$var.exp[,1], xlab = "PC", ylab = "explained variance", ylim = c(0,0.4))
# dev.off()

## Loading plot
# pdf("Ex2_LoadingPlot14Oct.pdf", width = 4*2, height = 4*2)
plot(pca.results$loadings[,1:2], col = 1, pch = 20, cex = 0.5,
     xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
     ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
     main = "PCA loadings",
     xlim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1),
     ylim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1))  
# dev.off()


## Score plot

# shapes for the plot
# mypch = rep(c(16,17), each = 6)

myfactors = data.frame("type" = substr(colnames(data2), start = 1, stop = 1),
                       "day" = substr(colnames(data2), start = 3, stop = 5)
)
myfactors$type
myfactors = data.frame(myfactors, "cond" = apply(myfactors, 1, paste, collapse = "_"))
# colors for the plot
mycol = as.character(myfactors$type)
mycol[mycol == 'S'] = "black"
mycol[mycol == 'E'] = "red2"


# pdf("Ex2_PCA.pdf", width = 5*2, height = 5)
par(mfrow = c(1,1))
# PC1 & PC2
rango = diff(range(pca.results$scores[,1:2]))
plot(pca.results$scores[,1:2], col = "white",
     xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
     ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
     main = "PCA scores",
     xlim = range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1),
     ylim = range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1))
points(pca.results$scores[,1], pca.results$scores[,2], col = mycol, cex = 1.5)  
legend("topright", legend=unique(as.character(myfactors$type)), col=unique(mycol), pch=21)


ggplot(data=melt((myARSyN[ ,sample(1:ncol(myARSyN), 12)])),aes(y=value, x=X2, group=X2, colour=X2))+geom_boxplot()

ggplot(data=melt((myARSyN[ ,sample(1:ncol(myARSyN), 12)])),aes(x=value, group=X2, colour=X2))+geom_density()

dim(myARSyN)
#[1] 6432  427

sanos<-myARSyN[,328:ncol(myARSyN)]
sanos<-cbind(gene=row.names(sanos),sanos)
enfermos<-myARSyN[,1:327]
enfermos<-cbind(gene=row.names(enfermos),enfermos)
#write.table(sanos, file = "Sanos.tab", sep="\t", quote=FALSE, row.names=FALSE)
#write.table(enfermos, file = "Enfermos.tab", sep="\t", quote=FALSE, row.names=FALSE)

#save(list=c("myRaw", "data2", "data3", "feature", "myARSyN"),
#  file="DataNormalized.RData", compress="xz")
load(file="DataNormalized.RData")


#### 3) DIFFERENTIAL EXPRESSION
## edgeR
##Realizado por Jesus
# # library("edgeR")
# # myRawF = filtered.data(data3, factor ="conditions")
# # #                          myfactors[,"type"])
# # myedgeR = DGEList(counts = myRawF, group = c(rep("E",327),rep("S",100)))
# # #                     myfactors[,"type"])
# # myedgeR = calcNormFactors(myedgeR)  
# # myedgeR = estimateCommonDisp(myedgeR)
# # myedgeR = estimateTagwiseDisp(myedgeR, trend = "movingave")
# # myedgeR = exactTest(myedgeR)
# # myedgeR = topTags(myedgeR, n = nrow(myRawF), sort.by = "logFC")[[1]]
# # dim(myedgeR)
# # 
# # myedgeR.deg = myedgeR[which(myedgeR[,"FDR"] < 1e-1), ] ##aprox 1000 genes
# # det = nrow(myedgeR.deg)
# # det
# # # write.table(myedgeR.deg, file = 'DEGsCasiFull_FDR1_lfc_1.txt', sep ="\t", quote = FALSE)
# # 
# # 
# # mydata = readData(data = data3, biotype = as.matrix(feature$biotype), chromosome = feature[,c("Chr", "GeneStart", "GeneEnd")], factors = "conditions", gc = feature$gc)
# # mynoiseqbio = noiseqbio(myedgeR, norm = "n", k = 0.5, lc = 0, factor = "type", conditions = NULL, r = 20, plot = TRUE, filter = 1)
# # mynoiseqbio.deg = degenes(mynoiseqbio, q=0.99999)
# # det = nrow(mynoiseqbio.deg)
# # ##visualizacion
# # pdf("Ex2_myNoiseqbioNoNorm.pdf", width = 7*3, height = 3*3)
# # par(mfrow = c(1,2))
# # DE.plot(mynoiseqbio, q = 0.95, graphic = "expr", log.scale = TRUE)
# # DE.plot(mynoiseqbio, q = 0.95, graphic = "MD", log.scale = TRUE)
# # dev.off()
# # write.csv(mynoiseqbio.deg, file = 'noiseqbioNoNorm.csv')

##Haciendo el an??lisis con lima
library("limma")
target<-data.frame(treatment=factor(substr(start=1, stop=1, colnames(myARSyN)), level=c("S", "E")),
  replicate=colnames(myARSyN))
table(target$treatment)
#   S   E 
# 100 327
dim(myARSyN)
# [1] 6432  427
design<-model.matrix(~1+treatment, data=target)
head(design)
#   (Intercept) treatmentE
# 1           1          1
# 427           1          0


fit <- lmFit(myARSyN, design)
head(fit$coefficients)
#       (Intercept) treatmentE
# A1BG    0.7676851  1.7098267
# A2M     9.2661820 -0.9901541
# AAAS    4.2051544  0.1241219
# AACS    5.1092921 -0.3309201
# AAGAB   3.9097905  0.7914879
# AAMP    6.9543866  0.1454249
##Enfermos-Sanos están en treatmentE
fit2 <- eBayes(fit)
fit2$fdr<-apply(fit2$"p.value", 2, p.adjust, method="fdr")

##Buscando los genes diferenciales---------------------------------
alphas<-c(0.05, 10^(-2:-40))
degCount<-sapply(alphas, function(alpha){
  table(fit2$fdr[,"treatmentE"]<alpha)
})
colnames(degCount)<-alphas
degCount
#       0.05 0.01 0.001 1e-04 1e-05 1e-06 1e-07 1e-08 1e-09 1e-10 1e-11 1e-12 1e-13 1e-14 1e-15 1e-16 1e-17 1e-18 1e-19
# FALSE 1087 1402  1740  2005  2271  2518  2715  2890  3065  3222  3360  3520  3645  3774  3895  4011  4107  4203  4277
# TRUE  5345 5030  4692  4427  4161  3914  3717  3542  3367  3210  3072  2912  2787  2658  2537  2421  2325  2229  2155
#       1e-20 1e-21 1e-22 1e-23 1e-24 1e-25 1e-26 1e-27 1e-28 1e-29 1e-30 1e-31 1e-32 1e-33 1e-34 1e-35 1e-36 1e-37 1e-38
# FALSE  4385  4454  4521  4605  4671  4732  4820  4908  4963  5023  5070  5103  5153  5196  5246  5285  5321  5364  5401
# TRUE   2047  1978  1911  1827  1761  1700  1612  1524  1469  1409  1362  1329  1279  1236  1186  1147  1111  1068  1031
#       1e-39 1e-40
# FALSE  5451  5489
# TRUE    981   943
vennDiagram(decideTests(fit2, adjust.method="fdr", p.value=1e-38))
##Tomemos alpha=1e-38 que son 1031 diferenciales y 5401 no diferenciales
fit2$deg<-cbind(limma38=fit2$fdr<1e-38)

##comparando con los de edgeR sin usar myARSyN
edgeRDeg<-read.table(file="DEGs_FDR10.txt", sep="\t", header=TRUE, as.is=TRUE)
fit2$deg<-cbind(fit2$deg, EdgeR=row.names(fit2$coefficients)%in%edgeRDeg$Gene)
table(limma=fit2$deg[,"treatmentE"], edger=fit2$deg[,"EdgeR"])
#        edger
# limma   FALSE TRUE
#   FALSE  5213  188
#   TRUE    543  488  Hay algunos que no me coninciden!!!! porque se uso data3 como base si filtrar los de baja expresion

##Comprobando supuestos a nivel global
ajustados<-fit2$coefficients%*%t(design)
residuos<-myARSyN-ajustados
ggplot(melt(as.data.frame(residuos)),aes(x=value))+geom_density()##Lo cumpple a nivel global

#HO:es normal H1:No es normal
testResult<-apply(residuos, 1, shapiro.test)
testResult<-as.data.frame(do.call(rbind,testResult))
testResult$data.name<-NULL
table(unlist(testResult$p.value)<0.01)
# FALSE  TRUE 
#   963  5469
#OJO!!! No cumplen los supuestos casi la mayor??a
table(unlist(testResult$p.value[fit2$deg[,"treatmentE"]])<0.01)
# FALSE  TRUE 
#   195   836 ##En los diferenciales poco y nada 
##Inspeccionando algunos genes
idResiduo<-which(fit2$deg[,"treatmentE"])[100]
qqnorm(residuos[idResiduo,])
qqline(residuos[idResiduo,],col="red")
#Pasando a NOISeqBIO


##Salida con las columnas Gene  logFC PValue  FDR


### A VER C??MO SACO ESTAS LINEAS
# TP = length(intersect(rownames(myedgeR.deg), simdata[[2]][,"gene"]))
# TP/nrow(simdata[[2]])  # Sensitivity = 0.969574
# (det - TP)/det  # False Discovery Rate = 0.1180812
###############################


## NOISeq
##### Este es para datos sin replicas biologicas,
##### pero nosotros tenemos muchas replicas, por 
##### lo tanto, no lo usamos
# mynoiseq = noiseq(mydata, norm = "tmm", k = 0.5, lc = 0, replicates = "biological", factor = "group", conditions = NULL, pnr = 0.2, nss = 5, v = 0.02)
# mynoiseq.deg = degenes(mynoiseq, q=0.9)
# head(mynoiseq.deg)
# det = nrow(mynoiseq.deg)
# TP = length(intersect(rownames(mynoiseq.deg), simdata[[2]][,"gene"]))
# TP/nrow(simdata[[2]])  # Sensitivity = 0.7910751
# (det - TP)/det  # False Discovery Rate = 0.02985075
###############################################

library("NOISeq")
## NOISeqBIO
### Este si, es para manejar varias replicas biologicas
### Sin normalizaci??n

# mydata<- readData(data=myARSyN, biotype=feature[], chromosome = mychroms, factors = myfactors, gc = myGC)
myfactors <- data.frame("type" = substr(colnames(myARSyN), start = 1, stop = 1),
                       "day" = substr(colnames(myARSyN), start = 3, stop = 5)
)
myfactors <- data.frame(myfactors, "cond" = apply(myfactors, 1, paste, collapse = "_"))
idGenes<-row.names(feature)%in%row.names(myARSyN)
mydata<- NOISeq::readData(data=myARSyN, 
                          length = feature[idGenes, c("rownames","length")], 
                          biotype = feature$mybiotypes[idGenes], 
                          chromosome = feature[idGenes, c("Chr", "GeneStart", "GeneEnd")], 
                          factors = myfactors, 
                          gc = feature[idGenes, c("rownames","gc")])

mynoiseqbio <- noiseqbio(mydata, norm = "n", k = 0.5, lc = 0, factor = "type", 
  conditions = NULL, r = 1000, plot = TRUE, filter = 0)

alphas<-c(0.05, 10^(-2:-40))
results<-as.data.frame(cbind(alpha=alphas, do.call(rbind, sapply(alphas, function(alpha){
  table(mynoiseqbio@results[[1]]$prob>(1-alpha))  
}))))
results[results[,"TRUE"]==nrow(myARSyN), "TRUE"]<-0
results
#    alpha FALSE TRUE
# 1  5e-02  1317 5115
# 2  1e-02  1944 4488
# 3  1e-03  2892 3540
# 4  1e-04  3646 2786
# 5  1e-05  4127 2305
# 6  1e-06  4340 2092
# 7  1e-07  4402 2030
# 8  1e-08  4438 1994
# 9  1e-09  4471 1961
# 10 1e-10  4492 1940
# 11 1e-11  4501 1931
# 12 1e-12  4509 1923
# 13 1e-13  4518 1914
# 14 1e-14  4535 1897
# 15 1e-15  4566 1866
# 16 1e-16  5138 1294 ##Usar este corte y todos o los primeros 1000 rankeados en este corte
# 17 1e-17  6432    0
# 18 1e-18  6432    0
# 19 1e-19  6432    0
# 20 1e-20  6432    0

mynoiseqbio.degChucho <- degenes(mynoiseqbio, q=1-(1e-10))
# [1] "1294 differentially expressed features"

det <- nrow(mynoiseqbio.degChucho)
##visualizacion
pdf("NoiseqbioNormFullChucho.pdf", width = 7*3, height = 3*3)
par(mfrow = c(1,2))
DE.plot(mynoiseqbio, q = 1-(1e-10), graphic = "expr", log.scale = FALSE)
DE.plot(mynoiseqbio, q = 1-(1e-10), graphic = "MD", log.scale = FALSE)
dev.off()
# [1] "1448 differentially expressed features"
write.csv(mynoiseqbio.degChucho, file = 'noiseqbioDEGNormalizadosFullChucho2.csv')

##Comparando los resultados de diferenciales
fit2$deg<-cbind(fit2$deg, NOISeqBIO=row.names(fit2$deg)%in%row.names(mynoiseqbio.deg))
fit2$deg<-cbind(fit2$deg, NOISeqBIO1000=row.names(fit2$deg)%in%row.names(mynoiseqbio.deg)[1:1000])
table(limma=fit2$deg[,"treatmentE"], edger=fit2$deg[,"EdgeR"])
#        edger
# limma   FALSE TRUE
#   FALSE  5213  188
#   TRUE    543  488
table(limma=fit2$deg[,"treatmentE"], NOISeqBIO=fit2$deg[,"NOISeqBIO"])
#        NOISeqBIO
# limma   FALSE TRUE
#   FALSE  4698  703
#   TRUE    440  591

## Solo con los primeros 1000
table(limma=fit2$deg[,"treatmentE"], NOISeqBIO1000=fit2$deg[,"NOISeqBIO1000"])
#        NOISeqBIO1000
# limma   FALSE TRUE
#   FALSE  4863  538
#   TRUE    569  462
table(edgeRDeg$Gene%in%row.names(mynoiseqbio.deg))
# FALSE  TRUE 
#   679   405 ##Solo se comparten 405 de los genes

###Hasta aca se realiz?? el an??lisis



# ### para el rpkm
# mynoiseqbioRPKM = noiseqbio(mydata, norm = "rpkm", k = 0.5, lc = 0, factor = "group", conditions = NULL, r = 30, plot = TRUE, filter = 1)
# mynoiseqbioRPKM.deg = degenes(mynoiseqbioRPKM, q=0.99999)
# det = nrow(mynoiseqbioRPKM.deg)
# ##visualizacion
# pdf("Ex2_myNoiseqbioRPKM.pdf", width = 7*3, height = 3*3)
# par(mfrow = c(1,2))
# DE.plot(mynoiseqbioRPKM, q = 0.99999, graphic = "expr", log.scale = TRUE)
# DE.plot(mynoiseqbioRPKM, q = 0.99999, graphic = "MD", log.scale = TRUE)
# dev.off()
# write.csv(mynoiseqbioRPKM.deg, file = 'noiseqbioDEGS_RPKM.csv')

# ### para el uqua
# mynoiseqbioUQUA = noiseqbio(mydata, norm = "uqua", k = 0.5, lc = 0, factor = "group", conditions = NULL, r = 30, plot = TRUE, filter = 1)
# mynoiseqbioUQUA.deg = degenes(mynoiseqbioUQUA, q=0.99999)
# det = nrow(mynoiseqbioUQUA.deg)
# ##visualizacion
# pdf("Ex2_myNoiseqbioUQUA.pdf", width = 7*3, height = 3*3)
# par(mfrow = c(1,2))
# DE.plot(mynoiseqbioUQUA, q = 0.99999, graphic = "expr", log.scale = TRUE)
# DE.plot(mynoiseqbioUQUA, q = 0.99999, graphic = "MD", log.scale = TRUE)
# dev.off()
# write.csv(mynoiseqbioUQUA.deg, file = 'noiseqbioDEGS_UQUA.csv')


# ### para el tmm
# mynoiseqbioTMM = noiseqbio(mydata, norm = "tmm", k = 0.5, lc = 0, factor = "group", conditions = NULL, r = 30, plot = TRUE, filter = 1)
# mynoiseqbioTMM.deg = degenes(mynoiseqbioTMM, q=0.99999)
# det = nrow(mynoiseqbioTMM.deg)
# ##visualizacion
# pdf("Ex2_myNoiseqbioTMM.pdf", width = 7*3, height = 3*3)
# par(mfrow = c(1,2))
# DE.plot(mynoiseqbioTMM, q = 0.99999, graphic = "expr", log.scale = TRUE)
# DE.plot(mynoiseqbioTMM, q = 0.99999, graphic = "MD", log.scale = TRUE)
# dev.off()
# write.csv(mynoiseqbioTMM.deg, file = 'noiseqbioDEGS_TMM.csv')
# TP = length(intersect(rownames(mynoiseqbio.deg), simdata[[2]][,"gene"]))
# TP/nrow(simdata[[2]])  # Sensitivity = 0.8864097
# (det - TP)/det  # False Discovery Rate = 0.02017937
##*****************************************************##



#### 4) Visualization of differential expression results
## Expression plot
## Manhattan plot
# DE.plot(mynoiseq, chromosomes = c("I", "II"), log.scale = TRUE, join = FALSE, q = 0.8, graphic = "chrom")


## DE per chromosomes and biotypes
# DE.plot(mynoiseqbio, chromosomes = NULL, q = 0.95, graphic = "distr")




##*****************************************************##




#### 5) FUNCTIONAL ENRICHMENT ANALYSIS


## Functional annotation (GO terms)

# GOannot <- read.delim("HS_GOannot.txt", as.is = TRUE, header = FALSE)
# head(GOannot)
# 
# 
# ## DE lists (up vs rest & down vs rest)   --> Step 1
# 
# mydeg1 = rownames(degenes(mynoiseqbio, q = 0.99, M = "up"))
# 
# mydeg2 = rownames(degenes(mynoiseqbio, q = 0.99, M = "down"))
# 
# mygenes = rownames(myRaw)
# 
# deg1.goseq = as.integer(mygenes %in% mydeg1)
# deg2.goseq = as.integer(mygenes %in% mydeg2)
# 
# names(deg1.goseq) = names(deg2.goseq) = mygenes
# 
# deg.goseq = list("up" = deg1.goseq, "down" = deg2.goseq)
# 
# 
# 
# 
# ## GOseq
# 
# library(goseq)
# 
# 
# mipwf <- walle <- enriched <- vector("list", length = 2)
# names(mipwf) <- names(walle) <- names(enriched) <- c("up", "down")
# 
# 
# for (i in 1:2) {
#   
#   mipwf[[i]] <- nullp(deg.goseq[[i]], bias.data = mylength)  # Step 2
#   
#   walle[[i]] <- goseq(pwf = mipwf[[i]], gene2cat = GOannot, method = "Wallenius", repcnt = 2000)  # Step 3
#   
#   enriched[[i]] <- walle[[i]]$category[p.adjust(walle[[i]]$over_represented_pvalue, method = "BH") < 0.01]
#   
# }  
# 
# enriched
# 
# head(walle[[1]])
# summary(walle[[1]][,2])
