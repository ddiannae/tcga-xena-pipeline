
#################################################################################
##Subtypes con pbcmc
#################################################################################
options(width=120)
library("NOISeq")
library("EDASeq")
load(file="FULLGC_FULLLength_TMM.RData")
# load(file="FULLGC_FULLLength_TMM_CPM10_ARSYN.RData")

library("pbcmc")
library("BiocParallel")

##FULLGC_FULLLength_TMM
M<-FULLGC_FULLLength_TMM$M[, FULLGC_FULLLength_TMM$Targets$Group=="E"]
genes<-FULLGC_FULLLength_TMM$Annot[, c("EntrezID", "Symbol.y", "EntrezID")]

##FULLGC_FULLLength_TMM_CPM10_ARSYN
# M<-FULLGC_FULLLength_TMM_CPM10_ARSYN$M[, FULLGC_FULLLength_TMM_CPM10_ARSYN$Targets$Group=="E"]
# genes<-FULLGC_FULLLength_TMM_CPM10_ARSYN$Annot[, c("EntrezID", "Symbol.y", "EntrezID")]

# names(genes)<-c("probe", "NCBI.gene.symbol", "EntrezGene.ID")
# object<-PAM50(exprs=M, annotation=genes)
# object

##Antes ARSyN
# A PAM50 molecular permutation classifier object
# Dimensions:
#             nrow ncol
# exprs      17215  780
# annotation 17215    3
# targets        0    0


##ARSyN
# A PAM50 molecular permutation classifier object
# Dimensions:
#             nrow ncol
# exprs      15281  780
# annotation 15281    3
# targets        0    0

# object<-filtrate(object, verbose=TRUE)
# object<-classify(object, std="median", verbose=TRUE)
# object<-permutate(object, nPerm=10000, pCutoff=0.01, where="fdr",
#   corCutoff=0.1, keep=TRUE, seed=1234567890, verbose=TRUE,
#   BPPARAM=MulticoreParam(workers=4, progressbar=TRUE))

##FULLGC_FULLLength_TMM
##MEDIAN con cuentas crudas
# object
#  Basal   Her2   LumA   LumB Normal 
#    178    113    249    156     84 
# table(permutation(object)$subtype$Permuted)
#     Assigned Not Assigned    Ambiguous 
#          435          295           50 

#     Assigned Not Assigned    Ambiguous 
#    55.769231    37.820513     6.410256 

# PCA con los Asignados y no Asignados
# ARSyN 
# 
# 
# Probar que pasa con log(M+1)

object2<-PAM50(exprs=log(M+1), annotation=genes)
object2<-filtrate(object2, verbose=TRUE)
object2<-classify(object2, std="median", verbose=TRUE)
object2<-permutate(object2, nPerm=10000, pCutoff=0.01, where="fdr",
  corCutoff=0.1, keep=TRUE, seed=1234567890, verbose=TRUE,
  BPPARAM=MulticoreParam(workers=15, progressbar=TRUE))
  
object2  
# $subtype
# 
#  Basal   Her2   LumA   LumB Normal 
#    170    115    254    170     71 
table(permutation(object2)$subtype$Permuted)/780*100
# 
#     Assigned    Ambiguous Not Assigned 
#          467           82          231 
# 
#     Assigned    Ambiguous Not Assigned 
#     59.87179     10.51282     29.61538


##Log(conteos+1)
aux2<-table(permuted=permutation(object2)$subtype$Permuted, PAM50=permutation(object2)$subtype$PAM50)
aux2<-cbind(aux2, Total=rowSums(aux2))
aux2<-rbind(aux2, Total=colSums(aux2))
aux2
#              Basal Her2 LumA LumB Normal Total
# Assigned       142   72  163   58     32   467
# Ambiguous        2   14   32   11     23    82
# Not Assigned    26   29   59  101     16   231
# Total          170  115  254  170     71   780

##Crudos
# aux<-table(permuted=permutation(object)$subtype$Permuted, PAM50=permutation(object)$subtype$PAM50)
# aux<-cbind(aux, Total=rowSums(aux))
# aux<-rbind(aux, Total=colSums(aux))
# aux
#              Basal Her2 LumA LumB Normal Total
# Assigned       144   70  141   43     37   435
# Ambiguous        5    6   20    2     17    50
# Not Assigned    29   37   88  111     30   295
# Total          178  113  249  156     84   780

# table(permutation(object)$subtype$Permuted)/780*100
#     Assigned Not Assigned    Ambiguous 
#    55.769231    37.820513     6.410256 

##FULLGC_FULLLength_TMM_CPM10_ARSYN.RData  
##NONE
# object
# A PAM50 molecular permutation classifier object
# Dimensions:
#            nrow ncol
# exprs        41  780
# annotation   41    3
# targets       0    0
# Classification: 
#             nrow ncol
# probability  780    5
# correlation  780    5
# $subtype
# 
#  Basal   Her2   LumA   LumB Normal 
#      0      2    556    222      0 

# table(permutation(object)$subtype$Permuted)
# 
#     Assigned Not Assigned 
#           35          745 
  
# Obtaining 10000 permutations for 780 subjects... done.
# Error in `rownames<-`(x, value) : 
#   attempt to set 'rownames' on an object with no dimensions
  
##MEDIAN
# $subtype
# 
#  Basal   Her2   LumA   LumB Normal
#    167    144    179    159    131
# table(permutation(object)$subtype$Permuted)
# 
#     Assigned Not Assigned    Ambiguous 
#           44          734            2 
##ROBUST
#  Basal   Her2   LumA   LumB Normal 
#    178    130    173    161    138 
#    
# table(permutation(object)$subtype$Permuted)
# 
#     Assigned Not Assigned    Ambiguous 
#           46          733            1    

# save(object2, file="pbcmc.RData", compress="xz")

##########################################################################
##Multidimesional PCA noise analysis
##########################################################################
library("NOISeq")
library("EDASeq")

#### 3) NOISE FILTERING
library("ggplot2")
library("reshape")

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
# [1] 1934 ##Filtrados por baja expression cpm=10 ##Esta
# [1] 2081 ##Filtrados por baja expression cpm=20
# [1] 2213 ##Filtrados por baja expression cpm=30
# [1] 2448 ##Filtrados por baja expression cpm=50

##Propagando a la anotacion
myfilterRaw<-list(M=myfilterRaw, 
  Annot=FULLGC_FULLLength_TMM$Annot[row.names(FULLGC_FULLLength_TMM$Annot)%in%row.names(myfilterRaw),],
  Targets=FULLGC_FULLLength_TMM$Targets)
stopifnot(nrow(myfilterRaw$M)==nrow(myfilterRaw$Annot))
stopifnot(row.names(myfilterRaw$M)==row.names(myfilterRaw$Annot))

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

## Loading plot
# pdf("PCALoadingRaw.pdf", width = 4*2, height = 4*2)
# plot(pca.results$loadings[,1:2], col = 1, pch = 20, cex = 0.5,
#      xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
#      ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
#      main = "PCA loadings",
#      xlim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1),
#      ylim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1))  
# dev.off()

## Score plot
scores<-data.frame(pca.results$scores[,1:3])
names(scores)<-paste("PC",1:3, sep="")
scores$PBCMC<-as.character(myfilterRaw$Targets$Group)
scores$PBCMC[scores$PBCMC=="E"]<-as.character(permutation(object2)$subtype$Permuted)
scores$PBCMC<-as.factor(scores$PBCMC)
levels(scores$PBCMC)
# [1] "Ambiguous"    "Assigned"     "Not Assigned" "S"           
levels(scores$PBCMC)<-c("Amb", "OK", "NA", "Healthy")
scores$PAM50<-as.character(scores$PBCMC)
scores$PAM50[scores$PAM50!="Healthy"]<-as.character(permutation(object2)$subtype$PAM50)

# PC1 & PC2
rango <- diff(range(pca.results$scores[,1:2]))
p<-ggplot(data=scores, aes(x=PC1, y=PC2, shape=PBCMC, colour=PAM50))
p<- p + geom_point(size=3) + scale_shape(solid = FALSE) 
p<- p + xlab(paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""))
p<- p + ylab(paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""))
p<- p + ggtitle("PCA scores")
p<- p + xlim(range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1))
p<- p + ylim(range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1))

# PC1 & PC3
rango <- diff(range(pca.results$scores[,c(1,3)]))
p1<-ggplot(data=scores, aes(x=PC1, y=PC3, shape=PBCMC, colour=PAM50))
p1<- p1 + geom_point(size=3) + scale_shape(solid = FALSE) 
p1<- p1 + xlab(paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""))
p1<- p1 + ylab(paste("PC 3 ", round(pca.results$var.exp[3,1]*100,0), "%", sep = ""))
p1<- p1 + ggtitle("PCA scores")
p1<- p1 + xlim(range(pca.results$scores[,c(1,3)]) + 0.02*rango*c(-1,1))
p1<- p1 + ylim(range(pca.results$scores[,c(1,3)]) + 0.02*rango*c(-1,1))

# PC1 & PC2
rango <- diff(range(pca.results$scores[,1:2]))
p2<-ggplot(data=subset(scores, PBCMC%in%c("OK", "Healthy")), aes(x=PC1, y=PC2, shape=PBCMC, colour=PAM50))
p2<- p2 + geom_point(size=3) + scale_shape(solid = FALSE) 
p2<- p2 + xlab(paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""))
p2<- p2 + ylab(paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""))
p2<- p2 + ggtitle("PCA scores")
p2<- p2 + xlim(range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1))
p2<- p2 + ylim(range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1))

# PC1 & PC3
rango <- diff(range(pca.results$scores[,c(1,3)]))
p3<-ggplot(data=subset(scores, PBCMC%in%c("OK", "Healthy")), aes(x=PC1, y=PC3, shape=PBCMC, colour=PAM50))
p3<- p3 + geom_point(size=3) + scale_shape(solid = FALSE) 
p3<- p3 + xlab(paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""))
p3<- p3 + ylab(paste("PC 3 ", round(pca.results$var.exp[3,1]*100,0), "%", sep = ""))
p3<- p3 + ggtitle("PCA scores")
p3<- p3 + xlim(range(pca.results$scores[,c(1,3)]) + 0.02*rango*c(-1,1))
p3<- p3 + ylim(range(pca.results$scores[,c(1,3)]) + 0.02*rango*c(-1,1))


pdf("PCAScoreRawPBCMCFull.pdf", width = 5*2, height = 5)
multiplot(p,p1, cols=2)
dev.off()

pdf("PCAScoreRawPBCMCOK.pdf", width = 5*2, height = 5)
multiplot(p2,p3, cols=2)
dev.off()

##Nos quedamos solo con los sanos y asignados c("OK", "Healthy")
myfilterRaw$Targets$PBCMC<-scores$PBCMC
myfilterRaw$Targets$PAM50<-scores$PAM50
myfilterRaw$Targets$Clean<-scores$PBCMC%in%c("OK", "Healthy")

##Probamos si ARSyN reduce el ruido########################################################
myARSyN <- ARSyN(data=log2(myfilterRaw$M[,scores$PBCMC%in%c("OK", "Healthy")] + 1), Covariates=t(as.matrix(as.factor(myfilterRaw$Targets$PAM50[scores$PBCMC%in%c("OK", "Healthy")]))))


pca.results <- PCA.GENES(t(myARSyN))
traditional.pca<-prcomp(t(myARSyN))
summary(traditional.pca)$importance[,1:10]
#                             PC1      PC2      PC3     PC4      PC5      PC6     PC7      PC8      PC9     PC10
# Standard deviation     46.77241 42.84308 6.375118 6.32114 6.255885 6.242717 6.20067 6.167552 6.115202 6.090169
# Proportion of Variance  0.23298  0.19548 0.004330 0.00426 0.004170 0.004150 0.00409 0.004050 0.003980 0.003950
# Cumulative Proportion   0.23298  0.42846 0.432790 0.43704 0.441210 0.445360 0.44945 0.453500 0.457490 0.461440


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
scores<-data.frame(pca.results$scores[,1:3])
names(scores)<-paste("PC",1:3, sep="")
scores$PAM50<-myfilterRaw$Targets$PAM50[myfilterRaw$Targets$Clean]
scores$PBCMC<-myfilterRaw$Targets$PBCMC[myfilterRaw$Targets$Clean]

# PC1 & PC2
rango <- diff(range(pca.results$scores[,1:2]))
p<-ggplot(data=scores, aes(x=PC1, y=PC2, shape=PBCMC, colour=PAM50))
p<- p + geom_point(size=3) + scale_shape(solid = FALSE) 
p<- p + xlab(paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""))
p<- p + ylab(paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""))
p<- p + ggtitle("PCA scores")
p<- p + xlim(range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1))
p<- p + ylim(range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1))

# PC1 & PC3
rango <- diff(range(pca.results$scores[,c(1,3)]))
p1<-ggplot(data=scores, aes(x=PC1, y=PC3, shape=PBCMC, colour=PAM50))
p1<- p1 + geom_point(size=3) + scale_shape(solid = FALSE) 
p1<- p1 + xlab(paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""))
p1<- p1 + ylab(paste("PC 3 ", round(pca.results$var.exp[3,1]*100,0), "%", sep = ""))
p1<- p1 + ggtitle("PCA scores")
p1<- p1 + xlim(range(pca.results$scores[,c(1,3)]) + 0.02*rango*c(-1,1))
p1<- p1 + ylim(range(pca.results$scores[,c(1,3)]) + 0.02*rango*c(-1,1))

pdf("PCAScoreRawPBCMCARSyN.pdf", width = 5*2, height = 5)
multiplot(p,p1, cols=2)
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
# [1] 15281   568

##Limpiando
FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN<-list(M=myARSyN, Annot=myfilterRaw$Annot, Targets=myfilterRaw$Targets[myfilterRaw$Targets$Clean,])
dim(FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$M)
# [1] 15281   568
stopifnot(nrow(FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$M)==nrow(FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$Annot))
stopifnot(all(row.names(FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$M)==row.names(FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$Annot)))

##Matrices de datos para las redes
#Todos= sanos | enfermos
M<-as.data.frame(FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$M)
M<-cbind(gene=as.character(FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$Annot$Symbol.y), M)
#sanos
sanos<-as.data.frame(FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$M[,FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$Targets$Group=="S"])
sanos<-cbind(gene=as.character(FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$Annot$Symbol.y), sanos)
#enfermos
enfermos<-as.data.frame(FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$M[,FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$Targets$Group=="E"])
enfermos<-cbind(gene=as.character(FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$Annot$Symbol.y), enfermos)
#Subtypes
Subtypes<-lapply(levels(as.factor(FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$Targets$PAM50)), function(x){
    aux<-as.data.frame(FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$M[,FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$Targets$PAM50==x])
    aux<-cbind(gene=as.character(FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$Annot$Symbol.y), aux)
    return(aux)
})
names(Subtypes)<-levels(as.factor(FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$Targets$PAM50))
#Symbolos
symbol<-as.character(FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$Annot$Symbol.y)
##Guardando
# write.table(M, file = "M_PBCMC_ARSYM.tab", sep="\t", quote=FALSE, row.names=FALSE)
# write.table(sanos, file = "Sanos_PBCMC_ARSYM.txt", sep="\t", quote=FALSE, row.names=FALSE)
# write.table(enfermos, file = "Enfermos_PBCMC_ARSYM.txt", sep="\t", quote=FALSE, row.names=FALSE)
# write.table(symbol, file = "ListaGenes_PBCMC_ARSYM.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
# save(FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN, file="FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN.RData", compress="xz")
# lapply(names(Subtypes), function(x){
#     write.table(Subtypes[x], file = paste(x, "_PBCMC_ARSYM.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
#     x
# })

#################################################################################
##3) DIFFERENTIAL EXPRESSION
#################################################################################
options(width=120)
library("NOISeq")
library("EDASeq")
load(file="FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN.RData")

##Haciendo el análisis con lima
library("limma")
table(FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$Targets$PAM50)
#   Basal Healthy    Her2    LumA    LumB  Normal 
#     142     101      72     163      58      32 
FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$Targets$PAM50<-relevel(as.factor(FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$Targets$PAM50), ref="Healthy")
design<-model.matrix(~1+PAM50, data=FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$Targets)
head(design)
#    (Intercept) PAM50Basal PAM50Her2 PAM50LumA PAM50LumB PAM50Normal
# S1           1          0         0         0         0           0
# S2           1          0         0         0         0           0
# S3           1          0         0         0         0           0
# S4           1          0         0         0         0           0
# S5           1          0         0         0         0           0
# S6           1          0         0         0         0           0

# y= mu + PAM50 + Error

M<-FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$M
row.names(M)<-FULLGC_FULLLength_TMM_CPM10_PBCMC_Clean_ARSYN$Annot$Symbol.y
fit <- lmFit(M, design)
head(fit$coefficients)
#         (Intercept) PAM50Basal   PAM50Her2   PAM50LumA  PAM50LumB PAM50Normal
# GTPBP6    11.761080 -0.2606990 -0.18049925 -0.07800278 -0.1509403 -0.06756924
# EFCAB12    6.977990 -2.1363589 -0.05588216  1.19221251  0.7903373  0.25355947
# A1BG       9.634291 -0.1179835  0.43911346  0.63485551  0.6735042  0.26481812
# GGACT      8.367021 -0.1604120 -0.20778695 -0.17245763 -0.2306454 -0.09643754
# A2ML1      6.834358  3.6830125  1.12862136 -0.72701445  0.1078384  0.14838072
# A2M       18.150604 -0.8611522 -1.45724861 -1.36560033 -1.7249991 -0.71156439

# Opcion 1)     Opcion 2)
# H0= lfc=0     H0= |lfc|<=l0
# H1= lfc!=0        H1= |lfc|>=l0
# 

##Con TREAT --------------------------------------
lfc<-seq(0, 3, by=0.5)
alphas<-c(0.05, 10^(-2:-10), 10^(seq(-20, -100, by=-10)))
B<-5
fitTreat<-lapply(lfc, function(x){
  ##Ajustando el modelo
  aux<-treat(fit, lfc=x)
  aux$fdr<-apply(aux$"p.value"[, -1], 2, p.adjust, method="fdr")
  ##Estadistico B
  p<-1-aux$fdr ##Probabilidad de salir diferencialºº
  aux$B<-log(p/(1-p))
  
  ##Diferenciales
  aux$dif<-data.frame(
    alpha=alphas,
    lfc=rep(x, length(alphas)),
    Dif=t(sapply(alphas, function(alpha){
        colSums(aux$fdr<alpha & (aux$B>5))
    })),
    Dif.TRUE=sapply(alphas, function(alpha){
        sum(rowSums(aux$fdr<alpha & (aux$B>5))>0)
    }),
    Noise=alphas*nrow(M)
  )
  
  return(aux)
})
names(fitTreat)<-paste("lfc", lfc, sep="")
difTreat<-do.call(rbind, lapply(fitTreat, function(x){x$dif}))
head(difTreat)
#        alpha lfc Dif.PAM50Basal Dif.PAM50Her2 Dif.PAM50LumA Dif.PAM50LumB Dif.PAM50Normal Dif.TRUE      Noise
# lfc0.1 5e-02   0          12749         12087         12144         12140            6501    14568 764.050000
# lfc0.2 1e-02   0          12749         12087         12144         12140            6501    14568 152.810000
# lfc0.3 1e-03   0          12218         11453         11526         11483            5266    14271  15.281000
# lfc0.4 1e-04   0          11684         10745         10863         10829            4152    13898   1.528100
# lfc0.5 1e-05   0          11212         10122         10296         10282            3338    13526   0.152810
# lfc0.6 1e-06   0          10765          9656          9798          9798            2746    13176   0.015281

library("ggplot2")
library("reshape2")         
difTreat2<-melt(difTreat, id.vars=c("alpha", "lfc"))         
main<-ggplot(difTreat2, aes(x=alpha, y=value, group=variable, colour=variable)) + geom_line() + facet_grid(. ~ lfc)
main<- main + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_x_log10()
main<- main + geom_hline(aes(yintercept=1000), color="red") + geom_hline(aes(yintercept=2000), color="blue")
pdf(file="DifGenes.pdf")
main
dev.off()

subplot<- main %+% subset(difTreat2, lfc == 1) 
pdf(file="Subplot.pdf")
subplot
dev.off()

library("grid")
pdf(file="DifComplete.pdf")
vp <- viewport(width=0.7, height=0.5, x=1, y=0.42, just = c("right", "bottom"))
print(main)
print(subplot, vp = vp)
dev.off()         

###ACA!!!!!


##Buscando para un lfc=1
subset(difTreat, lfc==1)
#         alpha lfc Dif.TRUE Dif.FALSE      Noise
# lfc1.1  5e-02   1     1431     13850 7.6405e+02
# lfc1.2  1e-02   1     1431     13850 1.5281e+02 Este es igual al anterior!!!! |lfc|>1 & fdr<0.01 & B>5
# lfc1.3  1e-03   1     1330     13951 1.5281e+01
# lfc1.4  1e-04   1     1254     14027 1.5281e+00
# lfc1.5  1e-05   1     1185     14096 1.5281e-01
# lfc1.6  1e-06   1     1116     14165 1.5281e-02
# lfc1.7  1e-07   1     1060     14221 1.5281e-03
# lfc1.8  1e-08   1      997     14284 1.5281e-04
# lfc1.9  1e-09   1      946     14335 1.5281e-05
# lfc1.10 1e-10   1      907     14374 1.5281e-06
# lfc1.11 1e-20   1      649     14632 1.5281e-16
# lfc1.12 1e-30   1      497     14784 1.5281e-26
# lfc1.13 1e-40   1      386     14895 1.5281e-36
# lfc1.14 1e-50   1      313     14968 1.5281e-46
# lfc1.15 1e-60   1      257     15024 1.5281e-56
# lfc1.16 1e-70   1      212     15069 1.5281e-66
# lfc1.17 1e-80   1      178     15103 1.5281e-76
# lfc1.18 1e-90   1      147     15134 1.5281e-86
# lfc1.19 1e-100  1      124     15157 1.5281e-96

##Usar B=5 es igual a usar
exp(5)/(exp(5)+1)
# [1] 0.9933071
1-(exp(5)/(exp(5)+1))
# [1] 0.006692851
##Tomemos |lfc|>1 & fdr<0.01 & B>5
fitTreat$lfc1$deg<-fitTreat$lfc1$fdr[,"GroupE", drop=FALSE]<0.01 & (fitTreat$lfc1$B>5)
table(fitTreat$lfc1$deg)
# FALSE  TRUE 
# 13850  1431 

##Comprobando supuestos a nivel global
ajustados<-fitTreat$lfc1$coefficients%*%t(design)
residuos<-M-ajustados
pdf(file="supuestos.pdf")
ggplot(melt(as.data.frame(residuos)),aes(x=value))+geom_density()
dev.off()

pdf(file="QQresiduos.pdf")
qqnorm(as.numeric(residuos))
qqline(as.numeric(residuos), col = 2)
dev.off()
##Lo cumple a nivel global una pinturita

#HO:es normal H1:No es normal
testResult<-apply(residuos, 1, shapiro.test)
testResult<-as.data.frame(do.call(rbind,testResult))
testResult$data.name<-NULL
table(unlist(testResult$p.value)<0.01)
# FALSE  TRUE 
#   670 14611 
#OJO!!! No cumplen los supuestos casi la mayoría de los genes

table(unlist(testResult$p.value[fitTreat$lfc1$deg[, "GroupE"]])<0.01)
# FALSE  TRUE 
#   195   836 ##En los diferenciales poco y nada 
##Inspeccionando algunos genes
idResiduo<-which(fitTreat$lfc1$deg[, "GroupE"])[10]
pdf(file="ResiduoGen10.pdf")
qqnorm(residuos[idResiduo,])
qqline(residuos[idResiduo,],col="red")
dev.off()

## heatmap Diferenciales
library("gplots")
pdf(file="Heatmap Enfermo vs sanos.pdf",  width=14, height=14)
heatmap.2(scale(t(M[fitTreat$lfc1$deg[, "GroupE"], ])), col=greenred, scale="none", trace="none", 
    main = "Enfermos vs Sanos = 1341 genes")
dev.off()

#############################################################################
## Generando salida
#############################################################################
genesFULL<-cbind(
  FULLGC_FULLLength_TMM_CPM10_ARSYN$Annot, 
  Coef=fitTreat$lfc1$coefficients,
  Diff=fitTreat$lfc1$deg[, "GroupE"], 
  "p.value"=do.call(cbind, lapply(fitTreat, function(x){x$"p.value"[, "GroupE"]})),
  FDR=do.call(cbind, lapply(fitTreat, function(x){x$fdr[, "GroupE"]})), 
  B=do.call(cbind, lapply(fitTreat, function(x){
    ##Estadistico B
    p<-1-x$fdr[, "GroupE"] ##Probabilidad de salir diferencial
    B<-log(p/(1-p))
    return(B)
  })),
  Exp=fitTreat$lfc1$coefficients%*%t(unique(design)))

head(genesFULL)
#        EntrezID Symbol.x Row Chr     Start       End    GC           Type     HGNCID Symbol.y Length Coef.(Intercept)
# 8225       8225        ?  28   X    304529    318819 57.63 protein_coding HGNC:30189   GTPBP6  14290        11.729239
# 90288     90288        ?  29   3 129401321 129428651 49.22 protein_coding HGNC:28061  EFCAB12  27330         6.972360
# 1             1     A1BG  30  19  58345178  58353499 55.83 protein_coding     HGNC:5     A1BG   8321         9.699121
# 87769     87769    A2LD1  33  13 100530164 100589528 45.71 protein_coding HGNC:25100    GGACT  59364         8.345224
# 144568   144568    A2ML1  34  12   8822472   8887001 44.25 protein_coding HGNC:23336    A2ML1  64529         6.933132
# 2             2      A2M  35  12   9067664   9116229 37.18 protein_coding     HGNC:7      A2M  48565        17.729564
#        Coef.GroupE  Diff p.value.lfc0 p.value.lfc0.5 p.value.lfc1 p.value.lfc1.5 p.value.lfc2 p.value.lfc2.5
# 8225    -0.1324185 FALSE 9.430800e-04   1.000000e+00    1.0000000              1            1              1
# 90288    0.2235230 FALSE 4.902726e-03   9.997453e-01    1.0000000              1            1              1
# 1        0.4780058 FALSE 5.266291e-24   6.839252e-01    1.0000000              1            1              1
# 87769   -0.1782946 FALSE 1.939596e-03   1.000000e+00    1.0000000              1            1              1
# 144568   0.4942618 FALSE 1.399856e-05   5.202184e-01    0.9999956              1            1              1
# 2       -0.9503528 FALSE 6.429964e-45   1.711047e-12    0.7816307              1            1              1
#        p.value.lfc3     FDR.lfc0  FDR.lfc0.5 FDR.lfc1 FDR.lfc1.5 FDR.lfc2 FDR.lfc2.5 FDR.lfc3    B.lfc0 B.lfc0.5 B.lfc1
# 8225              1 1.234682e-03 1.00000e+00        1          1        1          1        1  6.695707     -Inf   -Inf
# 90288             1 6.139858e-03 1.00000e+00        1          1        1          1        1  5.086795     -Inf   -Inf
# 1                 1 1.449724e-23 1.00000e+00        1          1        1          1        1       Inf     -Inf   -Inf
# 87769             1 2.492345e-03 1.00000e+00        1          1        1          1        1  5.992036     -Inf   -Inf
# 144568            1 2.022809e-05 1.00000e+00        1          1        1          1        1 10.808418     -Inf   -Inf
# 2                 1 3.131175e-44 1.14728e-11        1          1        1          1        1       Inf 25.19104   -Inf
#        B.lfc1.5 B.lfc2 B.lfc2.5 B.lfc3    Exp.S1    Exp.E1
# 8225       -Inf   -Inf     -Inf   -Inf 11.729239 11.596821
# 90288      -Inf   -Inf     -Inf   -Inf  6.972360  7.195883
# 1          -Inf   -Inf     -Inf   -Inf  9.699121 10.177127
# 87769      -Inf   -Inf     -Inf   -Inf  8.345224  8.166929
# 144568     -Inf   -Inf     -Inf   -Inf  6.933132  7.427394
# 2          -Inf   -Inf     -Inf   -Inf 17.729564 16.779211
  
# genesDif<-genesFULL[genesFULL$Diff, c("EntrezID", "Symbol.y", "Coef.GroupE", "Exp.S1", "Exp.E1", "p.value.GroupE",
#   "FDR.GroupE", "B")]
# names(genesDif)<-c("EntrezID", "Symbol", "lofFC", "Exp.S", "Exp.E", "p.value", "fdr", "B")  
# head(genesDif)
#       EntrezID Symbol     lofFC     Exp.S     Exp.E       p.value           fdr         B
# 10157    10157   AASS -1.241122 11.175543  9.934421  1.525122e-08  2.159906e-07 15.348031
# 10349    10349 ABCA10 -2.494541  9.483414  6.988873 3.574866e-129 7.382097e-127       Inf
# 23461    23461  ABCA5 -1.423001 11.689791 10.266789  5.261781e-11  8.535591e-10 20.881606
# 23460    23460  ABCA6 -2.332611 10.688962  8.356351 1.521264e-110 2.152448e-108       Inf
# 10350    10350  ABCA9 -3.148426 10.543619  7.395193 2.919128e-181 1.939443e-178       Inf
# 368        368  ABCC6 -1.439275  9.183789  7.744514  2.118443e-04  2.347493e-03  6.052057

# write.table(genesFULL, file="genesFULL_LFC.tab", sep="\t", row.names=FALSE, quote=FALSE)  
# write.table(genesDif, file="genesDif.tab", sep="\t", row.names=FALSE, quote=FALSE)  

save.image(file="GenesFULL.RData", compress="xz")

