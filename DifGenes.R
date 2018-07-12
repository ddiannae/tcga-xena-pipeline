
#################################################################################
##3) DIFFERENTIAL EXPRESSION
#################################################################################
options(width=120)
library("NOISeq")
library("EDASeq")
load(file="FULLGC_FULLLength_TMM_CPM10_ARSYN.RData")

##Haciendo el análisis con lima
library("limma")
table(FULLGC_FULLLength_TMM_CPM10_ARSYN$Targets$Group)
#   E   S 
# 780 101
FULLGC_FULLLength_TMM_CPM10_ARSYN$Targets$Group<-relevel(FULLGC_FULLLength_TMM_CPM10_ARSYN$Targets$Group, ref="S")
design<-model.matrix(~1+Group, data=FULLGC_FULLLength_TMM_CPM10_ARSYN$Targets)
head(design)
#    (Intercept) GroupE
# S1           1      0
# S2           1      0
# E779           1      1
# E780           1      1

# y= mu + Enfermos + Error

M<-FULLGC_FULLLength_TMM_CPM10_ARSYN$M
row.names(M)<-FULLGC_FULLLength_TMM_CPM10_ARSYN$Annot$Symbol.y
fit <- lmFit(M, design)
head(fit$coefficients)
#         (Intercept)     GroupE
# GTPBP6    11.729239 -0.1324185
# EFCAB12    6.972360  0.2235230
# A1BG       9.699121  0.4780058
# GGACT      8.345224 -0.1782946
# A2ML1      6.933132  0.4942618
# A2M       17.729564 -0.9503528

##Enfermos-Sanos están en GroupE
fit2 <- eBayes(fit)
fit2$fdr<-apply(fit2$"p.value", 2, p.adjust, method="fdr")

library("ggplot2")
library("reshape2")

p<-ggplot(as.data.frame(fit$coefficients), aes(x=GroupE))+geom_density()
pdf(file="LFCDensity.pdf")
p
dev.off()

nrow(M)
# [1] 15281

##Buscando los genes diferenciales con ebayes---------------------------------
alphas<-c(0.05, 10^(-2:-10), 10^(seq(-20, -100, by=-10)))
degCount<-sapply(alphas, function(alpha){
  table(fit2$fdr[,"GroupE"]<alpha)
})
colnames(degCount)<-alphas
degCount<-as.data.frame(t(degCount))
degCount$alpha<-alphas
degCount$Noise<-alphas*nrow(M)
degCount

# Opcion 1)		Opcion 2)
# H0= lfc=0		H0= |lfc|<=l0
# H1= lfc!=0		H1= |lfc|>=l0
# 
# Opcion 1)
#       FALSE  TRUE alpha      Noise
# 0.05   2223 13058 5e-02 7.6405e+02
# 0.01   2884 12397 1e-02 1.5281e+02
# 0.001  3678 11603 1e-03 1.5281e+01
# 1e-04  4277 11004 1e-04 1.5281e+00
# 1e-05  4873 10408 1e-05 1.5281e-01
# 1e-06  5361  9920 1e-06 1.5281e-02
# 1e-07  5772  9509 1e-07 1.5281e-03
# 1e-08  6145  9136 1e-08 1.5281e-04
# 1e-09  6514  8767 1e-09 1.5281e-05
# 1e-10  6841  8440 1e-10 1.5281e-06
# 1e-20  9238  6043 1e-20 1.5281e-16
# 1e-30 10775  4506 1e-30 1.5281e-26
# 1e-40 11861  3420 1e-40 1.5281e-36
# 1e-50 12570  2711 1e-50 1.5281e-46 
# 1e-60 13134  2147 1e-60 1.5281e-56
# 1e-70 13551  1730 1e-70 1.5281e-66
# 1e-80 13889  1392 1e-80 1.5281e-76
# 1e-90  14139  1142  1e-90 1.5281e-86
# 1e-100 14333   948 1e-100 1.5281e-96   ##Una potencia muy elevada

pdf(file="vennDiagram.pdf")
vennDiagram(decideTests(fit2, adjust.method="fdr", p.value=1e-100))
dev.off()

##Con TREAT --------------------------------------
lfc<-seq(0, 3, by=0.5)
B<-5
fitTreat<-lapply(lfc, function(x){
  ##Ajustando el modelo
  aux<-treat(fit, lfc=x)
  aux$fdr<-apply(aux$"p.value", 2, p.adjust, method="fdr")
  ##Estadistico B
  p<-1-aux$fdr[, "GroupE"] ##Probabilidad de salir diferencialºº
  aux$B<-log(p/(1-p))
  
  ##Diferenciales
  aux$dif<-data.frame(
    alpha=alphas,
    lfc=rep(x, length(alphas)),
    Dif=t(sapply(alphas, function(alpha){
      table(factor(aux$fdr[,"GroupE"]<alpha & (aux$B>5), levels=c(TRUE, FALSE)))
    })),
    Noise=alphas*nrow(M)
  )
   
  return(aux)
})
names(fitTreat)<-paste("lfc", lfc, sep="")
difTreat<-do.call(rbind, lapply(fitTreat, function(x){x$dif}))
head(difTreat)
#        alpha lfc Dif.TRUE Dif.FALSE      Noise
# lfc0.1 5e-02   0    12234      3047 764.050000
# lfc0.2 1e-02   0    12234      3047 152.810000
# lfc0.3 1e-03   0    11603      3678  15.281000
# lfc0.4 1e-04   0    11004      4277   1.528100
# lfc0.5 1e-05   0    10408      4873   0.152810
# lfc0.6 1e-06   0     9920      5361   0.015281


main<-ggplot(difTreat, aes(x=alpha, y=Dif.TRUE, group=lfc)) + geom_line() + facet_grid(. ~ lfc)
main<- main + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_x_log10()
main<- main + geom_hline(aes(yintercept=1000), color="red") + geom_hline(aes(yintercept=2000), color="blue")
pdf(file="DifGenes.pdf")
main
dev.off()

subplot<- main %+% subset(difTreat, lfc == 1) 
pdf(file="Subplot.pdf")
subplot
dev.off()

library("grid")
pdf(file="DifComplete.pdf")
vp <- viewport(width=0.7, height=0.5, x=1, y=0.42, just = c("right", "bottom"))
print(main)
print(subplot, vp = vp)
dev.off()         
         
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

#############################################################################
## GOSeq
#############################################################################

library("goseq")

mygenes = rownames(myRaw)
# 
# deg1.goseq = as.integer(mygenes %in% mydeg1)
# deg2.goseq = as.integer(mygenes %in% mydeg2)
# 
# names(deg1.goseq) = names(deg2.goseq) = mygenes
# 
# deg.goseq = list("up" = deg1.goseq, "down" = deg2.goseq)

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
deg<-
mipwf<- nullp(deg.goseq[[i]], bias.data = mylength)  # Step 2

# enriched
# 
# head(walle[[1]])
# summary(walle[[1]][,2])


