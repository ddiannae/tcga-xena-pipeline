#################################################################################
## Proyecto: Mecanismos de desregulación del metabolismo energético asociados a 
##           cáncer de mama
################################################################################
## 3) DIFFERENTIAL EXPRESSION
## Autor: Cristóbal Fresno
## Fecha: 2019/02/05
## Modificado por : Erandi Serrano
## https://github.com/CSB-IG/tcgarnaseqbc/blob/master/DifGenes.R
#################################################################################
#options(width=120)

# Instalar las bibliotecas necesarias
# NOIseq
#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#    BiocManager::install("NOISeq", version = "3.8")
# EDAseq
#if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#      BiocManager::install("EDASeq", version = "3.8")
#      BiocManager::install("Glimma", version = "3.8")
library("NOISeq")
library("EDASeq")
library("Glimma")

#Insatalar limma
#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#   BiocManager::install("limma", version = "3.8")
library("limma")
library("edgeR")
library("ggplot2")

setwd("/mnt/ddisk/transpipeline-data/breast-data/subtypes/")
DEGDATA <- "deg"
dir.create(DEGDATA)
RDATA <- "rdata"
load(paste(RDATA, "data_subtypes_arsyn.RData", sep="/"))

tnames <- c("LumA", "LumB", "Her2", "Basal")
Matrix <- data.subtypes.arsyn$M
Targets <- data.subtypes.arsyn$Targets
## Acomodemos por subtipos
M <- cbind(Matrix[, Targets$Subtype == "Healthy"], Matrix[, Targets$Subtype == tnames[1]],
           Matrix[, Targets$Subtype == tnames[2]], Matrix[, Targets$Subtype == tnames[3]],
           Matrix[, Targets$Subtype == tnames[4]])

#### Correcciones y acomodado de los Data Frames
Group <- factor(c(as.character(Targets[Targets$Subtype == "Healthy", "Subtype"]),
                  as.character(paste0("t", Targets[Targets$Subtype == tnames[1], "Subtype"])),
                  as.character(paste0("t", Targets[Targets$Subtype == tnames[2], "Subtype"])),
                  as.character(paste0("t", Targets[Targets$Subtype == tnames[3], "Subtype"])),
                  as.character(paste0("t", Targets[Targets$Subtype == tnames[4], "Subtype"]))))
targets <- data.frame(ID = colnames(M), Group = Group)
LDEG <- list(M = cpm(M, log = TRUE), Targets = targets)
# LDEG <- list(M = M, Targets = targets)

######### DISENIO
# ¿Para el diseño correcto ya no hay interceptor "Group-1"?
# El resultado es que muchos genes son diferencialmente expresados
design <- model.matrix(~1+Group, data = LDEG$Targets)
rownames(design) <- LDEG$Targets$ID
#colnames(design) <-c("Healthy", "Basal", "Her2", "LumA", "LumB")
M <- LDEG$M

######### AJUSTE
fit <- lmFit(M, design)
head(fit$coefficients)

#########################
eval <- FALSE
######### EVALUAR P-Vals y LFC
######  Buscando los genes diferenciales cortados por los alphas
if(eval) {
  gname <- paste0("Group", paste0("t", tnames[1])) # NOmbre del grupo en el experimento, por ejemplo: GrouptLA
  alphas <- c(0.05, 10^(-2:-10), 10^(seq(-20, -100, by=-10)))
  
  ##### [2] Select genes diferenciales por el metodo de TREAT
  lfc <- seq(0, 3, by=0.5)
  B <- 5
  fitTreat <- lapply(lfc, function(x){
    ##Ajustando el modelo
    aux <- treat(fit, lfc=x)
    aux$fdr <- apply(aux$p.value, 2, p.adjust, method="fdr")
    ##Estadistico B
    p <- 1-aux$fdr[, gname] ##Probabilidad de salir diferencialºº
    aux$B <- log(p/(1-p))
    ##Diferenciales
    aux$dif <- data.frame(
      alpha = alphas,
      lfc = rep(x, length(alphas)),
      Dif = t(sapply(alphas, function(alpha){
        table(factor(aux$fdr[,gname]<alpha & (aux$B>5), levels=c(TRUE, FALSE)))
      })),
      Noise = alphas*nrow(M)
    )
    return(aux)
  })
  names(fitTreat) <- paste("lfc", lfc, sep="")
  difTreat <- do.call(rbind, lapply(fitTreat, function(x){x$dif}))
  difTreat
  
  # Pval <- 1e-50             # Basal pvalue = 1e-50
  # LFCsel <- fitTreat$lfc0   # Basal lfc = 0.0
  # Pval <- 1e-02             # LumA pvalue = 1e-02
  # LFCsel <- fitTreat$lfc0.5 # LumA lfc = 0.5
  # Pval <- 1e-40             # LumB pvalue = 1e-40
  # LFCsel <- fitTreat$lfc0   # LumB lfc = 0.0
  Pval <- 1             # LumB pvalue = 1e-40
  LFCsel <- fitTreat$lfc0   # LumB lfc = 0.0
  
  LFCsel$deg <- LFCsel$fdr[, gname, drop=FALSE] < Pval 
  #& (LFCsel$B>5)
  table(LFCsel$deg)
  ##############################################################
  
} else {
  #### Obtener valores para todos los genes
  # contr.matrix <- makeContrasts(
  #   HealthyvsBasal = Basal-Healthy, 
  #   HealthyvsHer2 = Healthy-Her2, 
  #   HealthyvsLumA = Healthy-LumA, 
  #   HealthyvsLumB = Healthy-LumB, 
  #   levels = colnames(design))
  # LFCFit <- contrasts.fit(fit, contr.matrix)
  LFCFit <- eBayes(fit)
  LFCFit$fdr <- apply(LFCFit$p.value, 2, p.adjust, method="fdr")
  p <- 1-LFCFit$fdr ##Probabilidad de salir diferencialºº
  LFCFit$B <- log(p/(1-p))
}


#Glimma plot
dt <- decideTests(LFCFit, lfc = 0) #### NO ENTENDI
#La ocupamos despues para pintar las redes y por eso la guardamos...
write.table(dt, file = "deg/dt-ebayes-all.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = T )
summary(dt)
glMDPlot(path = DEGDATA, LFCFit, counts = M, groups = Group, status = dt)

lapply(tnames, function(name) {
  tname <- paste0("Groupt", name)
  genesFULL<-cbind(
    Gene = rownames(LFCFit$coefficients),
    Coef = LFCFit$coefficients[,c("(Intercept)", tname)],
    p.value = LFCFit$p.value[, tname],
    FDR = LFCFit$fdr[, tname],
    B = LFCFit$B[, tname]
  )
  
  fname <- paste("ebayes-",tname,".tsv",sep = "")
  write.table(genesFULL, file =  paste(DEGDATA, fname, sep="/"), quote = FALSE, sep = "\t", row.names = F, col.names = TRUE )
})


