##Remove self loops and keep only one-way link
sif<-read.delim(file="Enfermos1_all.sif", sep="\t", header=FALSE)
dim(sif)
# [1] 233294807         3

head(sif)
#     V1         V2      V3
# 1 A1BG 0.01390540  GTPBP6
# 2 A1BG 0.01090800 EFCAB12
# 3 A1BG 0.00660273   GGACT
# 4 A1BG 0.01640860   A2ML1
# 5 A1BG 0.01072390     A2M
# 6 A1BG 0.01407260  A4GALT

names(sif)<-c("From", "Interaction", "To")
sif$From<-as.character(sif$From)
sif$To<-as.character(sif$To)

##Self loops
ids<-sif$From == sif$To
table(ids)
# ids 
#     FALSE      TRUE
# 233294793        14

sif<-sif[!ids, ]
dim(sif)
# [1] 233294793         3

##Remove directional links
sif<-sif[order(sif$Interaction), ]
head(sif)
#              From Interaction      To
# 85169971    HCAR1  0.00243109  ZNF157
# 225431252  ZNF157  0.00243109   HCAR1
# 117714980    MIAT  0.00263478 TMEM38A
# 203140116 TMEM38A  0.00263478    MIAT
# 100906241  KIF18B  0.00274587   STIP1
# 190823309   STIP1  0.00274587  KIF18B
pares<-seq(2, nrow(sif), by=2)
aux<-sif$From[pares]
sif$From[pares]<-sif$To[pares]
sif$To[pares]<-aux
head(sif)
#             From Interaction      To
# 85169971   HCAR1  0.00243109  ZNF157
# 225431252  HCAR1  0.00243109  ZNF157
# 117714980   MIAT  0.00263478 TMEM38A
# 203140116   MIAT  0.00263478 TMEM38A
# 100906241 KIF18B  0.00274587   STIP1
# 190823309 KIF18B  0.00274587   STIP1
sif<-unique(sif)
dim(sif)
# [1] 174692160         3

head(sif)
#              From Interaction      To
# 85169971    HCAR1  0.00243109  ZNF157
# 117714980    MIAT  0.00263478 TMEM38A
# 100906241  KIF18B  0.00274587   STIP1
# 37736698      CGN  0.00279620   PSMC1
# 119102031   MMP17  0.00280426   PNMA2
# 69342486  FASTKD1  0.00281021  HILPDA

233294793/2
# [1] 116647396 

##OJO no es la mitad ¿Por que?
table(duplicated(sif$Interaction))
#     FALSE      TRUE 
#   1071587 173620573 
  
ids<-which(duplicated(sif$Interaction))
aux<-sif[ids[1]+(-2:2),]
aux
#              From Interaction      To
# 8214061   ANGPTL4  0.00318407   NELL2
# 42571529  COL11A1  0.00318540   NRDE2
# 67547253  PLEKHO2  0.00318540 FAM221A
# 132526313   NRDE2  0.00318540 COL11A1
# 147209353 FAM221A  0.00318540 PLEKHO2

unicos<-sif[!sif$Interaction%in%sif$Interaction[ids],]
duplicados<-sif[sif$Interaction%in%sif$Interaction[ids],]

rm("sif")

head(duplicados)
#              From Interaction      To
# 42571529  COL11A1  0.00318540   NRDE2
# 67547253  PLEKHO2  0.00318540 FAM221A
# 132526313   NRDE2  0.00318540 COL11A1
# 147209353 FAM221A  0.00318540 PLEKHO2
# 7977348   BHLHE22  0.00337844 ANAPC11
# 112062159   USP30  0.00337844   MAGI3

ids<-unique(duplicados$Interaction)
library("parallel")
intercambiar<-mclapply(ids, function(Interaction){
    aux<-which(duplicados$Interaction==Interaction)
    aux[1:(length(aux)%/%2)]
}, mc.cores=10)   
    
aux<-duplicados$From[intercambiar]
duplicados$From[intercambiar]<-duplicados$To[intercambiar]
duplicados$To[intercambiar]<-aux
duplicados<-unique(duplicados)



  
  