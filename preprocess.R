
#################################################################################
## Sanos
##  Fecha: 07/07/2015
##  Bajamos todos los datos de Breast Cancer
##  Filtros:
##	-Center/Platform
##		UNC(IlluminaHiSeq_RNASeq)
##	-Tumor/Normal
##		Normal-matched
##  Muestras: 101
#################################################################################
DATADIR <- '/pipeline/data/'
RDATA <- paste(DATADIR, "rdata", sep="")
dir.create(RDATA)
sanosFiles<-dir(full.names=TRUE, pattern = "*.trimmed.annotated.gene.quantification.txt",
                                 path= paste(DATADIR, "counts/healthy", sep=""))
sanos<-lapply(sanosFiles, read.delim, sep="\t")
length(sanos)
# [1] 101

head(sanos[[1]])
#          gene raw_counts median_length_normalized       RPKM
# 1 ?|100130426          2                 0.245098  0.0646161
# 2 ?|100133144        135                 5.572020  1.4731187
# 3 ?|100134869         88                 2.760351  0.7277216
# 4     ?|10357        170                13.364780  3.5234062
# 5     ?|10431       1502                67.147585 17.7092041
# 6    ?|136542          0                 0.000000  0.0000000

##Controlando el tamanio
tamanio<-unique(do.call(rbind,lapply(sanos, dim)))
tamanio
#       [,1] [,2]
# [1,] 20532    4
stopifnot(nrow(tamanio)==1)

## coinciden el orden de los genes
genes<-do.call(cbind,lapply(sanos, function(x)as.character(x[,1])))
genes<-t(unique(t(genes)))
dim(genes)
# [1] 20532     1
stopifnot(dim(genes)==c(tamanio[1,1], 1))

##Nos quedamos con las raw_counts
sanos<-do.call(cbind, lapply(sanos, function(x)x[,"raw_counts"]))
targets<-data.frame(File=sanosFiles, ID=paste("S", 1:length(sanosFiles), sep=""))
colnames(sanos)<-targets$ID

##Cambiamos la anotacion
genes<-do.call(rbind,sapply(genes[,1], strsplit, split="|", fixed=TRUE))
colnames(genes)<-c("Symbol", "EntrezID")

sanos<-list(Counts=sanos, Annot=genes, targets=targets)
save(sanos,  file=paste(RDATA, "sanosRaw.RData", sep="/"), compress="xz")

#################################################################################
## Enfermos
##  Fecha: 07/07/2015
##  Bajamos todos los datos de Breast Cancer
##  Filtros:
##	-Center/Platform
##		UNC(IlluminaHiSeq_RNASeq)
##	-Tumor/Normal
##		Tumor-matched
##  Muestras: 780
#################################################################################

DATADIR <- '/pipeline/data/'
enfermosFiles<-dir(full.names=TRUE, pattern = "*.trimmed.annotated.gene.quantification.txt",
                                 path= paste(DATADIR, "counts/cancer", sep=""))
enfermos<-lapply(enfermosFiles, read.delim, sep="\t")
length(enfermos)
# [1] 780

head(enfermos[[1]])
#          gene raw_counts median_length_normalized        RPKM
# 1 ?|100130426          0                0.0000000  0.00000000
# 2 ?|100133144         74                3.0604305  0.61684584
# 3 ?|100134869         14                0.4391468  0.08844058
# 4     ?|10357        176               13.8364780  2.78655364
# 5     ?|10431       3695              165.2298748 33.28009015
# 6    ?|136542          0                0.0000000  0.00000000

##Controlando el tamanio
tamanio<-unique(do.call(rbind,lapply(enfermos, dim)))
tamanio
#       [,1] [,2]
# [1,] 20532    4
stopifnot(nrow(tamanio)==1)

## coinciden el orden de los genes
genes<-do.call(cbind,lapply(enfermos, function(x)as.character(x[,1])))
genes<-t(unique(t(genes)))
dim(genes)
# [1] 20532     1
stopifnot(dim(genes)==c(tamanio[1,1], 1))

##Nos quedamos con las raw_counts
enfermos<-do.call(cbind, lapply(enfermos, function(x)x[,"raw_counts"]))
targets<-data.frame(File=enfermosFiles, ID=paste("E", 1:length(enfermosFiles), sep=""))
colnames(enfermos)<-targets$ID

##Cambiamos la anotacion
genes<-do.call(rbind,sapply(genes[,1], strsplit, split="|", fixed=TRUE))
colnames(genes)<-c("Symbol", "EntrezID")

enfermos<-list(Counts=enfermos, Annot=genes, targets=targets)
save(enfermos,  file=paste(RDATA, "enfermosRaw.RData", sep="/"), compress="xz")

#################################################################################
## Notacion
##	Database: Ensembl Genes 80
##	Dataset: Homo sapiens genes (GRCh38.p2)
##	Attributes:
##		(NO)-Ensemble Gene ID
##		(NO)-Ensemble Transcript ID
##		-Chromosome Name
##		-Gene start
##		-Gene end
##		-%GC content
##		-Gene type
##		-Entrez Gene ID
##		-HGNC Symbol
##		-HGNC ID(s)
##	Unique Results
#################################################################################
# options(width=150)

# annot<-read.delim(file="mart_export.txt", sep="\t")
annot<-read.delim(file="EntreZ_mart_export.txt", sep="\t")
dim(annot)
# [1] 218384      10
# [1] 66752     8

names(annot)<-c("E.Gene.ID", "E.Transcript.ID", "Chr", "Start", "End", "GC", "Type", 
    "EntrezID", "HGNCID", "Symbol")
names(annot)<-c("Chr", "Start", "End", "GC", "Type", "EntrezID", "HGNCID", "Symbol")
    
annot$Length<-abs(annot$End - annot$Start)
head(annot)
#   Chr Start  End    GC           Type EntrezID    HGNCID  Symbol Length
# 1  MT   577  647 40.85        Mt_tRNA       NA HGNC:7481   MT-TF     70
# 2  MT   648 1601 45.49        Mt_rRNA       NA HGNC:7470 MT-RNR1    953
# 3  MT  1602 1670 42.03        Mt_tRNA       NA HGNC:7500   MT-TV     68
# 4  MT  1671 3229 42.81        Mt_rRNA       NA HGNC:7471 MT-RNR2   1558
# 5  MT  3230 3304 38.67        Mt_tRNA       NA HGNC:7490  MT-TL1     74
# 6  MT  3307 4262 47.70 protein_coding     4535 HGNC:7455  MT-ND1    955

table(annot$Type)
#           3prime_overlapping_ncrna                          antisense                          IG_C_gene 
#                                 29                               5660                                 14 
#                    IG_C_pseudogene                          IG_D_gene                          IG_J_gene 
#                                  9                                 37                                 18 
#                    IG_J_pseudogene                          IG_V_gene                    IG_V_pseudogene 
#                                  3                                154                                190 
#                            lincRNA                           LRG_gene                       macro_lncRNA 
#                               7880                                562                                  1 
#                              miRNA                           misc_RNA                            Mt_rRNA 
#                               4962                               2468                                  2 
#                            Mt_tRNA                         non_coding             polymorphic_pseudogene 
#                                 22                                  3                                 79 
#               processed_pseudogene               processed_transcript                     protein_coding 
#                              10704                                770                              22466 
#                         pseudogene                           ribozyme                               rRNA 
#                                 43                                  8                                564 
#                             scaRNA                     sense_intronic                  sense_overlapping 
#                                 51                                932                                199 
#                             snoRNA                              snRNA                               sRNA 
#                               1017                               2041                                 20 
#                                TEC   transcribed_processed_pseudogene     transcribed_unitary_pseudogene 
#                               1057                                481                                  1 
# transcribed_unprocessed_pseudogene    translated_processed_pseudogene  translated_unprocessed_pseudogene 
#                                768                                  1                                  1 
#                          TR_C_gene                          TR_D_gene                          TR_J_gene 
#                                  5                                  3                                 73 
#                    TR_J_pseudogene                          TR_V_gene                    TR_V_pseudogene 
#                                  4                                106                                 30 
#                 unitary_pseudogene             unprocessed_pseudogene                           vaultRNA 
#                                170                               3143                                  1                        1
save(annot, file="EntreZ_annot.RData", compress="xz")

#################################################################################
##Merging 
##	M=sanos|enfermos
##	Annot=genes|biomart
##	targets=sanos+enfermos
#################################################################################

options(width=130)
# load(file="annot.RData")
load(file="EntreZ_annot.RData")
load(file="enfermosRaw.RData")
load(file="sanosRaw.RData")

M<-cbind(sanos$Counts, enfermos$Counts)
dim(M)
# [1] 20532   881

targets<-rbind(sanos$targets, enfermos$targets)
targets<-data.frame(targets, stringsAsFactors=FALSE)
targets$File<-as.character(targets$File)
targets$ID<-as.character(targets$ID)
dim(targets)
# [1] 881   2

##check M y targets
stopifnot(nrow(targets)==ncol(M))

##check Annot
##igual sanos y enfermos?
stopifnot(all(row.names(sanos$Annot)==row.names(enfermos$Annot)))
Annot<-sanos$Annot
##Hay EntrezID repetidos?
stopifnot(length(unique(as.numeric(Annot[, "EntrezID"])))==nrow(Annot))
Annot<-data.frame(Annot, stringsAsFactors=FALSE)
Annot$EntrezID<-as.integer(Annot$EntrezID)
Annot$Row<-1:nrow(Annot)

##Eliminando los chromosomas no convenvionales
annot<-as.data.frame(annot)
levels(annot$Chr)
{#   [1] "1"                                      "10"                                     "11"                                    
#   [4] "12"                                     "13"                                     "14"                                    
#   [7] "15"                                     "16"                                     "17"                                    
#  [10] "18"                                     "19"                                     "2"                                     
#  [13] "20"                                     "21"                                     "22"                                    
#  [16] "3"                                      "4"                                      "5"                                     
#  [19] "6"                                      "7"                                      "8"                                     
#  [22] "9"                                      "CHR_HG126_PATCH"                        "CHR_HG1362_PATCH"                      
#  [25] "CHR_HG142_HG150_NOVEL_TEST"             "CHR_HG151_NOVEL_TEST"                   "CHR_HG1832_PATCH"                      
#  [28] "CHR_HG2021_PATCH"                       "CHR_HG2030_PATCH"                       "CHR_HG2058_PATCH"                      
#  [31] "CHR_HG2062_PATCH"                       "CHR_HG2066_PATCH"                       "CHR_HG2095_PATCH"                      
#  [34] "CHR_HG2104_PATCH"                       "CHR_HG2128_PATCH"                       "CHR_HG2191_PATCH"                      
#  [37] "CHR_HG2217_PATCH"                       "CHR_HG2232_PATCH"                       "CHR_HG2233_PATCH"                      
#  [40] "CHR_HG2247_PATCH"                       "CHR_HG2249_PATCH"                       "CHR_HG2288_HG2289_PATCH"               
#  [43] "CHR_HG2291_PATCH"                       "CHR_HG986_PATCH"                        "CHR_HSCHR10_1_CTG1"                    
#  [46] "CHR_HSCHR10_1_CTG2"                     "CHR_HSCHR10_1_CTG4"                     "CHR_HSCHR11_1_CTG1_2"                  
#  [49] "CHR_HSCHR11_1_CTG5"                     "CHR_HSCHR11_1_CTG6"                     "CHR_HSCHR11_1_CTG7"                    
#  [52] "CHR_HSCHR11_1_CTG8"                     "CHR_HSCHR11_2_CTG1"                     "CHR_HSCHR11_2_CTG1_1"                  
#  [55] "CHR_HSCHR11_3_CTG1"                     "CHR_HSCHR1_1_CTG11"                     "CHR_HSCHR1_1_CTG3"                     
#  [58] "CHR_HSCHR1_1_CTG31"                     "CHR_HSCHR1_1_CTG32_1"                   "CHR_HSCHR12_1_CTG1"                    
#  [61] "CHR_HSCHR12_1_CTG2_1"                   "CHR_HSCHR12_2_CTG2"                     "CHR_HSCHR12_2_CTG2_1"                  
#  [64] "CHR_HSCHR12_3_CTG2"                     "CHR_HSCHR12_3_CTG2_1"                   "CHR_HSCHR12_4_CTG2"                    
#  [67] "CHR_HSCHR12_4_CTG2_1"                   "CHR_HSCHR12_5_CTG2"                     "CHR_HSCHR12_5_CTG2_1"                  
#  [70] "CHR_HSCHR12_6_CTG2_1"                   "CHR_HSCHR1_2_CTG3"                      "CHR_HSCHR1_2_CTG31"                    
#  [73] "CHR_HSCHR1_2_CTG32_1"                   "CHR_HSCHR13_1_CTG1"                     "CHR_HSCHR13_1_CTG3"                    
#  [76] "CHR_HSCHR1_3_CTG31"                     "CHR_HSCHR1_3_CTG32_1"                   "CHR_HSCHR14_1_CTG1"                    
#  [79] "CHR_HSCHR14_2_CTG1"                     "CHR_HSCHR14_3_CTG1"                     "CHR_HSCHR14_7_CTG1"                    
#  [82] "CHR_HSCHR1_4_CTG31"                     "CHR_HSCHR15_1_CTG1"                     "CHR_HSCHR15_1_CTG3"                    
#  [85] "CHR_HSCHR15_1_CTG8"                     "CHR_HSCHR15_2_CTG3"                     "CHR_HSCHR15_2_CTG8"                    
#  [88] "CHR_HSCHR15_3_CTG3"                     "CHR_HSCHR15_3_CTG8"                     "CHR_HSCHR15_4_CTG8"                    
#  [91] "CHR_HSCHR15_5_CTG8"                     "CHR_HSCHR16_1_CTG1"                     "CHR_HSCHR16_1_CTG3_1"                  
#  [94] "CHR_HSCHR16_2_CTG3_1"                   "CHR_HSCHR16_3_CTG1"                     "CHR_HSCHR16_4_CTG1"                    
#  [97] "CHR_HSCHR16_CTG2"                       "CHR_HSCHR17_10_CTG4"                    "CHR_HSCHR17_1_CTG1"                    
# [100] "CHR_HSCHR17_1_CTG2"                     "CHR_HSCHR17_1_CTG4"                     "CHR_HSCHR17_1_CTG5"                    
# [103] "CHR_HSCHR17_1_CTG9"                     "CHR_HSCHR17_2_CTG1"                     "CHR_HSCHR17_2_CTG2"                    
# [106] "CHR_HSCHR17_2_CTG5"                     "CHR_HSCHR17_3_CTG2"                     "CHR_HSCHR17_4_CTG4"                    
# [109] "CHR_HSCHR17_5_CTG4"                     "CHR_HSCHR17_6_CTG4"                     "CHR_HSCHR17_7_CTG4"                    
# [112] "CHR_HSCHR17_8_CTG4"                     "CHR_HSCHR18_1_CTG1_1"                   "CHR_HSCHR18_1_CTG2_1"                  
# [115] "CHR_HSCHR18_2_CTG2"                     "CHR_HSCHR18_2_CTG2_1"                   "CHR_HSCHR18_ALT21_CTG2_1"              
# [118] "CHR_HSCHR18_ALT2_CTG2_1"                "CHR_HSCHR19_1_CTG2"                     "CHR_HSCHR19_1_CTG3_1"                  
# [121] "CHR_HSCHR19_2_CTG2"                     "CHR_HSCHR19_3_CTG2"                     "CHR_HSCHR19_3_CTG3_1"                  
# [124] "CHR_HSCHR19_4_CTG2"                     "CHR_HSCHR19_4_CTG3_1"                   "CHR_HSCHR19_5_CTG2"                    
# [127] "CHR_HSCHR19KIR_ABC08_A1_HAP_CTG3_1"     "CHR_HSCHR19KIR_ABC08_AB_HAP_C_P_CTG3_1" "CHR_HSCHR19KIR_ABC08_AB_HAP_T_P_CTG3_1"
# [130] "CHR_HSCHR19KIR_FH05_A_HAP_CTG3_1"       "CHR_HSCHR19KIR_FH05_B_HAP_CTG3_1"       "CHR_HSCHR19KIR_FH06_A_HAP_CTG3_1"      
# [133] "CHR_HSCHR19KIR_FH06_BA1_HAP_CTG3_1"     "CHR_HSCHR19KIR_FH08_A_HAP_CTG3_1"       "CHR_HSCHR19KIR_FH08_BAX_HAP_CTG3_1"    
# [136] "CHR_HSCHR19KIR_FH13_A_HAP_CTG3_1"       "CHR_HSCHR19KIR_FH13_BA2_HAP_CTG3_1"     "CHR_HSCHR19KIR_FH15_A_HAP_CTG3_1"      
# [139] "CHR_HSCHR19KIR_FH15_B_HAP_CTG3_1"       "CHR_HSCHR19KIR_G085_A_HAP_CTG3_1"       "CHR_HSCHR19KIR_G085_BA1_HAP_CTG3_1"    
# [142] "CHR_HSCHR19KIR_G248_A_HAP_CTG3_1"       "CHR_HSCHR19KIR_G248_BA2_HAP_CTG3_1"     "CHR_HSCHR19KIR_GRC212_AB_HAP_CTG3_1"   
# [145] "CHR_HSCHR19KIR_GRC212_BA1_HAP_CTG3_1"   "CHR_HSCHR19KIR_LUCE_A_HAP_CTG3_1"       "CHR_HSCHR19KIR_LUCE_BDEL_HAP_CTG3_1"   
# [148] "CHR_HSCHR19KIR_RP5_B_HAP_CTG3_1"        "CHR_HSCHR19KIR_RSH_A_HAP_CTG3_1"        "CHR_HSCHR19KIR_RSH_BA2_HAP_CTG3_1"     
# [151] "CHR_HSCHR19KIR_T7526_A_HAP_CTG3_1"      "CHR_HSCHR19KIR_T7526_BDEL_HAP_CTG3_1"   "CHR_HSCHR19LRC_COX1_CTG3_1"            
# [154] "CHR_HSCHR19LRC_COX2_CTG3_1"             "CHR_HSCHR19LRC_LRC_I_CTG3_1"            "CHR_HSCHR19LRC_LRC_J_CTG3_1"           
# [157] "CHR_HSCHR19LRC_LRC_S_CTG3_1"            "CHR_HSCHR19LRC_LRC_T_CTG3_1"            "CHR_HSCHR19LRC_PGF1_CTG3_1"            
# [160] "CHR_HSCHR19LRC_PGF2_CTG3_1"             "CHR_HSCHR1_ALT2_1_CTG32_1"              "CHR_HSCHR20_1_CTG1"                    
# [163] "CHR_HSCHR20_1_CTG2"                     "CHR_HSCHR20_1_CTG3"                     "CHR_HSCHR20_1_CTG4"                    
# [166] "CHR_HSCHR21_2_CTG1_1"                   "CHR_HSCHR21_3_CTG1_1"                   "CHR_HSCHR21_4_CTG1_1"                  
# [169] "CHR_HSCHR21_5_CTG2"                     "CHR_HSCHR21_6_CTG1_1"                   "CHR_HSCHR21_8_CTG1_1"                  
# [172] "CHR_HSCHR2_1_CTG1"                      "CHR_HSCHR2_1_CTG15"                     "CHR_HSCHR2_1_CTG5"                     
# [175] "CHR_HSCHR2_1_CTG7"                      "CHR_HSCHR2_1_CTG7_2"                    "CHR_HSCHR22_1_CTG1"                    
# [178] "CHR_HSCHR22_1_CTG2"                     "CHR_HSCHR22_1_CTG3"                     "CHR_HSCHR22_1_CTG4"                    
# [181] "CHR_HSCHR22_1_CTG5"                     "CHR_HSCHR22_1_CTG6"                     "CHR_HSCHR22_1_CTG7"                    
# [184] "CHR_HSCHR22_2_CTG1"                     "CHR_HSCHR22_3_CTG1"                     "CHR_HSCHR22_4_CTG1"                    
# [187] "CHR_HSCHR22_5_CTG1"                     "CHR_HSCHR2_2_CTG1"                      "CHR_HSCHR2_2_CTG15"                    
# [190] "CHR_HSCHR2_2_CTG7"                      "CHR_HSCHR2_2_CTG7_2"                    "CHR_HSCHR2_3_CTG1"                     
# [193] "CHR_HSCHR2_3_CTG15"                     "CHR_HSCHR2_3_CTG7_2"                    "CHR_HSCHR2_4_CTG1"                     
# [196] "CHR_HSCHR3_1_CTG1"                      "CHR_HSCHR3_1_CTG2_1"                    "CHR_HSCHR3_1_CTG3"                     
# [199] "CHR_HSCHR3_2_CTG2_1"                    "CHR_HSCHR3_2_CTG3"                      "CHR_HSCHR3_3_CTG1"                     
# [202] "CHR_HSCHR3_3_CTG3"                      "CHR_HSCHR3_4_CTG2_1"                    "CHR_HSCHR3_4_CTG3"                     
# [205] "CHR_HSCHR3_5_CTG2_1"                    "CHR_HSCHR3_5_CTG3"                      "CHR_HSCHR3_6_CTG3"                     
# [208] "CHR_HSCHR3_7_CTG3"                      "CHR_HSCHR3_8_CTG3"                      "CHR_HSCHR3_9_CTG3"                     
# [211] "CHR_HSCHR4_1_CTG12"                     "CHR_HSCHR4_1_CTG4"                      "CHR_HSCHR4_1_CTG6"                     
# [214] "CHR_HSCHR4_1_CTG9"                      "CHR_HSCHR4_2_CTG12"                     "CHR_HSCHR4_3_CTG12"                    
# [217] "CHR_HSCHR4_4_CTG12"                     "CHR_HSCHR4_5_CTG12"                     "CHR_HSCHR4_6_CTG12"                    
# [220] "CHR_HSCHR4_7_CTG12"                     "CHR_HSCHR5_1_CTG1"                      "CHR_HSCHR5_1_CTG1_1"                   
# [223] "CHR_HSCHR5_1_CTG5"                      "CHR_HSCHR5_2_CTG1"                      "CHR_HSCHR5_2_CTG1_1"                   
# [226] "CHR_HSCHR5_2_CTG5"                      "CHR_HSCHR5_3_CTG1"                      "CHR_HSCHR5_3_CTG5"                     
# [229] "CHR_HSCHR5_4_CTG1"                      "CHR_HSCHR5_4_CTG1_1"                    "CHR_HSCHR5_5_CTG1"                     
# [232] "CHR_HSCHR5_6_CTG1"                      "CHR_HSCHR5_7_CTG1"                      "CHR_HSCHR6_1_CTG2"                     
# [235] "CHR_HSCHR6_1_CTG3"                      "CHR_HSCHR6_1_CTG4"                      "CHR_HSCHR6_1_CTG5"                     
# [238] "CHR_HSCHR6_1_CTG6"                      "CHR_HSCHR6_1_CTG7"                      "CHR_HSCHR6_1_CTG8"                     
# [241] "CHR_HSCHR6_1_CTG9"                      "CHR_HSCHR6_8_CTG1"                      "CHR_HSCHR6_MHC_APD_CTG1"               
# [244] "CHR_HSCHR6_MHC_COX_CTG1"                "CHR_HSCHR6_MHC_DBB_CTG1"                "CHR_HSCHR6_MHC_MANN_CTG1"              
# [247] "CHR_HSCHR6_MHC_MCF_CTG1"                "CHR_HSCHR6_MHC_QBL_CTG1"                "CHR_HSCHR6_MHC_SSTO_CTG1"              
# [250] "CHR_HSCHR7_1_CTG1"                      "CHR_HSCHR7_1_CTG4_4"                    "CHR_HSCHR7_1_CTG6"                     
# [253] "CHR_HSCHR7_2_CTG4_4"                    "CHR_HSCHR7_2_CTG6"                      "CHR_HSCHR7_3_CTG6"                     
# [256] "CHR_HSCHR8_1_CTG1"                      "CHR_HSCHR8_1_CTG6"                      "CHR_HSCHR8_1_CTG7"                     
# [259] "CHR_HSCHR8_2_CTG7"                      "CHR_HSCHR8_3_CTG1"                      "CHR_HSCHR8_3_CTG7"                     
# [262] "CHR_HSCHR8_4_CTG7"                      "CHR_HSCHR8_5_CTG1"                      "CHR_HSCHR8_5_CTG7"                     
# [265] "CHR_HSCHR8_7_CTG1"                      "CHR_HSCHR8_8_CTG1"                      "CHR_HSCHR8_9_CTG1"                     
# [268] "CHR_HSCHR9_1_CTG1"                      "CHR_HSCHR9_1_CTG2"                      "CHR_HSCHR9_1_CTG3"                     
# [271] "CHR_HSCHR9_1_CTG4"                      "CHR_HSCHR9_1_CTG5"                      "CHR_HSCHRX_1_CTG3"                     
# [274] "CHR_HSCHRX_2_CTG12"                     "CHR_HSCHRX_2_CTG3"                      "GL000008.2"                            
# [277] "GL000009.2"                             "GL000194.1"                             "GL000195.1"                            
# [280] "GL000205.2"                             "GL000213.1"                             "GL000216.2"                            
# [283] "GL000218.1"                             "GL000219.1"                             "GL000220.1"                            
# [286] "GL000224.1"                             "GL000225.1"                             "KI270442.1"                            
# [289] "KI270706.1"                             "KI270707.1"                             "KI270708.1"                            
# [292] "KI270711.1"                             "KI270713.1"                             "KI270714.1"                            
# [295] "KI270721.1"                             "KI270722.1"                             "KI270723.1"                            
# [298] "KI270724.1"                             "KI270726.1"                             "KI270727.1"                            
# [301] "KI270728.1"                             "KI270731.1"                             "KI270733.1"                            
# [304] "KI270734.1"                             "KI270741.1"                             "KI270743.1"                            
# [307] "KI270744.1"                             "KI270750.1"                             "KI270752.1"                            
# [310] "LRG_1"                                  "LRG_10"                                 "LRG_100"                               
# [313] "LRG_101"                                "LRG_102"                                "LRG_103"                               
# [316] "LRG_104"                                "LRG_105"                                "LRG_106"                               
# [319] "LRG_107"                                "LRG_108"                                "LRG_109"                               
# [322] "LRG_11"                                 "LRG_110"                                "LRG_111"                               
# [325] "LRG_112"                                "LRG_113"                                "LRG_114"                               
# [328] "LRG_115"                                "LRG_116"                                "LRG_117"                               
# [331] "LRG_118"                                "LRG_119"                                "LRG_12"                                
# [334] "LRG_120"                                "LRG_121"                                "LRG_122"                               
# [337] "LRG_123"                                "LRG_124"                                "LRG_125"                               
# [340] "LRG_126"                                "LRG_127"                                "LRG_128"                               
# [343] "LRG_129"                                "LRG_13"                                 "LRG_130"                               
# [346] "LRG_132"                                "LRG_133"                                "LRG_134"                               
# [349] "LRG_135"                                "LRG_136"                                "LRG_137"                               
# [352] "LRG_138"                                "LRG_139"                                "LRG_140"                               
# [355] "LRG_141"                                "LRG_142"                                "LRG_144"                               
# [358] "LRG_145"                                "LRG_146"                                "LRG_147"                               
# [361] "LRG_148"                                "LRG_149"                                "LRG_15"                                
# [364] "LRG_150"                                "LRG_151"                                "LRG_152"                               
# [367] "LRG_154"                                "LRG_155"                                "LRG_156"                               
# [370] "LRG_157"                                "LRG_158"                                "LRG_159"                               
# [373] "LRG_16"                                 "LRG_160"                                "LRG_161"                               
# [376] "LRG_162"                                "LRG_163"                                "LRG_164"                               
# [379] "LRG_165"                                "LRG_168"                                "LRG_169"                               
# [382] "LRG_17"                                 "LRG_170"                                "LRG_171"                               
# [385] "LRG_172"                                "LRG_173"                                "LRG_174"                               
# [388] "LRG_175"                                "LRG_176"                                "LRG_177"                               
# [391] "LRG_178"                                "LRG_179"                                "LRG_18"                                
# [394] "LRG_180"                                "LRG_182"                                "LRG_183"                               
# [397] "LRG_184"                                "LRG_185"                                "LRG_186"                               
# [400] "LRG_187"                                "LRG_188"                                "LRG_189"                               
# [403] "LRG_19"                                 "LRG_190"                                "LRG_191"                               
# [406] "LRG_192"                                "LRG_193"                                "LRG_194"                               
# [409] "LRG_195"                                "LRG_196"                                "LRG_197"                               
# [412] "LRG_198"                                "LRG_199"                                "LRG_2"                                 
# [415] "LRG_20"                                 "LRG_200"                                "LRG_201"                               
# [418] "LRG_202"                                "LRG_203"                                "LRG_204"                               
# [421] "LRG_205"                                "LRG_207"                                "LRG_208"                               
# [424] "LRG_209"                                "LRG_21"                                 "LRG_210"                               
# [427] "LRG_211"                                "LRG_212"                                "LRG_213"                               
# [430] "LRG_214"                                "LRG_215"                                "LRG_216"                               
# [433] "LRG_217"                                "LRG_218"                                "LRG_219"                               
# [436] "LRG_22"                                 "LRG_220"                                "LRG_221"                               
# [439] "LRG_226"                                "LRG_227"                                "LRG_228"                               
# [442] "LRG_229"                                "LRG_23"                                 "LRG_230"                               
# [445] "LRG_231"                                "LRG_234"                                "LRG_236"                               
# [448] "LRG_239"                                "LRG_24"                                 "LRG_241"                               
# [451] "LRG_242"                                "LRG_243"                                "LRG_245"                               
# [454] "LRG_246"                                "LRG_248"                                "LRG_249"                               
# [457] "LRG_25"                                 "LRG_250"                                "LRG_251"                               
# [460] "LRG_252"                                "LRG_253"                                "LRG_254"                               
# [463] "LRG_255"                                "LRG_256"                                "LRG_257"                               
# [466] "LRG_258"                                "LRG_26"                                 "LRG_260"                               
# [469] "LRG_261"                                "LRG_262"                                "LRG_263"                               
# [472] "LRG_264"                                "LRG_265"                                "LRG_266"                               
# [475] "LRG_267"                                "LRG_268"                                "LRG_269"                               
# [478] "LRG_27"                                 "LRG_270"                                "LRG_271"                               
# [481] "LRG_272"                                "LRG_273"                                "LRG_274"                               
# [484] "LRG_275"                                "LRG_276"                                "LRG_278"                               
# [487] "LRG_279"                                "LRG_28"                                 "LRG_280"                               
# [490] "LRG_281"                                "LRG_283"                                "LRG_284"                               
# [493] "LRG_285"                                "LRG_286"                                "LRG_287"                               
# [496] "LRG_288"                                "LRG_289"                                "LRG_29"                                
# [499] "LRG_290"                                "LRG_291"                                "LRG_292"                               
# [502] "LRG_293"                                "LRG_294"                                "LRG_295"                               
# [505] "LRG_296"                                "LRG_298"                                "LRG_299"                               
# [508] "LRG_3"                                  "LRG_30"                                 "LRG_300"                               
# [511] "LRG_301"                                "LRG_304"                                "LRG_306"                               
# [514] "LRG_307"                                "LRG_308"                                "LRG_309"                               
# [517] "LRG_31"                                 "LRG_310"                                "LRG_311"                               
# [520] "LRG_314"                                "LRG_316"                                "LRG_317"                               
# [523] "LRG_318"                                "LRG_319"                                "LRG_32"                                
# [526] "LRG_321"                                "LRG_322"                                "LRG_325"                               
# [529] "LRG_326"                                "LRG_327"                                "LRG_328"                               
# [532] "LRG_329"                                "LRG_33"                                 "LRG_330"                               
# [535] "LRG_331"                                "LRG_332"                                "LRG_333"                               
# [538] "LRG_334"                                "LRG_335"                                "LRG_336"                               
# [541] "LRG_337"                                "LRG_34"                                 "LRG_340"                               
# [544] "LRG_341"                                "LRG_343"                                "LRG_345"                               
# [547] "LRG_346"                                "LRG_347"                                "LRG_348"                               
# [550] "LRG_349"                                "LRG_35"                                 "LRG_350"                               
# [553] "LRG_351"                                "LRG_352"                                "LRG_353"                               
# [556] "LRG_354"                                "LRG_355"                                "LRG_356"                               
# [559] "LRG_357"                                "LRG_358"                                "LRG_359"                               
# [562] "LRG_36"                                 "LRG_361"                                "LRG_362"                               
# [565] "LRG_363"                                "LRG_364"                                "LRG_365"                               
# [568] "LRG_366"                                "LRG_368"                                "LRG_369"                               
# [571] "LRG_37"                                 "LRG_371"                                "LRG_372"                               
# [574] "LRG_373"                                "LRG_374"                                "LRG_375"                               
# [577] "LRG_376"                                "LRG_377"                                "LRG_378"                               
# [580] "LRG_379"                                "LRG_38"                                 "LRG_380"                               
# [583] "LRG_382"                                "LRG_383"                                "LRG_384"                               
# [586] "LRG_385"                                "LRG_386"                                "LRG_388"                               
# [589] "LRG_389"                                "LRG_39"                                 "LRG_390"                               
# [592] "LRG_391"                                "LRG_392"                                "LRG_393"                               
# [595] "LRG_394"                                "LRG_395"                                "LRG_396"                               
# [598] "LRG_397"                                "LRG_398"                                "LRG_399"                               
# [601] "LRG_4"                                  "LRG_40"                                 "LRG_400"                               
# [604] "LRG_401"                                "LRG_403"                                "LRG_404"                               
# [607] "LRG_405"                                "LRG_406"                                "LRG_408"                               
# [610] "LRG_409"                                "LRG_41"                                 "LRG_410"                               
# [613] "LRG_411"                                "LRG_413"                                "LRG_414"                               
# [616] "LRG_415"                                "LRG_416"                                "LRG_417"                               
# [619] "LRG_419"                                "LRG_42"                                 "LRG_421"                               
# [622] "LRG_422"                                "LRG_423"                                "LRG_424"                               
# [625] "LRG_426"                                "LRG_43"                                 "LRG_432"                               
# [628] "LRG_433"                                "LRG_434"                                "LRG_435"                               
# [631] "LRG_437"                                "LRG_439"                                "LRG_44"                                
# [634] "LRG_440"                                "LRG_442"                                "LRG_443"                               
# [637] "LRG_444"                                "LRG_445"                                "LRG_446"                               
# [640] "LRG_447"                                "LRG_448"                                "LRG_449"                               
# [643] "LRG_45"                                 "LRG_450"                                "LRG_451"                               
# [646] "LRG_452"                                "LRG_454"                                "LRG_455"                               
# [649] "LRG_456"                                "LRG_457"                                "LRG_458"                               
# [652] "LRG_46"                                 "LRG_460"                                "LRG_461"                               
# [655] "LRG_462"                                "LRG_463"                                "LRG_464"                               
# [658] "LRG_465"                                "LRG_466"                                "LRG_467"                               
# [661] "LRG_469"                                "LRG_47"                                 "LRG_470"                               
# [664] "LRG_471"                                "LRG_472"                                "LRG_473"                               
# [667] "LRG_474"                                "LRG_475"                                "LRG_476"                               
# [670] "LRG_48"                                 "LRG_482"                                "LRG_483"                               
# [673] "LRG_484"                                "LRG_485"                                "LRG_486"                               
# [676] "LRG_487"                                "LRG_488"                                "LRG_489"                               
# [679] "LRG_49"                                 "LRG_490"                                "LRG_491"                               
# [682] "LRG_492"                                "LRG_493"                                "LRG_494"                               
# [685] "LRG_495"                                "LRG_496"                                "LRG_497"                               
# [688] "LRG_498"                                "LRG_499"                                "LRG_5"                                 
# [691] "LRG_50"                                 "LRG_500"                                "LRG_501"                               
# [694] "LRG_502"                                "LRG_503"                                "LRG_504"                               
# [697] "LRG_505"                                "LRG_507"                                "LRG_509"                               
# [700] "LRG_51"                                 "LRG_510"                                "LRG_511"                               
# [703] "LRG_512"                                "LRG_513"                                "LRG_514"                               
# [706] "LRG_515"                                "LRG_516"                                "LRG_517"                               
# [709] "LRG_518"                                "LRG_519"                                "LRG_52"                                
# [712] "LRG_520"                                "LRG_521"                                "LRG_522"                               
# [715] "LRG_523"                                "LRG_524"                                "LRG_525"                               
# [718] "LRG_526"                                "LRG_527"                                "LRG_528"                               
# [721] "LRG_529"                                "LRG_53"                                 "LRG_530"                               
# [724] "LRG_531"                                "LRG_532"                                "LRG_533"                               
# [727] "LRG_534"                                "LRG_535"                                "LRG_536"                               
# [730] "LRG_537"                                "LRG_538"                                "LRG_539"                               
# [733] "LRG_54"                                 "LRG_540"                                "LRG_541"                               
# [736] "LRG_55"                                 "LRG_556"                                "LRG_557"                               
# [739] "LRG_56"                                 "LRG_57"                                 "LRG_58"                                
# [742] "LRG_59"                                 "LRG_6"                                  "LRG_60"                                
# [745] "LRG_607"                                "LRG_608"                                "LRG_609"                               
# [748] "LRG_61"                                 "LRG_610"                                "LRG_611"                               
# [751] "LRG_612"                                "LRG_613"                                "LRG_614"                               
# [754] "LRG_615"                                "LRG_616"                                "LRG_617"                               
# [757] "LRG_618"                                "LRG_62"                                 "LRG_620"                               
# [760] "LRG_621"                                "LRG_622"                                "LRG_623"                               
# [763] "LRG_625"                                "LRG_627"                                "LRG_629"                               
# [766] "LRG_63"                                 "LRG_631"                                "LRG_64"                                
# [769] "LRG_640"                                "LRG_642"                                "LRG_643"                               
# [772] "LRG_65"                                 "LRG_650"                                "LRG_652"                               
# [775] "LRG_653"                                "LRG_657"                                "LRG_659"                               
# [778] "LRG_66"                                 "LRG_661"                                "LRG_662"                               
# [781] "LRG_664"                                "LRG_665"                                "LRG_666"                               
# [784] "LRG_669"                                "LRG_67"                                 "LRG_670"                               
# [787] "LRG_672"                                "LRG_673"                                "LRG_674"                               
# [790] "LRG_675"                                "LRG_676"                                "LRG_683"                               
# [793] "LRG_684"                                "LRG_685"                                "LRG_687"                               
# [796] "LRG_689"                                "LRG_69"                                 "LRG_690"                               
# [799] "LRG_691"                                "LRG_692"                                "LRG_693"                               
# [802] "LRG_697"                                "LRG_7"                                  "LRG_70"                                
# [805] "LRG_700"                                "LRG_702"                                "LRG_71"                                
# [808] "LRG_715"                                "LRG_717"                                "LRG_719"                               
# [811] "LRG_72"                                 "LRG_720"                                "LRG_721"                               
# [814] "LRG_722"                                "LRG_723"                                "LRG_725"                               
# [817] "LRG_726"                                "LRG_727"                                "LRG_73"                                
# [820] "LRG_733"                                "LRG_734"                                "LRG_74"                                
# [823] "LRG_741"                                "LRG_742"                                "LRG_744"                               
# [826] "LRG_748"                                "LRG_749"                                "LRG_75"                                
# [829] "LRG_750"                                "LRG_751"                                "LRG_753"                               
# [832] "LRG_754"                                "LRG_755"                                "LRG_757"                               
# [835] "LRG_759"                                "LRG_76"                                 "LRG_760"                               
# [838] "LRG_767"                                "LRG_77"                                 "LRG_770"                               
# [841] "LRG_771"                                "LRG_776"                                "LRG_777"                               
# [844] "LRG_778"                                "LRG_78"                                 "LRG_787"                               
# [847] "LRG_788"                                "LRG_79"                                 "LRG_8"                                 
# [850] "LRG_80"                                 "LRG_81"                                 "LRG_83"                                
# [853] "LRG_84"                                 "LRG_85"                                 "LRG_86"                                
# [856] "LRG_88"                                 "LRG_89"                                 "LRG_90"                                
# [859] "LRG_91"                                 "LRG_92"                                 "LRG_93"                                
# [862] "LRG_94"                                 "LRG_96"                                 "LRG_97"                                
# [865] "LRG_98"                                 "LRG_99"                                 "MT"                                    
# [868] "X"                                      "Y"              
}
annot<-annot[annot$Chr%in%c(as.character(1:22), "X", "Y"),]
annot$Chr<-droplevels(annot$Chr)
table(annot$Chr)
#    1   10   11   12   13   14   15   16   17   18   19    2   20   21   22    3    4    5    6    7    8    9    X    Y 
# 5454 2326 3395 3060 1383 2292 2297 2727 3139 1213 3012 4177 1494  967 1402 3170 2644 3004 3009 2998 2468 2371 2538  590

##Sacamos los que no tiene EntrezID
annot<-annot[!is.na(annot$EntrezID), ]

##Sacamos los que no tiene Symbol | HGNCID 
annot<-annot[annot$HGNCID!="" | annot$Symbol!="", ]
head(annot)

dim(annot)
# [1] 218384     10 ##Con Ensembl Gene & Transcript
# [1] 66752     9   ##Sin Ensembl Gene & Transcript
# [1] 61130     9   ##Sin Chr no convencionales
# [1] 25011     9   ##Sin EntrezID
# [1] 23586     9   ##Sin Symbol & HGNCID

Annot<-merge(x=Annot, y=annot, by.x="EntrezID", by.y="EntrezID", 
  all=FALSE, all.x=TRUE, all.y=FALSE, sort=FALSE)
dim(Annot)
# [1] 20657    11 
dim(M)
# [1] 20532   881
nrow(Annot)-nrow(M)
# [1] 125 OJO hay duplicados

##Ojo hay duplicados
Annot<-Annot[order(Annot$Row, Annot$Chr),]
which(duplicated(Annot$Row))
ids<-head(which(duplicated(Annot$Row)))
# [1]  400  415  698 1002 1284 1415
Annot[sort(c(ids-1,ids)),]
#      EntrezID  Symbol.x  Row Chr     Start       End    GC           Type     HGNCID  Symbol.y Length
# 373      2334      AFF2  399   X 148500619 149000663 38.84 protein_coding  HGNC:3776      AFF2 500044
# 374      2334      AFF2  399   X 148501098 148991332 38.79 protein_coding  HGNC:3776      AFF2 490234
# 387    119016     AGAP4  413  10  45825594  45853875 41.11 protein_coding HGNC:23459     AGAP4  28281
# 388    119016     AGAP4  413  10  49982190  50010499 41.10 protein_coding HGNC:23466     AGAP6  28309
# 666    441425 ANKRD20A3  695   9  67859147  67902094 36.24 protein_coding HGNC:23665 ANKRD20A1  42947
# 667    441425 ANKRD20A3  695   9  66106815  66153135 36.11 protein_coding HGNC:31981 ANKRD20A3  46320
# 959     51326    ARL17A  998  17  46516702  46579722 43.35 protein_coding HGNC:24096    ARL17A  63020
# 960     51326    ARL17A  998  17  46274784  46361797 43.65 protein_coding HGNC:32387    ARL17B  87013
# 1238    84126     ATRIP 1279   3  48446710  48467645 48.87 protein_coding HGNC:12269     TREX1  20935
# 1239    84126     ATRIP 1279   3  48446710  48465716 47.47 protein_coding HGNC:33499     ATRIP  19006
# 1666    79741  C10orf68 1720  10  32446140  32574564 37.57 protein_coding HGNC:26533     CCDC7 128424
# 1667    79741  C10orf68 1720  10  32567723  32882874 37.18 protein_coding HGNC:26533     CCDC7 315151

##Que coincidan los symbol
ids<-which(duplicated(Annot$Row))
idsNoCoincidentes<-ids[Annot$Symbol.x[ids]!=Annot$Symbol.y[ids]]
head(Annot[idsNoCoincidentes,])
#      EntrezID Symbol.x  Row Chr    Start      End    GC           Type     HGNCID Symbol.y Length
# 388    119016    AGAP4  413  10 49982190 50010499 41.10 protein_coding HGNC:23466    AGAP6  28309
# 960     51326   ARL17A  998  17 46274784 46361797 43.65 protein_coding HGNC:32387   ARL17B  87013
# 1667    79741 C10orf68 1720  10 32567723 32882874 37.18 protein_coding HGNC:26533    CCDC7 315151
# 2894    55871    CBWD1 3039   9 65668805 65734041 36.05 protein_coding HGNC:24584    CBWD5  65236
# 2895    55871    CBWD1 3039   9 68232003 68300015 36.17 protein_coding HGNC:18519    CBWD3  68012
# 2896   150472    CBWD2 3040   9   121038   179147 36.26 protein_coding HGNC:17134    CBWD1  58109

Annot<-Annot[-idsNoCoincidentes, ]
dim(Annot)
# [1] 20590    11
nrow(Annot)-nrow(M)
# [1] 58 OJO hay duplicados

##No coinciden en el caso anterior tampoco
Annot<-Annot[order(Annot$Row, Annot$Chr),]
ids<-which(duplicated(Annot$Row))
ids<-sort(unique(c(ids-1, ids)))
idsNoCoincidentes<-ids[Annot$Symbol.x[ids]!=Annot$Symbol.y[ids]]
head(Annot[idsNoCoincidentes,])
#      EntrezID  Symbol.x  Row Chr    Start      End    GC           Type     HGNCID  Symbol.y Length
# 666    441425 ANKRD20A3  695   9 67859147 67902094 36.24 protein_coding HGNC:23665 ANKRD20A1  42947
# 1238    84126     ATRIP 1279   3 48446710 48467645 48.87 protein_coding HGNC:12269     TREX1  20935
# 2927   388389   CCDC103 3065  17 44899712 44905390 51.70 protein_coding HGNC:35153   FAM187A   5678
# 3655   548596    CKMT1A 3806  15 43593054 43604901 49.85 protein_coding  HGNC:1995    CKMT1B  11847
# 3657     1159    CKMT1B 3807  15 43692886 43699222 50.24 protein_coding HGNC:31736    CKMT1A   6336
# 3715   283971   CLEC18C 3863  16 69950705 69964452 57.08 protein_coding HGNC:30388   CLEC18A  13747

Annot<-Annot[-idsNoCoincidentes, ]
dim(Annot)
# [1] 20543    11
nrow(Annot)-nrow(M)
# [1] 11 OJO hay duplicados

##Nos quedamos con el de menor GC
Annot<-Annot[order(Annot$Row, Annot$Chr, Annot$GC),]
ids<-which(duplicated(Annot$Row))
Annot<-Annot[-ids,]

dim(Annot)
# [1] 20532    11
nrow(Annot)-nrow(M)
# [1] 0
##Listo el pollo
table(is.na(Annot$GC))
# FALSE  TRUE 
# 19449  1083 

##Estan todas las filas??
stopifnot(all(Annot$Row %in% 1:nrow(M)))

full<-list(M=M, Annot=Annot, Targets=targets) 
save(full, file="RawFull.RData", compress="xz")
