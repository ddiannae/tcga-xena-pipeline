###############################################################################
## CHROMATIN AND GENE REGULATION: FROM GENE TO GENOME FOLDING
## Practical session: In silico analysis of RNA-Seq data.
###############################################################################
## Data preparation:
##  -Get the data
## Author: 
##          Dr. Cristobal Fresno - cristobalfresno@gmail.com
## Date: 2016-12-12
###############################################################################
## Get the data
##      -TCGA Data Portal
##          -Disease: BRCA - Breast invasive carcinoma
##          -Platform: UNC (IlluminaHiSeq_RNASeq)
##          -Data Type: RNASeq
##          -Batch: 61
##          -Data Level: 3 - gene.quantification (only)
##          -Tumour/Normal: 
##              -Tumour-matched -> Tumour.tar.gz  (25 samples)
##              -Normal-matched -> Normal.tar.gz  (12 samples)
##      -Biomart:
##          -Database: Ensembl Genes 80
##          -Dataset: Homo sapiens genes (GRCh38.p2)
##              Attributes:
##                  -Chromosome Name
##                  -Gene start
##                  -Gene end
##                  -%GC content
##                  -Gene type
##                  -Entrez Gene ID
##                  -HGNC ID(s)
##                  -HGNC Symbol
###############################################################################
## 1) Read Normal Data
## 2) Read Tumour Data
## 3) Read Biomart Data
## 4) Merge count and annotation
###############################################################################
##Get the Work and Data dir
###############################################################################
DATADIR <- '/pipeline/data/'
args <- commandArgs(trailingOnly = TRUE)
DATADIR <- args[1]
RDATA <- paste(DATADIR, "rdata", sep="")
dir.create(RDATA)
cat('Data directory: ', DATADIR, '\n')

###############################################################################
##Usefull Libraries
###############################################################################
library("BiocParallel")
library("parallel")
#register(SnowParam(workers=detectCores()-1, progress=TRUE))#Windows
register(MulticoreParam(workers=detectCores()-1, progress=TRUE))#Linux
options(width=80)
###############################################################################
{## 1) Read Normal Data
##      -Exploring the first sample
##      -Check if all samples have the same size
##      -Check if the genes match positions
##      -Let's keep only the raw counts
##      -Let's change the annotation 
##      -Save clean data
##############################################################################
cat('Checking normal samples \n')
normalFiles<-dir(full.names=TRUE, pattern = "*.htseq.counts",
    path= paste(DATADIR, "counts/healthy", sep=""))
normal<-bplapply(normalFiles, read.delim, sep="\t", header=F, col.names=c("EnsemblID", "raw_counts"))

##Check if all samples have the same size
size<-unique(do.call(rbind,lapply(normal, dim)))
stopifnot(nrow(size)==1)
cat('Normal samples have the same size \n')

##Check if the genes match positions
genes<-do.call(cbind,lapply(normal, function(x)as.character(x[,1])))
genes<-t(unique(t(genes)))
stopifnot(dim(genes)==c(size[1,1], 1))
cat('Genes in normal samples match positions \n')

##Let's keep only the raw counts
normal<-do.call(cbind, lapply(normal, function(x)x[,"raw_counts"]))
targets<-data.frame(File=normalFiles, ID=paste("N", 1:length(normalFiles), sep=""))
colnames(normal)<-targets$ID

##Let's change the annotation 
genes<-do.call(rbind,sapply(genes[,1], strsplit, split=".", fixed=TRUE))
colnames(genes)<-c("EnsemblID", "version")

##Save clean data
normal<-list(Counts=normal, Annot=genes, targets=targets)
save(normal, file=paste(RDATA, "NormalRaw.RData", sep="/"), compress="xz")
cat('NormalRaw.RData saved \n')
}##############################################################################
{## 2) Read Tumour Data
##      -Exploring the first sample
##      -Check if all samples have the same size
##      -Check if the genes match positions
##      -Let's keep only the raw counts
##      -Let's change the annotation
##      -Save clean data
##############################################################################
cat('Checking tumor samples \n')
tumorFiles<-dir(full.names=TRUE, pattern = "*.htseq.counts",
    path=paste(DATADIR, "counts/cancer", sep=""))
tumor<-bplapply(tumorFiles, read.delim, sep="\t", header=F, col.names=c("EnsemblID", "raw_counts"))

##Check if all samples have the same size
size<-unique(do.call(rbind,lapply(tumor, dim)))
stopifnot(nrow(size)==1)
cat('Tumor samples have the same size \n')

##Check if the genes match positions
genes<-do.call(cbind,lapply(tumor, function(x)as.character(x[,1])))
genes<-t(unique(t(genes)))
stopifnot(dim(genes)==c(size[1,1], 1))
cat('Genes in tumor samples match positions \n')

##Lets keep only the raw counts
tumor<-do.call(cbind, lapply(tumor, function(x)x[,"raw_counts"]))
targets<-data.frame(File=tumorFiles, ID=paste("T", 1:length(tumorFiles), sep=""))
colnames(tumor)<-targets$ID

##Let's change the annotation
genes<-do.call(rbind,sapply(genes[,1], strsplit, split=".", fixed=TRUE))
colnames(genes)<-c("EnsemblID", "version")

tumor<-list(Counts=tumor, Annot=genes, targets=targets)
save(tumor, file=paste(RDATA, "TumorRaw.RData", sep="/"), compress="xz")
cat('TumorRaw.RData saved \n')
}##############################################################################
{## 3) Read Biomart Data
##      -Read the file
##      -Exploring the header
##      -Remove not conventional chromosomes
##      -Keep only annotated EntrezID individuals
##      -Remove the ones without Symbol | HGNCID
##      -Save clean data
###############################################################################
## Read the file
cat('Working with annotation file: Biomart_EnsemblG93_GRCh38_p12_NoSymbol.txt \n')
annot<-read.delim(file="/pipeline/Biomart_EnsemblG93_GRCh38_p12_NoSymbol.txt", sep="\t")

names(annot)<-c("EnsemblID", "Chr", "Start", "End", "GC", "Type")
annot$Length<-abs(annot$End - annot$Start)

## Exploring the header
head(annot)
#         EnsemblID                         Chr     Start       End    GC     Symbol           Type Length
# 1 ENSG00000283891                          15  55372940  55373034 41.05     MIR628          miRNA     94
# 2 ENSG00000251931                          12  59450673  59450772 33.00  RNU6-871P          snRNA     99
# 3 ENSG00000207766                          15  41691585  41691678 31.91     MIR626          miRNA     93
# 4 ENSG00000275323 CHR_HSCHR19LRC_LRC_J_CTG3_1  54201473  54208260 50.16 AC012314.7 protein_coding   6787
# 5 ENSG00000276678                           3  10287413  10287482 47.14     GHRLOS       misc_RNA     69
# 6 ENSG00000207260                           4 109992325 109992431 43.93   RNU6-35P          snRNA    106

##Remove not conventional chromosomes
annot<-as.data.frame(annot)
levels(annot$Chr)
annot<-annot[annot$Chr%in%c(as.character(1:22), "X", "Y"),]
annot$Chr<-droplevels(annot$Chr)
cat('Non conventional chromosomes removed \n')

uniq.annot <- length(unique(annot$EnsemblID)) == nrow(annot)
if(uniq.annot) {
  cat('Unique EnsemblIDs in annotation file\n')
} else {
  cat('Repeated EnsemblIDs in annotation file. \n')
  stop()
}

cat(paste('Annotation file. Final dimension: ', paste(dim(annot), collapse=", "), '\n'))
## Save clean data
save(annot, file=paste(RDATA, "annot.RData", sep="/"), compress="xz")
cat('annot.RData saved \n')
}##############################################################################
{##4) Merging count and annotation
##	    -M=normal|tumour
##      -targets=normal+tumor
##          -Check M y targets integrity
##	    -Annot=genes|annot
##          -Check if genes from normal and tumour match
##          -Are there repeated EntrezID?
##          -Add Biomart data
##          -Are there duplicated IDs?
##              -Keep the matching Symbol.x==Symbol.y
##              -Keep the matching Symbol.x==Symbol.y
##              -Keep the one with lowest GC content
##          -Did we lost any row in the process?
##          -Save the clean Data
##############################################################################
cat('Merging counts and annotations \n')
load(file=paste(RDATA, "annot.RData", sep="/"))
load(file=paste(RDATA, "TumorRaw.RData", sep="/"))
load(file=paste(RDATA, "NormalRaw.RData", sep="/"))

##M=normal|tumor
M<-cbind(normal$Counts, tumor$Counts)
cat(paste('Total number of features and samples: ', paste(dim(M), collapse=" ,"), '\n'))
# [1] 60488     6

##targets=normal+tumor
targets<-rbind(normal$targets, tumor$targets)
targets<-data.frame(targets, stringsAsFactors=FALSE)
targets$File<-as.character(targets$File)
targets$ID<-as.character(targets$ID)
dim(targets)
# [1] 6   2

##check M y targets integrity
stopifnot(nrow(targets)==ncol(M))
cat('Number of counts columns match sample number\n')

## Annot=genes|annot
## check if genes from normal and tumor match
stopifnot(all(row.names(normal$Annot)==row.names(tumor$Annot)))
cat('Genes from normal and tumor samples match \n')

Annot<-normal$Annot
##Are there repeated EnsemblID?
stopifnot(length(unique(Annot[, "EnsemblID"]))==nrow(Annot))
cat('No duplicated EnsemblIDs\n')

Annot<-data.frame(Annot, stringsAsFactors=FALSE)
Annot$Row<-1:nrow(Annot) ##Just to maintain the original order
head(Annot)
#                          EnsemblID version Row
# ENSG00000000003.13 ENSG00000000003      13   1
# ENSG00000000005.5  ENSG00000000005       5   2
# ENSG00000000419.11 ENSG00000000419      11   3
# ENSG00000000457.12 ENSG00000000457      12   4
# ENSG00000000460.15 ENSG00000000460      15   5
# ENSG00000000938.11 ENSG00000000938      11   6

##Add Biomart data
cat('Adding biomart data\n')
Annot<-merge(x=Annot, y=annot, by.x="EnsemblID", by.y="EnsemblID",
  all=FALSE, all.x=TRUE, all.y=FALSE, sort=FALSE)
cat(paste('Merged file. Final dimensions: ', paste(dim(Annot), collapse=", "), '.\n'))
# [1] 60488    10
dim(M)
# [1] 60488     6
extra.rows <- nrow(Annot)-nrow(M) 
cat(paste('There are ', extra.rows, ' extra rows in the counts matrix.\n'))
# [1] 0

##Are there duplicated IDs?
Annot<-Annot[order(Annot$Row, Annot$Chr),]
ids<-head(which(duplicated(Annot$Row)))
ids
if (identical(ids, integer(0))) {
  cat('There are no duplicated IDs \n')
} else {
  Annot[sort(c(ids-1,ids)),][1:3,]
  #     EntrezID Symbol.x Row Chr     Start       End    GC           Type
  # 373     2334     AFF2 399   X 148500619 149000663 38.84 protein_coding
  # 374     2334     AFF2 399   X 148501098 148991332 38.79 protein_coding
  # 387   119016    AGAP4 413  10  45825594  45853875 41.11 protein_coding
  #         HGNCID Symbol.y Length
  # 373  HGNC:3776     AFF2 500044
  # 374  HGNC:3776     AFF2 490234
  # 387 HGNC:23459    AGAP4  28281
  
  ##Keep the matching Symbol.x==Symbol.y
  ids<-which(duplicated(Annot$Row))
  ids
  # integer(0)
  
  noneMatchingSymbols<-ids[Annot$Symbol.x[ids]!=Annot$Symbol.y[ids]]
  head(Annot[noneMatchingSymbols,],3)
  #      EntrezID Symbol.x  Row Chr    Start      End    GC           Type
  # 388    119016    AGAP4  413  10 49982190 50010499 41.10 protein_coding
  # 960     51326   ARL17A  998  17 46274784 46361797 43.65 protein_coding
  # 1667    79741 C10orf68 1720  10 32567723 32882874 37.18 protein_coding
  #          HGNCID Symbol.y Length
  # 388  HGNC:23466    AGAP6  28309
  # 960  HGNC:32387   ARL17B  87013
  # 1667 HGNC:26533    CCDC7 315151
  # 
  Annot<-Annot[-noneMatchingSymbols, ]
  dim(Annot)
  # [1] 20590    11
  nrow(Annot)-nrow(M)
  # [1] 58 Still duplicated
  
  ##Keep the matching Symbol.x==Symbol.y
  Annot<-Annot[order(Annot$Row, Annot$Chr),]
  ids<-which(duplicated(Annot$Row))
  ids<-sort(unique(c(ids-1, ids)))
  noneMatchingSymbols<-ids[Annot$Symbol.x[ids]!=Annot$Symbol.y[ids]]
  head(Annot[noneMatchingSymbols,], 3)
  #      EntrezID  Symbol.x  Row Chr    Start      End    GC           Type
  # 666    441425 ANKRD20A3  695   9 67859147 67902094 36.24 protein_coding
  # 1238    84126     ATRIP 1279   3 48446710 48467645 48.87 protein_coding
  # 2927   388389   CCDC103 3065  17 44899712 44905390 51.70 protein_coding
  #          HGNCID  Symbol.y Length
  # 666  HGNC:23665 ANKRD20A1  42947
  # 1238 HGNC:12269     TREX1  20935
  # 2927 HGNC:35153   FAM187A   5678
  
  Annot<-Annot[-noneMatchingSymbols, ]
  dim(Annot)
  # [1] 20543    11
  nrow(Annot)-nrow(M)
  # [1] 11 Still duplicated
  
  ##Keep the one with lowest GC content
  Annot<-Annot[order(Annot$Row, Annot$Chr, Annot$GC),]
  ids<-which(duplicated(Annot$Row))
  Annot<-Annot[-ids,]
  
  dim(Annot)
  # [1] 20532    11
  nrow(Annot)-nrow(M)
  # [1] 0
  ##Finally!!!
  
}

gctable <- table(is.na(Annot$GC))
cat("There are", gctable[[1]], "entries with GC info and", gctable[[2]], "GC entries missing \n")
# FALSE  TRUE 
# 56963  3525

##Did we lost any row in the process?
stopifnot(all(Annot$Row %in% 1:nrow(M)))
stopifnot(all(1:nrow(M) %in% Annot$Row))

##Save the clean Data
full<-list(M=M, Annot=Annot, Targets=targets)
#save(full, file=paste(RDATA, "RawFull.RData", sep="/"), compress="xz")
cat("Saving RawFull.RData \n")
}##############################################################################
## GREAT JOB!!! YOU MADE IT TILL THE END!!!!
###############################################################################
