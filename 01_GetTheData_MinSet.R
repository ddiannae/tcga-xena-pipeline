load(file=paste(RDATA, "annot.RData", sep="/"))
load(file=paste(RDATA, "TumorRaw.RData", sep="/"))
load(file=paste(RDATA, "NormalRaw.RData", sep="/"))

##M=normal|tumor
M<-cbind(normal$Counts, tumor$Counts)
dim(M)

##targets=normal+tumor
targets <- rbind(normal$targets, tumor$targets)
targets <- data.frame(targets, stringsAsFactors=FALSE)
targets$File <- as.character(targets$File)
targets$ID <- as.character(targets$ID)
targets$Case <- as.character(targets$Case)
dim(targets)
# [1] 6   2

##check M y targets integrity
stopifnot(nrow(targets)==ncol(M))

## Annot=genes|annot
## check if genes from normal and tumor match
stopifnot(all(row.names(normal$Annot)==row.names(tumor$Annot)))
cat('Genes from normal and tumor samples match \n')

Annot<-normal$Annot
##Are there repeated EnsemblID?
stopifnot(length(unique(Annot[, "EnsemblID"]))==nrow(Annot))
Annot<-data.frame(Annot, stringsAsFactors=FALSE)
Annot$Row<-1:nrow(Annot) ##Just to maintain the original order

Annot<-merge(x=Annot, y=annot, by.x="EnsemblID", by.y="EnsemblID",
  all=FALSE, all.x=TRUE, all.y=FALSE, sort=FALSE)
dim(Annot)
# [1] 60488    10
dim(M)
# [1] 60488     6
nrow(Annot)-nrow(M)
# [1] 0

##Are there duplicated IDs?
Annot<-Annot[order(Annot$Row, Annot$Chr),]
ids<-head(which(duplicated(Annot$Row)))
ids
if (identical(ids, integer(0))) {
  cat('There are no duplicated IDs \n')
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
save(full, file=paste(RDATA, "RawFull.RData", sep="/"), compress="xz")
cat("Saving RawFull.RData \n")