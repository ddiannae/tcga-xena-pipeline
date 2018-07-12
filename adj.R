setwd("./Sanos10")
adj<-dir()[regexpr(pattern=".adj", dir())>0]


makeSIF <- function(x){
  # args - 
  #    x - m*m distance or correlation matrix
  # @returns data frame in SIF format 
  #
  sif <- as.data.frame(t(combn(as.character(rownames(x)), 2)))
  #print(sif)
  weight <- apply(sif, 1, indexDMatFromLookup, x)
  sif2 <- data.frame(sif, weight)
  return(sif2) 
}

indexDMatFromLookup <- function(lookup, x) {
  return(indexDMat(x, lookup[1], lookup[2]))
}

indexDMat <- function(x, i1,i2) {
  return(x[i1,i2])
}

