library(infotheo)
library(readr)
library(dplyr)
library(parallel)
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3 ) {
  stop("Incorrect number of arguments", call.=FALSE)
} else {
  MATRIX <- args[1]
  OUTFILE <- args[2]
  MCCORES <- args[3]
}

matrix <- read_tsv(MATRIX)
genes <- matrix$gene
matrix <- matrix %>% select(-gene) %>% discretize()

max <- length(genes)
#max <- 1000
lines <-  mclapply(X = 1:(max-1),  mc.cores = MCCORES, FUN = function(i) {
#lines <-  lapply(X = 1:(max-1), FUN = function(i) {
  gene1 <- unlist(matrix[i, ], use.names = F)
  mis <- lapply((i+1):max, function(j) {
    gene2 <-  unlist(matrix[j, ])
    return(round(mutinformation(gene1, gene2), 5))
  })
  line <- paste0(c(rep(" ", times = i), unlist(mis)), collapse = ",")
  return(line)
})

outFile <- file(OUTFILE)
writeLines(c(paste0(genes,  collapse = ","), unlist(lines)), outFile)
close(outFile)

