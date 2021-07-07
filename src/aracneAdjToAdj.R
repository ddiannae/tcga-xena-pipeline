library(data.table)
library(parallel)
library(plyr)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3 ) {
  stop("Incorrect number of arguments", call.=FALSE)
} else {
  ADJ_IN = args[1]
  ADJ_OUT = args[2]
  MCCORES = as.numeric(args[3])
}

cat("Loading file\n")
aracne_inter <- fread(ADJ_IN, nThread = MCCORES)
genes_row <- unlist(aracne_inter$V1)

cindex <- c(1, seq(2, ncol(aracne_inter), by=2))
genes <-  as.character(unlist(aracne_inter[1, ..cindex]))
fgene <- genes[!genes %in% genes_row]
cindex <- seq(3, ncol(aracne_inter), by=2)
aracne_inter_num <- as.matrix(aracne_inter[, ..cindex])

rownames(aracne_inter_num) <- genes_row
aracne_inter_num <- aracne_inter_num[sort(genes_row), ]
aracne_inter_num  <- rbind(t(aracne_inter_num[, 1]),aracne_inter_num)
rownames(aracne_inter_num)[1] <- fgene
aracne_inter_num  <- cbind(aracne_inter_num, as.numeric(aracne_inter_num[, ncol(aracne_inter_num)]))
colnames(aracne_inter_num) <- rownames(aracne_inter_num)
aracne_inter_num[upper.tri(aracne_inter_num, diag = TRUE)] <- NA
aracne_inter_num <- t(aracne_inter_num)
aracne_inter_num <- aracne_inter_num[1:(nrow(aracne_inter_num)-1),]
cat("Saving file\n")
fwrite(aracne_inter_num, file = ADJ_OUT, row.names = F, col.names = T, sep = ",", nThread = MCCORES)

# n <- nrow(aracne_inter_num)
# r <- lapply(X = 2:(n-1), 
#                        FUN = function(i){
#                          v = c(as.numeric(aracne_inter_num[i, 1:(i-1)]), 
#                                NA, 
#                                as.numeric(aracne_inter_num[i, i:ncol(aracne_inter_num)]))
#                        })
# cat("Parallel computations done\n")
# r <-ldply(r)
# r <- rbind(c(NA, aracne_inter_num[1,]), r,  c(aracne_inter_num[nrow(aracne_inter_num),], NA))
# r <- as.matrix(r)
# rt <- t(r) 
# diag(rt) <- 1
# diag(r) <- 1
# all(r == rt)


