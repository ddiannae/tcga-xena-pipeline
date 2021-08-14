log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(readr)
library(dplyr)

cat("Loading files \n")
load(snakemake@input[[1]])
outliers <- read_tsv(snakemake@input[["outliers"]])

cat(nrow(outliers), " outliers will be removed.\n")
targets <- full$targets %>% 
  filter(!id %in% (outliers %>% pull(sample)))

M <- full$M[, targets %>% pull(id)]

full <- list(annot = full$annot, M = M, targets = targets)

dimcancer <- nrow(full$targets %>% filter(group == "cancer"))
dimnormal <- nrow(full$targets %>% filter(group == "normal"))
cat("Final dimensions: ", dimcancer, " cancer and ", dimnormal,
" normal samples \n")

cat("Saving RData \n")
save(full, file=snakemake@output[[1]])
