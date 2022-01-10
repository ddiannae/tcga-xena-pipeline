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
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(limma)
library(edgeR)
library(Glimma)
library(readr)
library(dplyr)
library(janitor)
library(ggplot2)
library(ggthemes)

DEGDIR <- paste(snakemake@params[["deg_dir"]])
dir.create(DEGDIR)

TISSUE <- snakemake@params[["tissue"]]

cat("Loading files\n")
load(snakemake@input[["full"]])
load(snakemake@input[["annot"]])

full$targets <- full$targets %>% 
  dplyr::mutate(group = factor(group, levels=c("normal", "cancer")))

cat("Generating model matrix\n")
mm <- model.matrix(~1+group, data = full$targets)
rownames(mm) <- full$targets %>% pull(id)
full$M <- full$M[, rownames(mm)]

cat("Getting voom transformation\n")
y <- voom(full$M, mm, plot = F, save.plot = T)

p <- ggplot() +
  geom_point(aes(x = y$voom.xy$x, y = y$voom.xy$y), size = 0.5) + 
  geom_line(aes(x = y$voom.line$x, y = y$voom.line$y), size = 1, color = "red") +
  xlab(y$voom.xy$xlab) +
  ylab(y$voom.xy$ylab) + 
  theme_bw() +
  ggtitle(TISSUE)

cat("Saving voom plot\n")
png(file = paste0(DEGDIR, "/voom.png"), width = 800, height = 600)
print(p)
dev.off()

cat("Getting linear adjustment\n")
fit <- lmFit(y, mm)
lfc_fit <- eBayes(fit)

cat("Getting Glimma interactive plot plot\n")
dt <- decideTests(lfc_fit, lfc = 2) 

annot <- annot %>% filter(gene_id %in% rownames(full$M)) %>% 
  select(gene_id, ensembl_id, chr, start, end, gene_name)

glMDPlot(path = DEGDIR, lfc_fit, 
         counts = full$M, groups = full$targets %>% pull(group), 
         status = dt, anno = annot)

p <- ggplot() +
  geom_density(aes(x = lfc_fit$coefficients[,2]), size = 0.5) + 
  theme_bw() +
  xlab("lfc")
  ggtitle(TISSUE)

png(file = paste0(DEGDIR, "/density.png"), width = 800, height = 600)
print(p)
dev.off()
  
cat("Saving deg results\n")
topTable(lfc_fit, sort.by = "P", n = Inf, adjust.method="BH") %>%
  janitor::clean_names() %>% mutate(gene_id = rownames(.)) %>%
  inner_join(annot %>% select(gene_id, ensembl_id, chr, gene_name), by ="gene_id") %>% select(-gene_id) %>%
  select(ensembl_id, everything())  %>%
  mutate(adj_p_val = ifelse(adj_p_val > 0.05, 1, adj_p_val), 
                           log_fc = ifelse(adj_p_val > 0.05, 0, log_fc)) %>%
  write_tsv(snakemake@output[["deg_results"]])
