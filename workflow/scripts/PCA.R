log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(NOISeq)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(ggpubr)

w <- 1024
h <- 1024
p <- 24

PLOTSDIR <-paste(snakemake@params[["tissue_dir"]], "plots", 
                 snakemake@params[["plots_type"]], sep="/")
dir.create(PLOTSDIR, recursive = TRUE)

cat("Loading data")
load(snakemake@input[[1]])

cat("Performing PCA")
mydata <- NOISeq::readData(
  data = full$M, 
  length = full$annot %>% select(gene_id, length) %>% as.data.frame(), 
  biotype = full$annot %>% select(gene_id, gene_type) %>% as.data.frame(), 
  chromosome = full$annot %>% select(chr, start, end) %>% as.data.frame(), 
  factors = full$targets %>% select(group) %>% as.data.frame(),
  gc = full$annot %>% select(gene_id, gc) %>% as.data.frame())

pca_dat <- dat(mydata, type = "PCA", logtransf = F)
pca_results <- pca_dat@dat$result 
pca_factors <- pca_dat@dat$factors


cat("Getting variance plot")
pc_eigenvalues <- tibble(pc = factor(1:nrow(pca_results$var.exp)), 
                         pct = pca_results$var.exp[,1], 
                         pct_cum = pca_results$var.exp[,2])

g1 <- ggplot(data = pc_eigenvalues[1:50,], aes(x = pc)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained") +
  theme_hc(base_size = 20) +
  theme(axis.text.x = element_text(size=12)) +
  ggtitle("PCA Variance")

ggexport(g1, width = w, height = h/2, pointsize = p, 
         filename = snakemake@output[["variance"]])

cat("Getting scores plot")
pc_scores <- tibble(pc = factor(1:nrow(pca_results$var.exp)),
                    name = rownames(pca_factors),
                    group = factor(pca_factors$group, levels = c("normal", "cancer"), labels = c("Normal", "Cancer")), 
                    pc1 = pca_results$scores[,1],
                    pc2 = pca_results$scores[,2],
                    pc3 = pca_results$scores[,3])

color_pal <- c("#e3a098", "#a32e27")
xrange <- range(pc_scores %>% select(pc1))
yrange <- range(pc_scores %>% select(pc2, pc3))
  
g1 <- ggplot(pc_scores, aes(x=pc1, y=pc2, color=group)) +
  geom_point(size=2) +
  labs(x = "PC1", y = "PC2") +
  scale_color_manual(values = color_pal) +
  theme_hc(base_size = 20) +
  theme(legend.title = element_blank()) + 
  xlim(xrange) +
  ylim(yrange)

g2 <- ggplot(pc_scores, aes(x=pc1, y=pc3, color=group)) +
  geom_point(size=2) +
  labs(x = "PC1", y = "PC3") +
  scale_color_manual(values = color_pal) +
  theme_hc(base_size = 20) +
  theme(legend.title = element_blank()) +
  xlim(xrange) +
  ylim(yrange)

fig <- ggarrange(g1,  g2, nrow = 1)
fig <- annotate_figure(fig,
                top = text_grob("PCA Scores", size = 25))

ggexport(fig, width = w, height = h/2, pointsize = p, 
         filename = snakemake@output[["score"]])

cat("Getting loading plot")
pc_loadings <- tibble(gene = rownames(pca_results$loadings), 
                      pc1 = pca_results$loadings[, 1],
                      pc2 = pca_results$loadings[, 2])

top_genes <- pc_loadings %>% 
  pivot_longer(matches("pc"), names_to = "pc", values_to = "loading")%>% 
  group_by(pc) %>% 
  arrange(desc(abs(loading))) %>% 
  slice(1:10) %>% 
  pull(gene) %>% 
  unique()

top_loadings <- pc_loadings %>% 
  filter(gene %in% top_genes)

g1 <- ggplot(top_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = pc1, yend = pc2), 
               arrow = arrow(length = unit(0.1, "in")),
               color = "brown") +
  geom_text(aes(x = pc1, y = pc2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  ggtitle("PCA Loadings") +
  theme_hc(base_size = 15) 

ggexport(g1, width = w/2, height = h/2, pointsize = p, 
         filename = snakemake@output[["loading"]])
