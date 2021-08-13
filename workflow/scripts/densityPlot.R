log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(readr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggthemes)
TISSUE <- snakemake@params[["tissue"]]
w <- 1024
h <- 1024
p <- 24

cat("Loading data \n")
load(snakemake@input[[1]])

cat("Calculating mean \n")
melted_data <- melt(log2(full$M+1))
means <- colMeans(full$M)
mm <- mean(means)
sd <- sd(means)

down_outliers <- which(means < mm-2*sd)
up_outliers <- which(means > mm+2*sd)

if(length(down_outliers) > 0 | length(up_outliers) > 0 ) {
  cat("Outliers found, saving samples names \n")
} else {
  cat("Outliers not found\n")
}

outliers <- data.frame(sample = c(names(down_outliers), names(up_outliers)),
           mean = c(means[down_outliers], means[up_outliers]))
outliers <- outliers %>% mutate(sd_from_mean = (mean - mm)/sd )

write_tsv(outliers, file=snakemake@output[["outliers"]])

pl <- ggplot(data = melted_data, aes(x=value, group=Var2, colour=Var2)) + 
  geom_density(show.legend = F) + 
  theme_hc(base_size = 20) +
  xlab("log2(expr+1)") + 
  ggtitle(TISSUE)

png(file = snakemake@output[["density"]], width = w, height = h/2, pointsize = p)
print(pl)
dev.off()

cat("Density plot for all samples generated\n")


