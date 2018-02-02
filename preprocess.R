WORKDIR <- '/home/diana/Workspace/rnapipeline/breast_cancer_pipeline/data/counts/cancer'
#args <- commandArgs(trailingOnly = TRUE)
#WORKDIR <- args[1]
cat('Setting WORKDIR to:', WORKDIR, '\n' )
setwd(WORKDIR)

fns <- system2('ls', args = '*.counts', stdout = T)
fns
counts.df <- NULL
for (fn in fns) {
  df <- read.table(fn, check.names=F, sep='\t', col.names=c("Geneid", "Counts"))
  colnames(df)[2] <- strsplit(basename(fn), '\\.')[[1]][1] # beautify sample name
  if (is.null(counts.df)) {
    counts.df <- df
  } else {
    counts.df <- merge(counts.df, df, by="Geneid", all.x=T, sort=F)
  }
}

cat('Finished reading featureCount into data.frame with shape: ', dim(counts.df), '\n' )
