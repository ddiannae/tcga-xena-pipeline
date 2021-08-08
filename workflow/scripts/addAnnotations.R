log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")
## Get the Work and Data dir
library(readr)
library(dplyr)
library(rtracklayer)

RDATADIR <- paste(snakemake@params[["tissue_dir"]], "rdata", sep="/")
dir.create(RDATADIR)

IS_XENA <- snakemake@params[["is_xena"]]
##############################################################################
## Get annotation data
## 1. Read the original annotation data
## 2. Read the new annotation data
## 3. Keep only genes/protein-coding present in original and new data
## 4. Query bioMart ensembl 80 for GC content 
## 5. Save data
###############################################################################
if(!IS_XENA) {
  cat("Getting annotation file \n")
  cat("Reading original file \n")
  ## Annotation file used for RNA-seq pipeline according to:
  ## https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files
  annot <- rtracklayer::import('input/gencode.v22.annotation.gtf.gz')
  annot <- as.data.frame(annot)
  
  ## Only protein coding genes
  annot <- annot %>% dplyr::select(gene_id, seqnames, start, end, width, type, 
                                   gene_type, gene_name) %>% 
    filter(type == "gene" & gene_type == "protein_coding")
  
  annot <- annot %>% mutate(ensembl_id = unlist(lapply(strsplit(gene_id, "\\."), "[[", 1))) %>%
    select(-gene_type)

} else {
  annot <- read_tsv(snakemake@input[["xena_annot"]],
                    col_names  = c("gene_id", "gene_name", "seqnames", "start", "end", "strand"), 
                    skip = 1)
  annot <- annot %>% dplyr::mutate(ensembl_id = unlist(lapply(strsplit(gene_id, "\\."), "[[", 1)), 
                                   width = end - start)
  join_by <- "ensembl_id"
}

cat("Reading new file \n")
## Newest gencode file. April, 2021.
annot_new <-  rtracklayer::import('input/gencode.v37.annotation.gtf.gz')
annot_new <- as.data.frame(annot_new)
annot_new <- annot_new %>% dplyr::select(gene_id, gene_name, type, gene_type) %>% 
  dplyr::filter(type == "gene" & gene_type == "protein_coding")

annot_new <- annot_new %>% dplyr::mutate(ensembl_id = unlist(lapply(strsplit(gene_id, "\\."), "[[", 1)))

## Keep genes that remain in the newest annotation file
## but get the newest names and keep only conventional chromosomes
## remove duplicates
annot <- annot %>% dplyr::select(-gene_name) %>%
  inner_join(annot_new %>% dplyr::select(ensembl_id, gene_name, gene_type), by = "ensembl_id") %>%
  filter(seqnames %in% paste0("chr", c(as.character(1:22), "X", "Y"))) %>% distinct()

cat('Annotation file new/old merge: ', paste(dim(annot), collapse=", "), '\n')

## We need GC content per gene for normalization. 
## Query Biomart 80 (accoring to gencode.v22.annotation.gtf, version 79 was used
## but it is no longer accessible in the website)
## http://may2015.archive.ensembl.org/biomart
biomart <- read_tsv("input/Biomart_Ensembl80_GRCh38_p2.txt", 
                    col_names = c("ensembl_id", "version", "gc"), skip = 1)

## Get only genes matching ensemblID
biomart <- biomart %>% mutate(gene_id = paste(ensembl_id, version, sep="."))

if(!IS_XENA) {
  annot <- annot %>%
    inner_join(biomart %>% dplyr::select(gene_id, gc), by = "gene_id")
} else {
  annot <- annot %>%
    inner_join(biomart %>% dplyr::select(ensembl_id, gc), by = "ensembl_id")
}

annot <- annot %>% mutate(chr = gsub("chr", "", seqnames)) %>%
  dplyr::rename(length = width) %>%
  dplyr::select(-seqnames) %>%   dplyr::select(gene_id, chr, everything())

cat('Annotation file. Final dimension: ', paste(dim(annot), collapse=", "), '\n')
save(annot, file=snakemake@output[["annot_rdata"]], compress="xz")
write_tsv(annot, snakemake@output[["annot_tsv"]])
cat('annot.RData saved \n')

##############################################################################
## Merging count and annotation
## 1. Build counts matrix. M=normal|tumour
## 2. Build targets matrix. targets=normal+tumor
## 3. Check M y targets integrity
## 4. Filter by annotation file 
## 5. Save the clean data
##############################################################################
{
  cat('Merging counts and annotations \n')
 
  normal_samples <- list(matrix=read_tsv(snakemake@input[["normal_matrix"]]),
                         targets=read_tsv(snakemake@input[["normal_targets"]])) 
  cancer_samples <- list(matrix=read_tsv(snakemake@input[["cancer_matrix"]]),
                         targets=read_tsv(snakemake@input[["cancer_targets"]])) 
  ## Raw counts
  M <- normal_samples$matrix %>% inner_join(cancer_samples$matrix, by = "gene_id")
  cat('Total number of features and samples: ', paste(dim(M), collapse=" ,"), '\n')
  
  ## Samples
  targets <- bind_rows(normal_samples$targets, cancer_samples$targets)
  
  ## Check M y targets integrity. Remove gene_ids col
  stopifnot(nrow(targets) == (ncol(M)-1))
  cat('Number of counts columns match sample number\n')
  
  ## Filter counts by annotation data
  cat('Adding biomart data\n')
  M <- M %>% semi_join(annot, by = "gene_id")
  annot <- annot %>% semi_join(M, by = "gene_id")
  
  cat('Total number of annotated (genes/protein-coding) features:', nrow(M), '\n')
  cat('Total number of samples:', ncol(M)-1, '\n')
  
  no_gc <- sum(is.na(annot$gc))
  cat("There are",no_gc, "entries with no GC info \n")
  
  ids <- annot %>% filter(!is.na(gc) & !is.na(length)) %>%
    select(gene_id) %>% unlist(use.names = F)
  
  M <- M %>% filter(gene_id %in% ids) %>% arrange(gene_id)
  annot <- annot %>% filter(gene_id %in% ids) %>% arrange(gene_id)
  cat("Non GC and lenght annotated genes removed.\n")
  
  ## Save it as a matrix
  ids <- M$gene_id
  MM <- M %>% select(-gene_id) %>% as.matrix() 
  rownames(MM) <- ids
  MM <- MM[,targets$id]
  
  ## Make sure they are factors
  targets$group <- factor(targets$group, levels=c("cancer", "normal"))

  ##Save clean data
  cat('Saving raw full data \n')
  full <- list(M = MM, annot = annot, targets = targets)
  
  save(full, file=snakemake@output[["raw_rdata"]], compress="xz")
  cat("raw_full.RData saved \n")
}
