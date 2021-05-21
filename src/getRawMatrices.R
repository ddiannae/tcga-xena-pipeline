###############################################################################
##  LOSS OF LONG RANGE CO-EXPRESSION IN CANCER
##  Analysis of RNA-Seq data.
###############################################################################
## Diana Garcia - diana.gco@gmail.com
## Date: April, 2021
## 
## Original code. Dr. Cristobal Fresno - cristobalfresno@gmail.com
## Date:  2016-12-12
###############################################################################
## 1) Read Normal Data
## 2) Read Cancer Data
## 3) Read Annotation Data
## 4) Merge count and annotation
###############################################################################

## Get the Work and Data dir
library(readr)
library(dplyr)
library(rtracklayer)

args <- commandArgs(trailingOnly = T)

if (length(args) < 2 ) {
  stop("Incorrect number of arguments", call.=FALSE)
} else {
  TISSUE = args[1]
  DATADIR = args[2]
}
DATADIR <- paste(DATADIR, TISSUE, sep="/")
RAWDIR <- paste(DATADIR, "raw", sep="/")
RDATADIR <- paste(DATADIR, "rdata", sep="/")
dir.create(RDATADIR)
###############################################################################
## Get and check count matrices
## 1. Find files
## 2. Read data
## 4. Check sample sizes
## 5. Check genes order
## 6. Build target matrix
###############################################################################
{
  getMatrix <- function(type) {
    cat(paste0("Checking ", type, " samples \n"))
    
    ## Find the files
    files_to_read <- list.files(path = paste0(RAWDIR,"/", TISSUE, "-", type, "-rna"), 
                                pattern = "\\.htseq.counts.gz$", full.names = T, recursive = T)
    
    ## Read the data
    all_files <- lapply(files_to_read, function(file) {
      data <- read_tsv(file, col_names = c("gene_id", "raw_counts"))
      data <- data %>% filter(!startsWith(gene_id, "_" ))
      return(data)
    })
    
    ## Check samples sizes
    size <- unique(do.call(rbind,lapply(all_files, dim)))
    stopifnot(nrow(size)==1)
    cat(paste0(type, " samples have the same size \n"))
    
    ## Check genes order in samples
    genes <- do.call(cbind, lapply(all_files, function(x) dplyr::select(x, "gene_id")))
    genes <- t(unique(t(genes)))
    stopifnot(dim(genes)==c(size[1,1], 1))
    cat(paste0("Genes in ", type, " samples match positions \n"))
    
    ## Build targets matrix
    targets <- data.frame(id = paste(TISSUE, type, 1:length(files_to_read), sep = "_"), 
                          file = unlist(lapply(strsplit(files_to_read, "/"), "[[", 10)),
                          file_id = unlist(lapply(strsplit(files_to_read, "/"), "[[", 11)),
                          group = type, stringsAsFactors = FALSE)
    
    ## Rename columns in counts matrix
    matrix <- bind_cols(lapply(all_files, function(x) dplyr::select(x, "raw_counts")))
    colnames(matrix) <- targets$id
    
    matrix <- matrix %>% mutate(gene_id = genes[,1]) %>% 
      dplyr::select(gene_id, everything())
    
    return(list(targets = targets, matrix = matrix))
    
  }
  
  ## Check normal samples
  normal_samples <- getMatrix("normal")
  ## Check normal samples
  cancer_samples <- getMatrix("cancer")
  
  cat("Saving normal samples \n")
  save(normal_samples, file=paste(RDATADIR, "raw_normal.RData", sep="/"), compress="xz")
  write_tsv(normal_samples$matrix, paste0(DATADIR, "/", TISSUE, "-normal-matrix.tsv"))
  write_tsv(normal_samples$targets, paste0(DATADIR,"/", TISSUE, "-normal-samples.tsv"))
  
  cat("Saving cancer samples \n")
  save(cancer_samples, file=paste(RDATADIR, "raw_cancer.RData", sep="/"), compress="xz")
  write_tsv(cancer_samples$matrix, paste0(DATADIR, "/", TISSUE, "-cancer-matrix.tsv"))
  write_tsv(cancer_samples$targets, paste0(DATADIR,"/", TISSUE, "-cancer-samples.tsv"))
  cat('raw_normal.RData, raw_cancer.RData and matrices saved \n')
}
##############################################################################
## Get annotation data
## 1. Read the original annotation data
## 2. Read the new annotation data
## 3. Keep only genes/protein-coding present in original and new data
## 4. Query bioMart ensembl 80 for GC content 
## 5. Save data
###############################################################################
{
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
  
  annot <- annot %>% mutate(ensembl_id = unlist(lapply(strsplit(gene_id, "\\."), "[[", 1)))
  
  cat("Reading new file \n")
  ## Newest gencode file. April, 2021.
  annot_new <-  rtracklayer::import('input/gencode.v37.annotation.gtf.gz')
  annot_new <- as.data.frame(annot_new)
  annot_new <- annot_new %>% dplyr::select(gene_id, gene_name, type, gene_type) %>% 
    filter(type == "gene" & gene_type == "protein_coding")

  annot_new <- annot_new %>% mutate(ensembl_id = unlist(lapply(strsplit(gene_id, "\\."), "[[", 1)))
  
  ## Keep genes that remain in the newest annotation file
  ## but get the newest names and keep only conventional chromosomes
  ## remove duplicates
  annot <- annot %>% dplyr::select(-gene_name) %>%
    inner_join(annot_new %>% dplyr::select(ensembl_id, gene_name), by = "ensembl_id") %>%
    filter(seqnames %in% paste0("chr", c(as.character(1:22), "X", "Y"))) %>% distinct()
  
  cat('Annotation file new/old merge: ', paste(dim(annot), collapse=", "), '\n')
  
  ## We need GC content per gene for normalization. 
  ## Query Biomart 80 (accoring to gencode.v22.annotation.gtf, version 79 was used
  ## but it is no longer accessible in the website)
  ## http://may2015.archive.ensembl.org/biomart
  biomart <- read_tsv("input/Biomart_Ensembl80_GRCh38_p2.txt", 
                      col_names = c("ensembl_id", "version", "gc"), skip = 1)
  
  ## Get only genes matching ensemblID and version
  biomart <- biomart %>% mutate(gene_id = paste(ensembl_id, version, sep="."))
  annot <- annot %>%
    inner_join(biomart %>% dplyr::select(gene_id, gc), by = "gene_id")
   
  annot <- annot %>% mutate(chr = gsub("chr", "", seqnames)) %>%
    dplyr::rename(length = width) %>%
    dplyr::select(-seqnames) %>%   dplyr::select(gene_id, chr, everything())
  
  cat('Annotation file. Final dimension: ', paste(dim(annot), collapse=", "), '\n')
  
  save(annot, file=paste(RDATADIR, "annot.RData", sep="/"), compress="xz")
  write_tsv(annot, paste0(DATADIR,"/", TISSUE, "-annotation.tsv"))
  cat('annot.RData saved \n')
}
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
  
  no_gc <- sum(is.na(annot$percentage_gc_content))
  cat("There are",no_gc, "entries with no GC info \n")
  
  ids <- annot %>% filter(!is.na(gc) & !is.na(length)) %>%
    select(gene_id) %>% unlist(use.names = F)
  
  M <- M %>% filter(gene_id %in% ids) %>% arrange(gene_id)
  annot <- annot %>% filter(gene_id %in% ids) %>% arrange(gene_id)
  cat("Non GC and lenght annotated genes removed.\n")
  
  rownames(annot) <- annot$gene_id
  
  ## Save it as a matrix
  ids <- M$gene_id
  MM <- M %>% select(-gene_id) %>% as.matrix() 
  rownames(MM) <- ids
  MM <- MM[,targets$id]
  
  ## Make sure they are factors
  targets$group <- factor(targets$group, levels=c("cancer", "normal"))
  rownames(targets) <- targets$id
  ##Save clean data
  cat('Saving raw full data \n')
  full <- list(M = MM, annot = annot, targets = targets)
  
  save(full, file=paste(RDATADIR, "raw_full.RData", sep="/"), compress="xz")
  write_tsv(M, paste0(DATADIR, "/", TISSUE, "-matrix.tsv"))
  write_tsv(targets, paste0(DATADIR, "/", TISSUE, "-samples.tsv"))
  
  cat("raw_full.RData and matrices saved \n")
}
