library(data.table)
library(readr)
library(janitor)
library(dplyr)
library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 8) {
  stop("Incorrect number of arguments", call.=FALSE)
} else {
  XENA_COUNTS <- args[1]
  XENA_SAMPLES <- args[2]
  XENA_ANNOT <- args[3]
  TISSUE <- args[4]
  NORMAL_TISSUE <-args[5]
  PRIMARY_DISEASE <- args[6]
  DATADIR <- args[7]
  MCCORES <- as.numeric(args[8])
}
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
  matrix_samples <- read_tsv(XENA_SAMPLES) %>% clean_names()
  matrix_counts <-  fread(XENA_COUNTS, nThread = MCCORES)
  
  getSamples <- function(tissue, primary_disease, type) {
    extended_type <- ifelse(type == "normal", "Normal Tissue", "Primary Tumor")
    camelTissue <- paste0(toupper(substring(TISSUE, 1, 1)), substring(TISSUE, 2))
    
    targets <- matrix_samples %>% dplyr::filter(primary_site == camelTissue &
                                           primary_disease_or_tissue == primary_disease &
                                           sample_type == extended_type)
    matrix <- matrix_counts %>% dplyr::select_if(names(.) %in% c("sample", targets$sample)) %>%
      select(sample, everything())
    
    targets <- targets %>% dplyr::filter(sample %in% names(matrix)) %>%
      dplyr::mutate(id = paste(tissue, type, 1:nrow(.), sep = "_"), group = type) %>% 
      dplyr::rename(file = sample) %>% select(id, file, group)
    
    matrix <- matrix %>% dplyr::select(sample, targets$file) 
    colnames(matrix) <- c("gene_id", targets$id)
    
    cat(paste0("Matrices for ", type, " ready.\n"))
    return(list(targets = targets, matrix = matrix))
  }
  
  saveMatrix <- function(samples, tissue, type) {
    cat(paste0("Saving ", type, " ", tissue, " samples \n"))
    save(samples, file=paste0(RDATADIR, "/raw_",type, ".RData"), compress="xz")
    write_tsv(samples$matrix, paste0(DATADIR, "/", tissue, "-", type, "-matrix.tsv"))
    write_tsv(samples$targets, paste0(DATADIR,"/", tissue, "-", type, "-samples.tsv"))
  }
  
  normal_samples <- getSamples(TISSUE, NORMAL_TISSUE, "normal")
  cancer_samples <- getSamples(TISSUE, PRIMARY_DISEASE, "cancer")
  saveMatrix(normal_samples, TISSUE, "normal")
  saveMatrix(cancer_samples, TISSUE, "cancer")
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
  annot <- read_tsv(XENA_ANNOT,
                    col_names  = c("gene_id", "gene_name", "seqnames", "start", "end", "strand"), 
                    skip = 1)
  annot <- annot %>% dplyr::mutate(ensembl_id = unlist(lapply(strsplit(gene_id, "\\."), "[[", 1)), 
                            width = end - start)
  
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
  annot <- annot %>%
    inner_join(biomart %>% dplyr::select(ensembl_id, gc), by = "ensembl_id")
  
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
  
  save(full, file=paste(RDATADIR, "raw_full.RData", sep="/"), compress="xz")
  write_tsv(M, paste0(DATADIR, "/", TISSUE, "-matrix.tsv"))
  write_tsv(targets, paste0(DATADIR, "/", TISSUE, "-samples.tsv"))
  
  cat("raw_full.RData and matrices saved \n")
}
