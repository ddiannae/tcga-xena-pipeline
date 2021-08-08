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
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")
## Get the Work and Data dir
library(data.table)
library(readr)
library(dplyr)
library(janitor)
library(rtracklayer)

TISSUE <- snakemake@params[["tissue"]]
TYPE <-  snakemake@params[["type"]]
IS_XENA <- snakemake@params[["is_xena"]]
MCCORES <- snakemake@params[["mccores"]]
# ###############################################################################
# ## Get and check count matrices
# ## 1. Find files
# ## 2. Read data
# ## 4. Check sample sizes
# ## 5. Check genes order
# ## 6. Build target matrix
# ###############################################################################

if(!IS_XENA) {
  
  RAWDIR <- paste(snakemake@params[["tissue_dir"]], "raw", sep="/")
  
  cat(paste0("Checking ", TYPE, " samples \n"))
     
  ## Find the files
  files_to_read <- list.files(path = paste0(RAWDIR,"/",  TYPE), 
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
  cat(paste0(TYPE, " samples have the same size \n"))
  
  ## Check genes order in samples
  genes <- do.call(cbind, lapply(all_files, function(x) dplyr::select(x, "gene_id")))
  genes <- t(unique(t(genes)))
  stopifnot(dim(genes)==c(size[1,1], 1))
  cat(paste0("Genes in ", TYPE, " samples match positions \n"))
  
  ## Build targets matrix
  targets <- data.frame(id = paste(TISSUE, TYPE, 1:length(files_to_read), sep = "_"),
                        file = unlist(lapply(strsplit(files_to_read, "/"), "[[", 10)),
                        file_id = unlist(lapply(strsplit(files_to_read, "/"), "[[", 9)),
                        group = TYPE, stringsAsFactors = FALSE)
  
  ## Rename columns in counts matrix
  matrix <- bind_cols(lapply(all_files, function(x) dplyr::select(x, "raw_counts")))
  colnames(matrix) <- targets$id
  
  matrix <- matrix %>% mutate(gene_id = genes[,1]) %>%
    dplyr::select(gene_id, everything())
  
} else {
  XENA_COUNTS <- snakemake@input[[1]]
  XENA_SAMPLES <- snakemake@input[[2]]
  PRIMARY <-snakemake@params[["primary"]]
  
  matrix_samples <- read_tsv(XENA_SAMPLES) %>% clean_names()
  matrix_counts <-  fread(XENA_COUNTS, nThread = MCCORES)
  
  extended_type <- ifelse(TYPE == "normal", "Normal Tissue", "Primary Tumor")
  camelTissue <- paste0(toupper(substring(TISSUE, 1, 1)), substring(TISSUE, 2))
    
  targets <- matrix_samples %>% dplyr::filter(primary_site == camelTissue &
                                                primary_disease_or_tissue == PRIMARY &
                                                sample_type == extended_type)
    
  ### Expected counts from XENA are in log2(expected_count+1)
  ### get them back to expected_count for downstream pipeline
  matrix <- matrix_counts %>% dplyr::select_if(names(.) %in% c("sample", targets$sample)) %>%
    dplyr::select(sample, everything())%>% dplyr::mutate(across(-sample, ~ .x^2-1))
  matrix[matrix < 0] <- 0
    
  targets <- targets %>% dplyr::filter(sample %in% names(matrix)) %>%
    dplyr::mutate(id = paste(TISSUE, TYPE, 1:nrow(.), sep = "_"), group = TYPE) %>% 
    dplyr::rename(file = sample) %>% select(id, file, group)
    
  matrix <- matrix %>% dplyr::select(sample, targets$file) 
  colnames(matrix) <- c("gene_id", targets$id)
}

cat(paste0("Matrices for ", TYPE, " ready.\n"))
cat("Saving matrix\n")
write_tsv(matrix, snakemake@output[["matrix"]])
cat("Saving samples\n")
write_tsv(targets, snakemake@output[["samples"]])