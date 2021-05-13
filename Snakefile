## Snakefile for ASCAT2 files from GDC
##
## Tissue type just like in GDC, lowercase is fine
### DONE_IMAGES = ["esophagus", "skin", "breast", "kidney","lung", "liver", "prostate", "pancreas", "bladder", "brain", "colorectal", "uterus", "thyroid"]  
#TISSUES = ["breast", "prostate", "pancreas", "bladder", "skin", "brain", "liver", "esophagus",  "lung", "kidney", "colorectal", "uterus", "thyroid"]
TISSUES = [ "prueba"]
DATADIR ="/datos/ot/diana/regulacion-trans"
# Adjust the number of cores according to the machine and number of tissues
MCCORES = 70
files = [] 
for t in TISSUES:
  ## Example: data/breast/manifests/breast-cancer-rna_counts.txt"
  files.append(DATADIR+ "/" + t + "/rdata/no_no_tmm_norm_cpm10_arsyn.RData")
  #files.append(DATADIR+ "/" + t + "/rdata/normalization_results.tsv")

rule all:
  input:
    files

#rule get_heatmap:
#  input: 
#    DATADIR+"/{tissue}/{tissue}-{type}-ascat-matrix.tsv"
#  output:
#    FIGDIR+"/{tissue}/{tissue}-{type}-ascat-heatmap.png"
#  shell:
#    """
#    mkdir -p {FIGDIR}/{wildcards.tissue}
#    Rscript src/getHeatmap.R {wildcards.tissue} {biomart} {input} {FIGDIR}/{wildcards.tissue}
#    """
# 

rule user_normalization:
  input:
    DATADIR+"/{tissue}/rdata/mean10_proteinCoding.RData"
  output:
    DATADIR+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_arsyn.RData"
  shell:
    "Rscript src/userNormalization.R {wildcards.tissue} {DATADIR} {wildcards.step1} {wildcards.step2} {wildcards.step3} > {DATADIR}/{wildcards.tissue}/{wildcards.step1}_{wildcards.step2}_{wildcards.step3}normalization_plots.log" 

rule normalization_plots:
  input:
    DATADIR+"/{tissue}/rdata/normalization_results.tsv"
  output:
    DATADIR+"/{tissue}/plots/normalization_plots.pdf"
  shell:
    "Rscript src/normalizationPlots.R {wildcards.tissue} {DATADIR} {MCCORES} > {DATADIR}/{wildcards.tissue}/normalization_plots.log" 

rule normalization_test:
  input:
    DATADIR+"/{tissue}/rdata/mean10_proteinCoding.RData",
  output:
    DATADIR+"/{tissue}/rdata/normalization_results.tsv"
  shell:
    "Rscript src/normalizationTest.R {wildcards.tissue} {DATADIR} {MCCORES} > {DATADIR}/{wildcards.tissue}/normalization_test.log"

rule filter_low_expressin:
  input:
    DATADIR+"/{tissue}/rdata/raw_full.RData",
    DATADIR+"/{tissue}/plots/pca_score_raw.png"
  output:
    DATADIR+"/{tissue}/rdata/mean10_proteinCoding.RData",
  shell:
    "Rscript src/filterLowExpression.R {wildcards.tissue} {DATADIR} > {DATADIR}/{wildcards.tissue}/filter_low.log"

rule qc:
  input:
    DATADIR+"/{tissue}/rdata/raw_full.RData"
  output:
    DATADIR+"/{tissue}/plots/pca_score_raw.png"
  shell:
    "Rscript src/QC.R {wildcards.tissue} {DATADIR} > {DATADIR}/{wildcards.tissue}/qc.log"
    
## We need to run these two together because the output of the download_files
## tasks depends on the manifest and there is no easy way to specify this on 
## snakemake
rule download_files_and_get_ascat_matrix:
  input:
    ## Manifest file
    normal=DATADIR+"/{tissue}/manifests/{tissue}-normal-rna_counts.txt",
    cancer=DATADIR+"/{tissue}/manifests/{tissue}-cancer-rna_counts.txt",
  output: 
    DATADIR+"/{tissue}/{tissue}-matrix.tsv",
    DATADIR+"/{tissue}/{tissue}-samples.tsv",
    DATADIR+"/{tissue}/rdata/raw_full.RData"
  shell:
    """
    mkdir -p {DATADIR}/{wildcards.tissue}/raw/{wildcards.tissue}-normal-rna
    mkdir -p {DATADIR}/{wildcards.tissue}/raw/{wildcards.tissue}-cancer-rna
    ./bin/gdc-client download -d {DATADIR}/{wildcards.tissue}/raw/{wildcards.tissue}-normal-rna -m {input.normal} --retry-amount 3
    ./bin/gdc-client download -d {DATADIR}/{wildcards.tissue}/raw/{wildcards.tissue}-cancer-rna -m {input.cancer} --retry-amount 3
    Rscript src/getRawMatrices.R {wildcards.tissue} {DATADIR} > {DATADIR}/{wildcards.tissue}/get-matrix.log
    """
rule get_manifest:
  output:
    ## Example: data/breast/manifests/breast-cancer-rna_counts.txt"
    DATADIR+"/{tissue}/manifests/{tissue}-{type}-rna_counts.txt"
  shell:
    """
    mkdir -p {DATADIR}/{wildcards.tissue}/manifests
    python src/queryGDC.py {wildcards.tissue} {wildcards.type} {DATADIR}/{wildcards.tissue}/manifests false 
    """
