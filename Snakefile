## Snakefile for ASCAT2 files from GDC
##
## Tissue type just like in GDC, lowercase is fine
#TISSUES = ["prostate", "pancreas", "bladder", "skin", "brain", "liver", "esophagus", "breast", "lung", "kidney", "colorectal", "uterus", "thyroid"]
TISSUES = ["esophagus"]
DATADIR ="/home/diana/Workspace/regulaciontrans-data"
FIGDIR = "figures"
files = [] 
for t in TISSUES:
  ## Example: data/breast/manifests/breast-cancer-rna_counts.txt"
    files.append(DATADIR+ "/" + t + "/plots/PCAScore_raw.png")
    files.append(DATADIR+ "/" + t + "/plots/normalization_plots.pdf")

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
rule normalization_test:
  input:
    DATADIR+"/{tissue}/rdata/raw_full.RData"
  output:
    DATADIR+"/{tissue}/plots/normalization_plots.pdf"
  shell:
    "Rscript src/normalizationTest.R {wildcards.tissue} {DATADIR}"

rule qc:
  input:
    DATADIR+"/{tissue}/rdata/raw_full.RData"
  output:
    DATADIR+"/{tissue}/plots/PCAScore_raw.png"
  shell:
    "Rscript src/QC.R {wildcards.tissue} {DATADIR}"
    
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
    Rscript src/getRawMatrices.R {wildcards.tissue} {DATADIR}
    """
rule get_manifest:
  output:
    ## Example: data/breast/manifests/breast-cancer-rna_counts.txt"
    DATADIR+"/{tissue}/manifests/{tissue}-{type}-rna_counts.txt",
  shell:
    """
    mkdir -p {DATADIR}/{wildcards.tissue}/manifests
    python src/queryGDC.py {wildcards.tissue} {wildcards.type} {DATADIR}/{wildcards.tissue}/manifests false 
    """
