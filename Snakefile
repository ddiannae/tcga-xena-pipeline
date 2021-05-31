## Snakefile for ASCAT2 files from GDC
##
## Tissue type just like in GDC, lowercase is fine
### DONE_IMAGES = ["esophagus", "skin", "breast", "kidney","lung", "liver", "prostate", "pancreas", "bladder", "brain", "colorectal", "uterus", "thyroid"]  
#
#TISSUES = ["breast", "prostate", "pancreas", "bladder", "skin", "brain", "liver", "esophagus",  "lung", "kidney", "colorectal", "uterus", "thyroid"]
#"bladder", "brain", "breast", "colorectal", "esophagus", "kidney", "liver", "lung", "pancreas", "prostate", "skin", "thyroid", "uterus"
TISSUES = ["bladder"]
DATADIR ="/datos/ot/diana/regulacion-trans"
ARACNEHOME ="./bin/ARACNE"

# Adjust the umber of cores according to the machine and number of tissues
MCCORES = 70
files = [] 
for t in TISSUES:
  ## Example:data/breast/manifests/breast-cancer-rna_counts.txt"
  files.append(DATADIR+ "/" + t + "/gc-full_length-full_tmm_networks/normal_network_1e-8.sort")
  files.append(DATADIR+ "/" + t + "/gc-full_length-full_tmm_networks/cancer_network_1e-8.sort")
  #files.append(DATADIR+"/" + t+ "/manifests/"+t+"-cancer-rna_counts.txt")
  #files.append(DATADIR+"/" + t+ "/manifests/"+t+"-normal-rna_counts.txt")
  #files.append(DATADIR+"/" + t + "/plots/normalization_plots.pdf")
  #files.apped(DATADIR+ "/" + t + "/rdata/normalization_results.tsv")

rule all:
  input:
    files

    
rule get_sif_sorted:
  input:
    DATADIR+"/{tissue}/{step1}_{step2}_{step3}_networks/{type}_network_{pval}.sif",
  output:
    DATADIR+"/{tissue}/{step1}_{step2}_{step3}_networks/{type}_network_{pval}.sort",
  shell:
    """
    sort -r -k3,3 {input} > {output}
    """

rule get_sif:
  input:
    DATADIR+"/{tissue}/{step1}_{step2}_{step3}_networks/{type}_network_{pval}.adj"
  output:
    DATADIR+"/{tissue}/{step1}_{step2}_{step3}_networks/{type}_network_{pval}.sif"
  shell:
    """
    python src/python/adj2sif.py {input} > {output}
    """

rule get_adj_matrix:
  input:
    DATADIR+"/{tissue}/{step1}_{step2}_{step3}_networks/{type}_adj_{pval}",
  output:
    DATADIR+"/{tissue}/{step1}_{step2}_{step3}_networks/{type}_network_{pval}.adj",
  shell:
    """
    python src/python/joinadj.py {input} {output}
    """

rule run_aracne_all:
  input: 
    cancer=DATADIR+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_arsyn_cancer.tsv",
    normal=DATADIR+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_arsyn_normal.tsv",
    genelist=DATADIR+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_genelist.txt"
  output:
    normal=directory(DATADIR+"/{tissue}/{step1}_{step2}_{step3}_networks/normal_adj_{pval}"),
    cancer=directory(DATADIR+"/{tissue}/{step1}_{step2}_{step3}_networks/cancer_adj_{pval}")
  shell:
    """
    export ARACNEHOME={ARACNEHOME}
    mkdir -p {output.normal}
    python src/python/aracne-par.py {input.normal} {input.genelist} {MCCORES} {wildcards.pval} {output.normal} > {DATADIR}/{wildcards.tissue}/log/aracne_normal.log
    mkdir -p {output.cancer}
    python src/python/aracne-par.py {input.cancer} {input.genelist} {MCCORES} {wildcards.pval} {output.cancer} > {DATADIR}/{wildcards.tissue}/log/aracne_cancer.log
    """

rule user_normalization:
  input:
    DATADIR+"/{tissue}/rdata/mean10_proteinCoding.RData"
  output:
    DATADIR+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_arsyn.RData",
    DATADIR+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_arsyn_cancer.tsv",
    DATADIR+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_arsyn_normal.tsv",
    DATADIR+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_genelist.txt"
  shell:
    "Rscript src/userNormalization.R {wildcards.tissue} {DATADIR} {wildcards.step1} {wildcards.step2} {wildcards.step3} > {DATADIR}/{wildcards.tissue}/log/{wildcards.step1}_{wildcards.step2}_{wildcards.step3}normalization_plots.log" 

rule normalization_plots:
  input:
    DATADIR+"/{tissue}/rdata/normalization_results.tsv"
  output:
    DATADIR+"/{tissue}/plots/normalization_plots.pdf"
  shell:
    "Rscript src/normalizationPlots.R {wildcards.tissue} {DATADIR} {MCCORES} > {DATADIR}/{wildcards.tissue}/log/normalization_plots.log" 

rule normalization_test:
  input:
    DATADIR+"/{tissue}/rdata/mean10_proteinCoding.RData",
  output:
    DATADIR+"/{tissue}/rdata/normalization_results.tsv"
  shell:
    "Rscript src/normalizationTest.R {wildcards.tissue} {DATADIR} {MCCORES} > {DATADIR}/{wildcards.tissue}/log/normalization_test.log"

rule filter_low_expression:
  input:
    DATADIR+"/{tissue}/rdata/raw_full.RData",
    DATADIR+"/{tissue}/plots/pca_score_raw.png"
  output:
    DATADIR+"/{tissue}/rdata/mean10_proteinCoding.RData",
  shell:
    "Rscript src/filterLowExpression.R {wildcards.tissue} {DATADIR} > {DATADIR}/{wildcards.tissue}/log/filter_low.log"

rule qc:
  input:
    DATADIR+"/{tissue}/rdata/raw_full.RData"
  output:
    DATADIR+"/{tissue}/plots/pca_score_raw.png"
  shell:
    "Rscript src/QC.R {wildcards.tissue} {DATADIR} > {DATADIR}/{wildcards.tissue}/log/qc.log"
    
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
    Rscript src/getRawMatrices.R {wildcards.tissue} {DATADIR} > {DATADIR}/{wildcards.tissue}/log/get-matrix.log
    """
rule get_manifest:
  output:
    ## Example: data/breast/manifests/breast-cancer-rna_counts.txt"
    DATADIR+"/{tissue}/manifests/{tissue}-{type}-rna_counts.txt"
  shell:
    """
    mkdir -p {DATADIR}/{wildcards.tissue}/log
    mkdir -p {DATADIR}/{wildcards.tissue}/manifests
    python src/queryGDC.py {wildcards.tissue} {wildcards.type} {DATADIR}/{wildcards.tissue}/manifests false 
    """
