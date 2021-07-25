
configfile: "config.yaml"

include: "rules/common.smk"
include: "rules/xena.smk"

# Adjust the umber of cores according to the machine and number of tissues
MCCORES = 78
files = [] 
for t in config["xena_tissues"]:
	files.append(config["datadir"]+"/"+t["name"]+"/plots/pca_score_raw.png")
	#files.append(config["datadir"]+ "/" + t + "/gc-upper_length-no_tmm_networks/cancer_network_1.adj")
	#files.append(config["datadir"]+ "/" + t + "/gc-upper_length-no_tmm_networks/normal_network_1.adj")

#for t in config["tissues"]:
#  files.append(config["datadir"]+"/"+t+"/"+t+"-matrix.tsv")

def getRawMatrixInput(wildcards):
	print(wildcards.tissue)
	if wildcards.tissue not in config["xena"]:
		return [config["datadir"]+"/{wildcards.tissue}/manifests/{wildcards.tissue}-normal-rna_counts.txt", config["datadir"]+"/{wildcards.tissue}/manifests/{wildcards.tissue}-cancer-rna_counts.txt"]
	else:
		return [config["datadir"]+"/xena/counts.gz",config["datadir"]+"/xena/samples.txt.gz", config["datadir"]+"/xena/annot.tsv"]

def getNormalTissue(wildcards):
	return [x["normal"] for x in config["xena_tissues"] if x["name"] == wildcards.tissue][0]

def getCancerTissue(wildcards):
	return [x["cancer"] for x in config["xena_tissues"] if x["name"] == wildcards.tissue][0]

rule all:
	input:
		files

rule run_infotheo:
	input: 
		config["datadir"]+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_arsyn_{type}.tsv"
	output:
		config["datadir"]+"/{tissue}/{step1}_{step2}_{step3}_mi/{type}_mi_matrix.adj"
	shell:
		"""
		mkdir -p {config['datadir']}/{wildcards.tissue}/{wildcards.step1}_{wildcards.step2}_{wildcards.step3}_mi
		Rscript src/getMIMatrix.R {input} {output} {MCCORES}
		"""

rule get_sif:
	input:
		config["datadir"]+"/{tissue}/{step1}_{step2}_{step3}_networks/{type}_network_{pval}.adj"
	output:
		config["datadir"]+"/{tissue}/{step1}_{step2}_{step3}_networks/{type}_network_{pval}.sif"
	shell:
		"""
		python src/python/adj2sif.py {input} {output}
		"""

rule get_adj:
	input:
		config["datadir"]+"/{tissue}/{step1}_{step2}_{step3}_networks/{type}_network_{pval}.aracne_adj"
	output:
		config["datadir"]+"/{tissue}/{step1}_{step2}_{step3}_networks/{type}_network_{pval}.adj"
	shell:
		"""
		Rscript src/aracneAdjToAdj.R {input} {output} {MCCORES}
		"""

rule get_adj_aracne_matrix:
	input:
		config["datadir"]+"/{tissue}/{step1}_{step2}_{step3}_networks/{type}_adj_{pval}"
	output:
		config["datadir"]+"/{tissue}/{step1}_{step2}_{step3}_networks/{type}_network_{pval}.aracne_adj"
	shell:
		"""
		python src/python/joinadj.py {input} {output}
		"""

rule run_aracne_all:
	input: 
		matrix=config["datadir"]+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_arsyn_{type}.tsv",
		genelist=config["datadir"]+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_genelist.txt"
	output:
		directory(config["datadir"]+"/{tissue}/{step1}_{step2}_{step3}_networks/{type}_adj_{pval}")
	shell:
		"""
		export ARACNEHOME={config['aracnehome']}
		mkdir -p {output}
		python src/python/aracne-par.py {input.matrix} {input.genelist} {MCCORES} {wildcards.pval} {output} > {config['datadir']}/{wildcards.tissue}/log/aracne_{wildcards.type}_{wildcards.pval}.log
		"""

rule user_normalization:
	input:
		config["datadir"]+"/{tissue}/rdata/mean10_proteinCoding.RData"
	output:
		config["datadir"]+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_arsyn.RData",
		config["datadir"]+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_arsyn_cancer.tsv",
		config["datadir"]+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_arsyn_normal.tsv",
		config["datadir"]+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_genelist.txt"
	shell:
		"Rscript src/userNormalization.R {wildcards.tissue} {config['datadir']} {wildcards.step1} {wildcards.step2} {wildcards.step3} > {config['datadir']}/{wildcards.tissue}/log/{wildcards.step1}_{wildcards.step2}_{wildcards.step3}normalization_plots.log" 

rule normalization_plots:
	input:
		config["datadir"]+"/{tissue}/rdata/normalization_results.tsv"
	output:
		config["datadir"]+"/{tissue}/plots/normalization_plots.pdf"
	shell:
		"Rscript src/normalizationPlots.R {wildcards.tissue} {config['datadir']} {MCCORES} > {config['datadir']}/{wildcards.tissue}/log/normalization_plots.log" 

rule normalization_test:
	input:
		config["datadir"]+"/{tissue}/rdata/mean10_proteinCoding.RData",
	output:
		config["datadir"]+"/{tissue}/rdata/normalization_results.tsv"
	shell:
		"Rscript src/normalizationTest.R {wildcards.tissue} {config['datadir']} {MCCORES} > {config['datadir']}/{wildcards.tissue}/log/normalization_test.log"

rule filter_low_expression:
	input:
		config["datadir"]+"/{tissue}/rdata/raw_full.RData",
		config["datadir"]+"/{tissue}/plots/pca_score_raw.png"
	output:
		config["datadir"]+"/{tissue}/rdata/mean10_proteinCoding.RData",
	params:
		tissuedir=config["datadir"]+"/{wildcards.tissue}",
		logfile=config["datadir"]+"/{wildcards.tissue}/log/filter_low.log"
	run:
		shell(f'Rscript src/filterLowExpression.R {wildcards.tissue} {params.tissuedir} > {params.logfile}')

rule qc:
	input:
		config["datadir"]+"/{tissue}/rdata/raw_full.RData"
	output:
		config["datadir"]+"/{tissue}/plots/pca_score_raw.png"
	params:
		tissuedir=config["datadir"]+"/{wildcards.tissue}",
		logfile=config["datadir"]+"/{wildcards.tissue}/log/qc.log"
	run:
		shell(f'Rscript src/QC.R {wildcards.tissue} {params.tissuedir} > {params.logfile}')
		
## We need to run these two together because the output of the download_files
## tasks depends on the manifest and there is no easy way to specify this on 
## snakemake
rule download_files_and_get_ascat_matrix:
	input: getRawMatrixInput
	output: 
		config["datadir"]+"/{tissue}/{tissue}-matrix.tsv",
		config["datadir"]+"/{tissue}/{tissue}-samples.tsv",
		config["datadir"]+"/{tissue}/rdata/raw_full.RData"
	params:
		normaldir=config["datadir"]+"/{wildcards.tissue}/raw/{wildcards.tissue}-normal-rna",
		cancerdir=config["datadir"]+"/{wildcards.tissue}/raw/{wildcards.tissue}-cancer-rna",
		tissuedir=config["datadir"]+"/{wildcards.tissue}",
		logdir=config["datadir"]+"/{wildcards.tissue}/log",
		logfile=config["datadir"]+"/{wildcards.tissue}/log/get-matrix.log",
		normaltissue=getNormalTissue,
		cancertissue=getCancerTissue
	run: 
		if wildcards.tissue not in config["xena"]:
			shell(f'mkdir -p {params.normaldir}')
			shell(f'mkdir -p {params.cancerdir}')
			shell(f'./bin/gdc-client download -d {params.normaldir} -m {input[0]} --retry-amount 3')
			shell(f'./bin/gdc-client download -d {params.cancerdir} -m {input[1]} --retry-amount 3')
			shell(f'Rscript src/getRawMatrices.R {wildcards.tissue} {params.tissuedir} > {params.logfile}')
		else:
			shell(f'mkdir -p {params.tissuedir}')
			shell(f'mkdir -p {params.logdir}')
			shell(f'Rscript src/getRawXENAMatrices.R {input[0]} {input[1]} {input[2]} {wildcards.tissue} "{params.normaltissue}" "{params.cancertissue}" {params.tissuedir} {config["mccores"]} > {params.logfile}') 

rule get_manifest:
	output:
		## Example: data/breast/manifests/breast-cancer-rna_counts.txt"
		config["datadir"]+"/{tissue}/manifests/{tissue}-{type}-rna_counts.txt"
	params:
		logdir=config["datadir"]+"/{wildcards.tissue}/log",
		manidir=config["datadir"]+"/{wildcards.tissue}/manifests"
	shell:
		"""
		mkdir -p {params.logdir}
		mkdir -p {params.manidir}
		python src/queryGDC.py {wildcards.tissue} {wildcards.type} {params.manidir} false 
		"""
