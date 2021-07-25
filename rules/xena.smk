rule get_xena_counts:
    output:
        config["datadir"]+"/"+config["xenadir"]+"/counts.gz"
    params:
        xenadir=get_xena_dir
    shell:
        "mkdir -p {params.xenadir};"
        "wget https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_gene_expected_count.gz -O {output}"

rule get_xena_samples:
    output:
        config["datadir"]+"/"+config["xenadir"]+"/samples.txt.gz"
    params:
        xenadir=get_xena_dir
    shell:
        "mkdir -p {params.xenadir};"
        "wget https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGTEX_phenotype.txt.gz -O {output}"

rule get_xena_annotations:
    output:
        config["datadir"]+"/"+config["xenadir"]+"/annot.tsv"
    params:
        xenadir=get_xena_dir
    shell:
        "mkdir -p {params.xenadir};"
        "wget https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/probeMap%2Fgencode.v23.annotation.gene.probemap -O {output}"

