rule aracne:
    input:
        config["datadir"]+"/{tissue}/results/{step1}_{step2}_{step3}_norm_{arsyn}{type}.tsv" 
    output:
        config["datadir"]+"/{tissue}/correlation/{step1}_{step2}_{step3}_{arsyn}{type}_network.adj"
    singularity:
        config["aracne_singularity"]
    params:
        get_tissue_dir
    threads: 39
    log:
        config['datadir']+"/{tissue}/log/{step1}_{step2}_{step3}{arsyn}_{type}_aracne.log"
    script:
        "../scripts/aracne_matrix.py"

rule get_pearson_matrix:
    input: 
        config["datadir"]+"/{tissue}/results/{step1}_{step2}_{step3}_norm_{arsyn}{type}.tsv", 
        annot=config["datadir"]+"/{tissue}/rdata/annot.RData"
    output:
        plot=config["datadir"]+"/{tissue}/correlation/{step1}_{step2}_{step3}_norm_{arsyn}{type}_heatmap.png",
        csv=config["datadir"]+"/{tissue}/correlation/{step1}_{step2}_{step3}_norm_{arsyn}{type}_pearson.tsv"
    log:
        config['datadir']+"/{tissue}/log/{step1}_{step2}_{step3}{arsyn}_{type}_pearson.log"
    script:
        "../scripts/getPearsonMatrix.R"
