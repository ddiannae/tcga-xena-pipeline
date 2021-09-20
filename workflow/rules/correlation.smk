rule aracne:
    input:
        expr_matrix=config["datadir"]+"/{tissue}/results/{step1}_{step2}_{step3}_{arsyn}_{type}.tsv",
        density=config["datadir"]+"/{tissue}/plots/{step1}_{step2}_{step3}_{arsyn}/density.png"
    output:
        config["datadir"]+"/{tissue}/correlation/{step1}_{step2}_{step3}_{arsyn}_{type}_mi.adj"
    singularity:
        config["aracne_singularity"]
    params:
        get_tissue_dir,
        "{type}", 
        "{arsyn}"
    threads: 39
    log:
        config["datadir"]+"/{tissue}/log/{step1}_{step2}_{step3}_{arsyn}_{type}_aracne.log"
    script:
        "../scripts/aracne_matrix.py"

rule get_pearson_matrix:
    input: 
        config["datadir"]+"/{tissue}/results/{step1}_{step2}_{step3}_{arsyn}_{type}.tsv", 
        annot=config["datadir"]+"/{tissue}/rdata/annot.RData"
    output:
        plot=config["datadir"]+"/{tissue}/correlation/{step1}_{step2}_{step3}_{arsyn}_{type}_heatmap.png",
        csv=config["datadir"]+"/{tissue}/correlation/{step1}_{step2}_{step3}_{arsyn}_{type}_pearson.tsv"
    log:
        config['datadir']+"/{tissue}/log/{step1}_{step2}_{step3}_{arsyn}_{type}_pearson.log"
    script:
        "../scripts/getPearsonMatrix.R"
