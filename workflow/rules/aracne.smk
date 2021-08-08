rule aracne:
    input:
        config["datadir"]+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_{arsyn}{type}.tsv", 
    output:
        config["datadir"]+"/{tissue}/networks/{step1}_{step2}_{step3}_{arsyn}{type}_network.adj"
    singularity:
        config["aracne_singularity"]
    params:
        get_tissue_dir
    threads: 78
    log:
        config['datadir']+"/{tissue}/log/{step1}_{step2}_{step3}{arsyn}_{type}_aracne.log"
    script:
        "../scripts/aracne_matrix.py"
