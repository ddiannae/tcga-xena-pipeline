rule aracne:
    input:
        config["datadir"]+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_arsyn_{type}.tsv", 
    output:
        config["datadir"]+"/{tissue}/networks/{step1}_{step2}_{step3}_{type}_network.adj"
    singularity:
        config["aracne_singularity"]
    params:
        get_tissue_dir,
        config["mccores"]
    log:
        config['datadir']+"/{tissue}/log/{step1}_{step2}_{step3}_{type}_aracne.log"
    script:
        "../scripts/aracne_matrix.py"
