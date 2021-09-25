rule get_deg:
    input:
        full=config["datadir"]+"/{tissue}/rdata/{step1}_{step2}_{step3}_{arsyn}_full.RData",
        annot=config["datadir"]+"/{tissue}/rdata/annot.RData"
    output:
        deg_results=config["datadir"]+"/{tissue}/deg/{step1}_{step2}_{step3}_{arsyn}_deg_results.tsv"
    params:
        tissue="{tissue}",
        deg_dir=config["datadir"]+"/{tissue}/deg"
    log:
        config['datadir']+"/{tissue}/log/{step1}_{step2}_{step3}_{arsyn}_deg.log"
    script:
        "../scripts/deg.R"
     
