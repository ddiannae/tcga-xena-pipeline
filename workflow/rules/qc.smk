rule qc:
    input:
        config["datadir"]+"/{tissue}/rdata/raw_full.RData"
    output:
        config["datadir"]+"/{tissue}/plots/pca_score_raw.png"
    params:
        tissue_dir=get_tissue_dir
    log:
        config["datadir"]+"/{tissue}/log/qc.log"
    script:
        "../scripts/QC.R"
