rule qc:
    input:
        config["datadir"]+"/{tissue}/rdata/raw_full.RData"
    output:
        config["datadir"]+"/{tissue}/plots/{plots_type}/pca_score.png"
    params:
        tissue_dir=get_tissue_dir,
        plots_type="{plots_type}"
    log:
        config["datadir"]+"/{tissue}/log/{plots_type}_qc.log"
    script:
        "../scripts/NOISeqPlots.R"
