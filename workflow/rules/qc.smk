rule qc:
    input:
        config["datadir"]+"/{tissue}/rdata/raw_full.RData"
    output:
        config["datadir"]+"/{tissue}/plots/{plots_type}/rna_composition.png"
    params:
        tissue_dir=get_tissue_dir,
        plots_type="{plots_type}"
    log:
        config["datadir"]+"/{tissue}/log/{plots_type}_qc.log"
    script:
        "../scripts/NOISeqPlots.R"

rule pca:
    input:
        config["datadir"]+"/{tissue}/rdata/{plots_type}_full.RData"
    output:
        score=config["datadir"]+"/{tissue}/plots/{plots_type}/pca_score.png",
        loading=config["datadir"]+"/{tissue}/plots/{plots_type}/pca_loading.png",
        variance=config["datadir"]+"/{tissue}/plots/{plots_type}/pca_variance.png"
    params:
        tissue_dir=get_tissue_dir,
        plots_type="{plots_type}"
    log:
        config["datadir"]+"/{tissue}/log/{plots_type}_pca.log"
    script:
        "../scripts/PCA.R"
