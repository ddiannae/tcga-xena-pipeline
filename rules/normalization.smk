rule normalization_plots:
    input:
        config["datadir"]+"/{tissue}/rdata/normalization_results.tsv"
    output:
        config["datadir"]+"/{tissue}/plots/normalization_plots.pdf"
    params:
        tissue_dir=get_tissue_dir
    log:
        config['datadir']+"/{tissue}/log/normalization_plots.log"
    script:
        "../scripts/normalizationPlots.R"

rule normalization_test:
    input:
        config["datadir"]+"/{tissue}/rdata/mean10_protein_coding.RData",
    output:
        config["datadir"]+"/{tissue}/rdata/normalization_results.tsv"
    params:
        tissue_dir=get_tissue_dir,
        mccores=config["mccores"]
    log:
        config['datadir']+"/{tissue}/log/normalization_test.log"
    script:
        "../scripts/normalizationTest.R"

rule filter_low_expression:
    input:
        config["datadir"]+"/{tissue}/rdata/raw_full.RData",
        config["datadir"]+"/{tissue}/plots/pca_score_raw.png"
    output:
        config["datadir"]+"/{tissue}/rdata/mean10_protein_coding.RData"
    log:
        config["datadir"]+"/{tissue}/log/filter_low.log"
    script:
        "../scripts/filterLowExpression.R"
