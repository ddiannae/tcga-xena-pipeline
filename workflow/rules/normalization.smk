rule user_normalization:
    input:
        config["datadir"]+"/{tissue}/rdata/mean10_protein_coding.RData"
    output:
        norm_cpm10=config["datadir"]+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10.RData",
        normal_matrix=config["datadir"]+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_normal.tsv",
        cancer_matrix=config["datadir"]+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_cancer.tsv",
        gene_list=config["datadir"]+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_genelist.txt"
    params:
        tissue_dir=get_tissue_dir,
        step1="{step1}",
        step2="{step2}",
        step3="{step3}"
    log:
        config['datadir']+"/{tissue}/log/{step1}_{step2}_{step3}_normalization.log"
    script:
        "../scripts/userNormalization.R"

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
        config["datadir"]+"/{tissue}/plots/raw/pca_score.png"
    output:
        config["datadir"]+"/{tissue}/rdata/mean10_protein_coding.RData"
    log:
        config["datadir"]+"/{tissue}/log/filter_low.log"
    script:
        "../scripts/filterLowExpression.R"
