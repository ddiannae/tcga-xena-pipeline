rule normalization_plots:
    input:
        config["datadir"]+"/{tissue}/results/normalization_results.tsv"
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
        config["datadir"]+"/{tissue}/rdata/raw.RData",
    output:
        results=config["datadir"]+"/{tissue}/results/normalization_results.tsv",
        gc_norms=config["datadir"]+"/{tissue}/rdata/gc_norms.RData",
        ln_norms=config["datadir"]+"/{tissue}/rdata/ln_norms.RData"
    params:
        tissue_dir=get_tissue_dir,
    threads: 25
    log:
        config['datadir']+"/{tissue}/log/normalization_test.log"
    script:
        "../scripts/normalizationTest.R"
 
rule user_normalization:
    input:
        config["datadir"]+"/{tissue}/rdata/raw.RData"
    output:
        norm_rdata=config["datadir"]+"/{tissue}/rdata/{step1}_{step2}_{step3}_no-arsyn_full.RData",
        normal_matrix=config["datadir"]+"/{tissue}/results/{step1}_{step2}_{step3}_no-arsyn_normal.tsv",
        cancer_matrix=config["datadir"]+"/{tissue}/results/{step1}_{step2}_{step3}_no-arsyn_cancer.tsv",
        gene_list=config["datadir"]+"/{tissue}/results/{step1}_{step2}_{step3}_no-arsyn_genelist.txt"
    params:
        tissue_dir=get_tissue_dir,
        step1="{step1}",
        step2="{step2}",
        step3="{step3}",
        xena=is_xena_tissue
    log:
        config['datadir']+"/{tissue}/log/{step1}_{step2}_{step3}_normalization.log"
    script:
        "../scripts/userNormalization.R"

rule arsyn:
    input:
        config["datadir"]+"/{tissue}/rdata/{step1}_{step2}_{step3}_no-arsyn_full.RData",
        config["datadir"]+"/{tissue}/plots/{step1}_{step2}_{step3}_no-arsyn/density.png"
    output:
        arsyn_rdata=config["datadir"]+"/{tissue}/rdata/{step1}_{step2}_{step3}_si-arsyn_full.RData",
        normal_matrix=config["datadir"]+"/{tissue}/results/{step1}_{step2}_{step3}_si-arsyn_normal.tsv",
        cancer_matrix=config["datadir"]+"/{tissue}/results/{step1}_{step2}_{step3}_si-arsyn_cancer.tsv",
        gene_list=config["datadir"]+"/{tissue}/results/{step1}_{step2}_{step3}_si-arsyn_genelist.txt"
    params:
        tissue_dir=get_tissue_dir,
        step1="{step1}",
        step2="{step2}",
        step3="{step3}"
    log:
        config['datadir']+"/{tissue}/log/{step1}_{step2}_{step3}_arsyn.log"
    script:
        "../scripts/runArsyn.R"
