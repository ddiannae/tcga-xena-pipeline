rule get_gdc_files:
    input:
        config["datadir"]+"/{tissue}/manifests/{type}-rna_counts.txt"
    output:
        directory(config["datadir"]+"/{tissue}/raw/{type}"),
        config["datadir"]+"/{tissue}/raw/{type}/downloads_done.txt"
    shell:
        """
        mkdir -p {output[0]};
        ./bin/gdc-client download -d {output[0]} -m {input[0]} --retry-amount 3;
        touch {output[1]}
        """
    
rule get_manifest:
    output:
        ## Example: data/breast/manifests/breast-cancer-rna_counts.txt"
        config["datadir"]+"/{tissue}/manifests/{type}-files.tsv",
        config["datadir"]+"/{tissue}/manifests/{type}-rna_counts.txt",
    params:
        tissue_dir=get_tissue_dir,
        tissue="{tissue}",
        type="{type}"
    log:
        config['datadir']+"/{tissue}/log/query_gdc_{type}.log"
    script:
        "../scripts/queryGDC.py"
