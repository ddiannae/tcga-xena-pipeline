def get_xena_dir(wildcards):
    return f'{config["datadir"]}/{config["xenadir"]}'

def get_tissue_dir(wildcards):
    return f'{config["datadir"]}/{wildcards.tissue}'
