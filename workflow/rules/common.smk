import glob

def get_output_files(wildcards):
    files = []
    if config["end"] == "qc":
        for t in config["xena_tissues"]:
            files.append(config["datadir"]+"/"+t["name"]+"/rdata/raw.RData")
        for t in config["tissues"]:
            files.append(config["datadir"]+"/"+t+"/rdata/raw.RData")
    elif config["end"] == "norm_test":
        for t in config["xena_tissues"]:
            files.append(config["datadir"]+"/"+t["name"]+"/plots/normalization_plots.pdf")
        for t in config["tissues"]:
            files.append(config["datadir"]+"/"+t+"/plots/normalization_plots.pdf")
    return files

def get_expr_matrix(wildcards):
        return [file for file in glob.glob(f'{config["datadir"]}/{wildcards.tissue}/rdata/*_*_*_norm_cpm10_arsyn_{wildcards.type}.tsv')][0]


def get_xena_dir(wildcards):
    return f'{config["datadir"]}/{config["xenadir"]}'

def get_tissue_dir(wildcards):
    return f'{config["datadir"]}/{wildcards.tissue}'

def is_xena_tissue(wildcards):
    return wildcards.tissue in config["xena"]
    
def get_raw_matrix_input(wildcards):
    if is_xena_tissue(wildcards):
        xena_dir = get_xena_dir(wildcards)
        return [xena_dir +"/counts.gz", xena_dir +"/samples.txt.gz"]
    else:
        return f'{config["datadir"]}/{wildcards.tissue}/raw/{wildcards.type}/downloads_done.txt'
        
def get_annot_input(wildcards):
    preffix = f'{config["datadir"]}/{wildcards.tissue}/results/{wildcards.tissue}'
    input = {
    "normal_targets": f'{preffix}-normal-samples.tsv',
    "cancer_targets": f'{preffix}-cancer-samples.tsv',
    "normal_matrix":  f'{preffix}-normal-matrix.tsv',
    "cancer_matrix":  f'{preffix}-cancer-matrix.tsv'
    }
    
    if is_xena_tissue(wildcards):
        xena_dir = get_xena_dir(wildcards)
        input["xena_annot"] =  f'{xena_dir}/annot.tsv'
        return input
    else:
        return input

def get_xena_primary(wildcards):
    if is_xena_tissue(wildcards):
        if wildcards.type == "normal":
            return get_normal_tissue(wildcards)
        else:
            return get_cancer_tissue(wildcards)
    else:
        return ""

def get_normal_tissue(wildcards):
    return [x["normal"] for x in config["xena_tissues"] if x["name"] == wildcards.tissue][0]

def get_cancer_tissue(wildcards):
    return [x["cancer"] for x in config["xena_tissues"] if x["name"] == wildcards.tissue][0]
