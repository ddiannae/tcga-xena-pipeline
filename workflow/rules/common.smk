import glob

def get_output_files(wildcards):
    files = []
    if config["end"] == "qc":
        for t in config["tissues"]:
            files.append(config["datadir"]+"/"+t["name"]+"/rdata/raw.RData")
    elif config["end"] == "norm_test":
        for t in config["tissues"]:
            files.append(config["datadir"]+"/"+t["name"]+"/plots/normalization_plots.pdf")
    elif config["end"] == "correlation":
        for t in config["tissues"]:
            files.append(config["datadir"]+"/"+t["name"]+"/correlation/"+t["step1"]+"_"+t["step2"]+"_"+t["step3"]+"_si-arsyn_cancer_mi.adj")
            files.append(config["datadir"]+"/"+t["name"]+"/correlation/"+t["step1"]+"_"+t["step2"]+"_"+t["step3"]+"_si-arsyn_normal_mi.adj")
            #files.append(config["datadir"]+"/"+t["name"]+"/correlation/"+t["step1"]+"_"+t["step2"]+"_"+t["step3"]+"_si-arsyn_normal_pearson.tsv")
            #files.append(config["datadir"]+"/"+t["name"]+"/correlation/"+t["step1"]+"_"+t["step2"]+"_"+t["step3"]+"_si-arsyn_cancer_pearson.tsv")
    elif config["end"] == "deg":
        for t in config["tissues"]:
            files.append(config["datadir"]+"/"+t["name"]+"/deg/"+t["step1"]+"_"+t["step2"]+"_"+t["step3"]+"_si-arsyn_deg_results.tsv")
            print(config["datadir"]+"/"+t["name"]+"/deg/"+t["step1"]+"_"+t["step2"]+"_"+t["step3"]+"_si-arsyn_deg_results.tsv")
    return files

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
    return [x["normal"] for x in config["tissues"] if x["name"] == wildcards.tissue][0]

def get_cancer_tissue(wildcards):
    return [x["cancer"] for x in config["tissues"] if x["name"] == wildcards.tissue][0]

def get_xena_extended_type(wildcards):
    return [x[f'sample_type_{wildcards.type}'] if f'sample_type_{wildcards.type}' in x else None for x in config["tissues"] if x["name"] == wildcards.tissue][0]

def get_tissue_name(wildcards):
    return [x[f'tissue_name_{wildcards.type}'] if f'tissue_name_{wildcards.type}' in x else x["name"] for x in config["tissues"] if x["name"] == wildcards.tissue][0]

