def get_xena_dir(wildcards):
    return f'{config["datadir"]}/{config["xenadir"]}'

def get_tissue_dir(wildcards):
    return f'{config["datadir"]}/{wildcards.tissue}'

def get_output_files(wildcards):
    files = []
    for t in config["xena_tissues"]:
        files.append(config["datadir"]+"/"+t["name"]+"/rdata/mean10_protein_coding.RData")
    return files
        #files.append(config["datadir"]+ "/" + t + "/gc-upper_length-no_tmm_networks/cancer_network_1.adj")
        #files.append(config["datadir"]+ "/" + t + "/gc-upper_length-no_tmm_networks/normal_network_1.adj")
    
    #for t in config["tissues"]:
    #  files.append(config["datadir"]+"/"+t+"/"+t+"-matrix.tsv")

def getRawMatrixInput(wildcards):
    if wildcards.tissue not in config["xena"]:
        return [config["datadir"]+"/{wildcards.tissue}/manifests/{wildcards.tissue}-normal-rna_counts.txt", config["datadir"]+"/{wildcards.tissue}/manifests/{wildcards.tissue}-cancer-rna_counts.txt"]
    else:
        return [config["datadir"]+"/xena/counts.gz",config["datadir"]+"/xena/samples.txt.gz", config["datadir"]+"/xena/annot.tsv"]

def getNormalTissue(wildcards):
    return [x["normal"] for x in config["xena_tissues"] if x["name"] == wildcards.tissue][0]

def getCancerTissue(wildcards):
    return [x["cancer"] for x in config["xena_tissues"] if x["name"] == wildcards.tissue][0]
