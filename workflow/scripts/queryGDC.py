import requests
import json
import sys
import pandas as pd
from io import StringIO
import pathlib

logdir = snakemake.params[0] + "/log"
manifestdir = snakemake.params[0] + "/manifests"
pathlib.Path(logdir).mkdir(parents=True, exist_ok=True)
pathlib.Path(manifestdir).mkdir(parents=True, exist_ok=True)

sys.stderr = open(snakemake.log[0], "w")

primary_site = snakemake.params[1]
sample_type =snakemake.params[2]
sst = sample_type

if (sample_type == "cancer"):
    sample_type = "primary tumor"
elif (sample_type == "normal"):
    sample_type = "solid tissue normal"
else:
    sys.exit("Incorrect sample type")

print("Getting data for: " + primary_site + ", " + sample_type)

# Fields for the query
fields = [
    "cases.case_id",
    "file_name"
    ]

fields = ",".join(fields)

# Endpoints used
files_endpt = "https://api.gdc.cancer.gov/files"
manifest_endpt = "https://api.gdc.cancer.gov/manifest"

# Filters to get RNA seq data
filters = {
    "op": "and",
    "content":[
        {
        "op": "in",
        "content":{
            "field": "cases.project.primary_site",
            "value": [primary_site]
            }
        },
        {
        "op": "in",
        "content":{
            "field": "cases.samples.sample_type",
            "value": [sample_type]
            }
        },
        {
        "op": "in",
        "content":{
            "field": "files.analysis.workflow_type",
            "value": ["HTSeq - Counts"]
            }
        },
        {
        "op": "in",
        "content":{
            "field": "files.data_type",
            "value": ["Gene Expression Quantification"]  
            }
        }
    ]
}

params = {
    "filters": json.dumps(filters),
    "fields": fields,
    "format": "tsv",
    "size": "2000"
    }

# Getting data
response = requests.get(files_endpt, params = params)
data = response.content.decode("utf-8")

print("RNA-seq query done")

if(not data.strip()):
    sys.exit("No records found for " + primary_site + ", " + sample_type)

df = pd.read_csv(StringIO(data), sep ="\t")

df.columns = [col + "_rna" for col in df.columns] 
        
# Save file ids and file names just in case
df.to_csv(snakemake.output[0], sep="\t", index=False)  

# Getting RNASeq files manifest for future download
params = { "ids" : df["id_rna"].tolist() }
response = requests.post(manifest_endpt, data= json.dumps(params), headers = {"Content-Type":
    "application/json"})

print("Got RNASeq manifest")

f = open(snakemake.output[1], "w")
f.write(response.content.decode("utf-8"))
f.close()

print("RNASeq manifest written")
