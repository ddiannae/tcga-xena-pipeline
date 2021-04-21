import requests
import json
import sys
import pandas as pd
from io import StringIO

primary_site = sys.argv[1]
sample_type = sys.argv[2]
mdir = sys.argv[3]

if(len(sys.argv) > 4):
    ascat = json.loads(sys.argv[4]);
else:
    ascat = False;

sst = sample_type

if (sample_type == "tumor"):
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

if(ascat and sst == "tumor"):
    df_rna = pd.read_csv(StringIO(data), sep ="\t")
    df_rna.drop_duplicates(subset=["cases.0.case_id"], inplace=True)
    
    # Change filters to get ASCAT cnvs data
    filters["content"][2]["content"]["value"] = "ASCAT2"
    filters["content"][3]["content"]["value"] = "Gene Level Copy Number"
    
    params = {
        "filters": json.dumps(filters),
        "fields": fields,
        "format": "tsv",
        "size": "2000"
    }
    
    response = requests.get(files_endpt, params = params)
    data = response.content.decode("utf-8")
    
    print("ASCAT query done")
    
    df_ascat = pd.read_csv(StringIO(data), sep ="\t")
    df_ascat.drop_duplicates(subset=["cases.0.case_id"], inplace=True)
    
    # We need cases that have both files
    df = pd.merge(df_ascat, df_rna, "inner", on="cases.0.case_id",
            suffixes=("_ascat", "_rna"))
     
    # Getting ASCAT2 files manifest for future download
    params = { "ids" : df["id_ascat"].tolist() }
    response = requests.post(manifest_endpt, data= json.dumps(params), headers = {"Content-Type":
        "application/json"})
    
    print("Got ASCAT2 manifest")
    
    f = open(mdir + "/" + "-".join([primary_site.lower(),
        sst.lower(), "ascat.txt"]), "w")
    f.write(response.content.decode("utf-8"))
    f.close()
    
    print("ASCAT2 manifest written")
            
else:
    print("No ASCAT2 query")
    if(not data.strip()):
        sys.exit("No records found for " + primary_site + ", " + sample_type)
    
    df = pd.read_csv(StringIO(data), sep ="\t")
    df.drop_duplicates(subset=["cases.0.case_id"], inplace=True)

    df.columns = [col + "_rna" for col in df.columns] 
        
# Save file ids and file names just in case
df.to_csv(mdir + "/" + "-".join([primary_site.lower(), sst.lower()]) + "-files.tsv", sep="\t", index=False)  

# Getting RNASeq files manifest for future download
params = { "ids" : df["id_rna"].tolist() }
response = requests.post(manifest_endpt, data= json.dumps(params), headers = {"Content-Type":
    "application/json"})

print("Got RNASeq manifest")

f = open(mdir + "/" + "-".join([primary_site.lower(),
    sst.lower(), "rna_counts.txt"]), "w")
f.write(response.content.decode("utf-8"))
f.close()

print("RNASeq manifest written")
