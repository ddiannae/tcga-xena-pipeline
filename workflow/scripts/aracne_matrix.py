from MultiAracne import Aracne
import sys
import pathlib

sys.stdout = open(snakemake.log[0], "w")

exp_matrix = snakemake.input[0]
outmatrix = snakemake.output[0]
outdir = snakemake.params[0] + "/correlation/output_" + snakemake.params[1] + "_" + snakemake.params[2]
procs = int(snakemake.threads)

print(f"Saving files in {outdir}")
pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
    
ma = Aracne(exp_matrix)
ma.run(processes=procs, outdir=outdir, pval=1)
ma.join_adj(outdir, outdir+"/matrix.adj")
Aracne.adj_to_matrix(outdir+"/matrix.adj", outmatrix)
