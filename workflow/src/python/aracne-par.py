import pandas as pd
from multiprocessing import Pool
from functools import partial
import os, sys


def f(x, fexp, outd, pval):
    cmd = os.environ.get('ARACNEHOME') + "/aracne2 -H " + os.environ.get('ARACNEHOME')  + " -i " + fexp + " -p " + pval + " -h "
    cmd = cmd + x + " -o " + outdir + "/" + x + ".adj" 
    cmd = cmd + " > " + outdir + "/salida-" + x + ".log" 
    os.system(cmd)
    print(cmd)

if __name__ == '__main__':
    fname = sys.argv[1]
    fgenlist = sys.argv[2]
    procs = int(sys.argv[3])
    pvalue = sys.argv[4]
    outdir = sys.argv[5]
    
    print("ParAracne using " + str(procs) + " processors")
    print("fname " + fname + " fgenlist " + fgenlist + " outdir " + outdir)

    df = pd.read_csv(fgenlist, sep='\t')
    genes = df.iloc[:, 0].tolist()

    p = Pool(procs)
    fparam=partial(f, fexp=fname, outd=outdir, pval=pvalue)
    result = p.map(fparam, genes)
    p.close()
    p.join()

    os.mkdir(outdir + "/log") 
    os.system("mv " + outdir + "/*.log "  + outdir + "/log")
