import pandas as pd
import os, sys

adjfolder = sys.argv[1]
outfile = sys.argv[2]
file_names = [fn for fn in os.listdir(adjfolder) if fn.endswith(".adj")]

os.chdir(adjfolder)
wlines = []
for fname in file_names:
  fh = open(fname, "r")
  lines = fh.read().splitlines()
  last_line = lines[-1]
  wlines.append(last_line)
  fh.close()

with open(outfile, 'w') as f:
    for l in wlines:
        if not l.startswith(">"):
            f.write("%s\n" % l)

