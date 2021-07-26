import sys
import pandas as pd

# Takes ARACNE adj and builds an undirected sif
adjmat = sys.argv[1]
out = sys.argv[2]
lines = [line.rstrip('\n') for line in open(adjmat)]
sources = []
targets = []
mis = []

for l in lines:
  vals = l.split()
  name = vals[0]
  for i in range(1, len(vals), 2):
    st = [name, vals[i]]
    st.sort()
    sources.append(st[0])
    targets.append(st[1])
    mis.append(vals[i+1])

dict = {'source': sources, 'mi':mis, 'target' : targets}
df = pd.DataFrame(dict)
df.sort_values(by = ["mi"], ascending = False, inplace = True)
df.drop_duplicates(inplace = True)
df.to_csv(out, sep="\t", index=False)