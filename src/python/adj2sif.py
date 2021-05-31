
# Convierte una matriz de adyacencia de ARACNe en un archivo 
# de interacciones SIF con la matriz triangular 
import sys

adjmat = sys.argv[1]
lines = [line.rstrip('\n') for line in open(adjmat)]

c = 0
for l in lines:
  vals = l.split()
  name = vals[0]
  for i in range(2*c + 1, len(vals), 2):
    print(name + "\t" + vals[i] + "\t" + vals[i+1])
  c = c + 1 #Elemento en la triangular
