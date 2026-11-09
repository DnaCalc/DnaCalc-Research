from ask import *
import sys
loans=[(0.1,400,100000.0),(0.0371,300,250000.0),(0.21,150,5000.0),(0.007,360,1e6)]
with open("d2.txt","w") as o:
  for r,n,pv in loans:
    rows=[[r,n,pv,k,k,0] for k in range(1,n+1)]
    a=ask("CUMPRINC",rows)["rows"]
    pm=ask("PMT",[[r,n,pv,0,0]])["rows"][0]["bits"]
    for row,x in zip(rows,a):
        o.write(" ".join(h(float(v)) for v in row)+" "+pm+" "+x["bits"]+"\n")
