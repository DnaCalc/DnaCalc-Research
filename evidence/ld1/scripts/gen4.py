from ask import *
import random
random.seed(5)
with open("d4.txt","w") as o:
  for r,n,pv in [(1.0,20,1.0),(1.0,40,12345.67),(3.0,25,777.77),(1.0,30,0.1),(0.5,60,1000.0),(0.75,40,3.3)]:
    rows=[[r,n,pv,k,k,0] for k in range(1,n+1)]
    a=ask("CUMPRINC",rows)["rows"]
    pm=ask("PMT",[[r,n,pv,0,0]])["rows"][0]["bits"]
    for row,x in zip(rows,a): o.write(" ".join(h(float(v)) for v in row)+" "+pm+" "+x["bits"]+"\n")
