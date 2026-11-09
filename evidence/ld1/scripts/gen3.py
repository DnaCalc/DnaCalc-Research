from ask import *
with open("d3.txt","w") as o:
  for r,n,pv in [(1.0,20,1-2**-20),(1.0,53,1.0),(3.0,10,1-4**-10)]:
    rows=[[r,n,pv,k,k,0] for k in range(1,n+1)]
    a=ask("CUMPRINC",rows)["rows"]
    pm=ask("PMT",[[r,n,pv,0,0]])["rows"][0]["bits"]
    for row,x in zip(rows,a): o.write(" ".join(h(float(v)) for v in row)+" "+pm+" "+x["bits"]+"\n")
