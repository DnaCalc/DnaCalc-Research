from ask import *
loans=[(0.1,400,100000.0),(0.0371,300,250000.0),(0.21,150,5000.0),(0.007,360,1e6),(1.0,40,12345.67),(0.75,40,3.3),(0.1013,352,79905.9)]
with open("d5.txt","w") as o:
  for r,n,pv in loans:
    rows=[[r,n,pv,k,k,1] for k in range(1,n+1)]
    a=ask("CUMPRINC",rows)["rows"]
    pm=ask("PMT",[[r,n,pv,0,1]])["rows"][0]["bits"]
    pm0=ask("PMT",[[r,n,pv,0,0]])["rows"][0]["bits"]
    a0=ask("CUMPRINC",[[r,n,pv,k,k,0] for k in range(1,n+1)])["rows"]
    pp=ask("PPMT",[[r,k,n,pv,0,1] for k in range(1,n+1)])["rows"]
    for row,x,y,z in zip(rows,a,a0,pp): o.write(" ".join(h(float(v)) for v in row)+" "+pm+" "+x["bits"]+" "+pm0+" "+y["bits"]+" "+z["bits"]+"\n")
