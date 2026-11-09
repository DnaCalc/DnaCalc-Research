from ask import *
rows=[]
for n in [352,352.25,352.5,352.999,353]:
    for s,e in [(1,1),(1,10),(50,152),(352,352)]:
        rows.append([0.1013,n,79905.9,s,e,0])
for s,e in [(1.5,1.5),(1.9,10.2),(50.7,152.7)]:
    rows.append([0.1013,352,79905.9,s,e,0])
a=ask("CUMPRINC",rows)["rows"]
b=ask("CUMIPMT",rows)["rows"]
for r,x,y in zip(rows,a,b): print(r[1],r[3],r[4],x.get("bits",x),f(x["bits"]) if "bits" in x else "", y.get("bits",y))
p=ask("PMT",[[0.1013,n,79905.9,0,0] for n in [352,352.5,353]])["rows"]
print(p)
