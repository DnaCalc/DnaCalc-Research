from ask import *
import random, json
random.seed(7)
cases=[]
while len(cases)<400:
    iss=random.randint(36000,44000)
    freq=random.choice([1,2,4]); basis=random.choice([0,1,2,3,4])
    # first coupon on a schedule from maturity: pick maturity then step back
    yrs=random.randint(2,30)
    mat=iss+random.randint(400,365*yrs)
    rate=round(random.uniform(0.0,0.12),4)
    pr=round(random.uniform(60,130),2)
    red=random.choice([100,100,100,105,95])
    cases.append([iss,mat,freq,basis,rate,pr,red])
# derive first coupon via COUPNCD of issue date (schedule from maturity)
fc=ask("COUPNCD",[[c[0],c[1],c[2],c[3]] for c in cases])["rows"]
rows=[]
for c,x in zip(cases,fc):
    if x["type"]!="number": continue
    first=f(x["bits"])
    if first>=c[1]: continue
    sett=random.randint(c[0]+1,int(first)-1) if first-c[0]>2 else None
    if sett is None: continue
    rows.append([sett,c[1],c[0],first,c[4],c[5],c[6],c[2],c[3]])
y=ask("ODDFYIELD",rows)["rows"]
ok=[(r,x) for r,x in zip(rows,y) if x["type"]=="number"]
pr=ask("ODDFPRICE",[[r[0],r[1],r[2],r[3],r[4],f(x["bits"]),r[6],r[7],r[8]] for r,x in ok])["rows"]
out=[]
for (r,x),p in zip(ok,pr):
    if p["type"]!="number": continue
    out.append(dict(args=r,y=x["bits"],p=p["bits"],res=f(p["bits"])-r[5]))
json.dump(out,open("oy1.json","w"))
import statistics
res=[abs(o["res"]) for o in out]
print(len(rows),len(ok),len(out))
res.sort(); print(res[-10:], sum(1 for x in res if x==0), sum(1 for x in res if x<1e-12))
