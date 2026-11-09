from ask import *
import json, math
o=json.load(open("oy1.json"))
rows=[]
for x in o:
    a=x["args"]; y=f(x["y"]); hh=1e-6
    rows.append([a[0],a[1],a[2],a[3],a[4],y+hh,a[6],a[7],a[8]])
    rows.append([a[0],a[1],a[2],a[3],a[4],y-hh,a[6],a[7],a[8]])
p=ask("ODDFPRICE",rows)["rows"]
steps=[]
for i,x in enumerate(o):
    d=(f(p[2*i]["bits"])-f(p[2*i+1]["bits"]))/2e-6
    steps.append((x["res"]/(-d), x["y"], x["args"]))
s=sorted(abs(t[0]) for t in steps)
print("max |y err|", s[-8:])
h=[0]*20
for v in s:
    if v>0: h[min(19,int(-math.log10(v)))]+=1
print(h)
json.dump(steps,open("oy2.json","w"))
