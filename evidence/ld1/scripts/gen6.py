from ask import *
import random
random.seed(4)
with open("d6.txt","w") as o:
  for (r,n,pv) in [(1.0,40,12345.67),(1.0,60,987.65),(3.0,25,777.77),(1.0,50,0.3),(0.1,400,100000.0),(0.0371,300,250000.0),(0.75,40,3.3)]:
    pm=ask("PMT",[[r,n,pv,0,0]])["rows"][0]["bits"]
    ranges=set()
    while len(ranges)<150:
        s=random.randint(1,n); e=random.randint(s,n)
        if e>s: ranges.add((s,e))
    ranges=sorted(ranges)
    got=ask("CUMPRINC",[[r,n,pv,s,e,0] for s,e in ranges])["rows"]
    for (s,e),g in zip(ranges,got): o.write(" ".join(h(float(v)) for v in [r,n,pv,s,e,0])+" "+pm+" "+g["bits"]+"\n")
