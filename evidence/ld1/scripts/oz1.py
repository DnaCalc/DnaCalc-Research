from ask import *
import math
A=[39763,44256,39736,39873,0.0,None,100,2,0]
p=[f(x["bits"]) for x in ask("ODDFPRICE",[A[:5]+[y]+A[6:] for y in [0.05,0.1,0.2]])["rows"]]
T=[-math.log(pp/100)/math.log(1+y/2) for pp,y in zip(p,[0.05,0.1,0.2])]
print("prices",p,"T",T)
prs=[60+i*0.37 for i in range(100)]
ys=ask("ODDFYIELD",[A[:5]+[pr]+A[6:] for pr in prs])["rows"]
import json
json.dump(dict(T=T,prs=prs,ys=[y.get("bits") for y in ys]),open("oz1.json","w"))
Tm=T[0]
for pr,y in list(zip(prs,ys))[:10]:
    Y=f(y["bits"]); root=2*((100/pr)**(1/Tm)-1)
    print(pr, Y, root, (Y-root))
