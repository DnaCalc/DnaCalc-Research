import sys, struct, itertools, random
sys.path.insert(0,"/home/koderbot/projects/DnaCalc/DnaCalc-Research/racer")
from fractions import Fraction as Q
from intervals import B64, E64
def f(hx): return struct.unpack(">d", bytes.fromhex(hx))[0]
rows=[l.split() for l in open(sys.argv[1])]
random.seed(1); rows=random.sample(rows,120)
consts=dict(x.split("=") for x in sys.argv[2:])
C={k:Q(f(v)) for k,v in consts.items()}
data=[(Q(f(r[0])),Q(f(r[2]))) for r in rows]
def R(x,p): return B64.rn(x) if p==53 else E64.rn(x)
def ap(op,a,b):
    if op=="*": return a*b
    if op=="/": return a/b
    if op=="\\": return b/a
names=list(C)
ops=["*","/","\\"]
results=[]
for depth in (1,2,3):
    for cs in itertools.product(names,repeat=depth):
        for os_ in itertools.product(ops,repeat=depth):
            for precs in itertools.product((53,64),repeat=depth-1):
                ok=0
                for pv,P in data:
                    x=pv
                    for i,(c,o) in enumerate(zip(cs,os_)):
                        x=ap(o,x,C[c])
                        x=R(x,precs[i]) if i<depth-1 else B64.rn(x)
                    if x==P: ok+=1
                if ok>=len(data)*0.9: results.append((ok,cs,os_,precs))
results.sort(reverse=True)
for r in results[:15]: print(r)
print("total>=90%",len(results))
