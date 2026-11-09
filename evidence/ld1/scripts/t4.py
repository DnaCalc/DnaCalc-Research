from ask import *
import sys, random, math
sys.path.insert(0,"/home/"+__import__("getpass").getuser()+"/projects/DnaCalc/DnaCalc-Research/racer")
from intervals import B64,E64
from fractions import Fraction as Q
random.seed(2)
loans=[(0.05/12,360,200000.0),(0.1,30,1000.0),(0.013,100,5e5),(0.2,50,777.0)]
for (r,n,pv) in loans:
    R=Q(r);v=1/(1+R);vn=v**n
    ex=lambda s,e: sum(-Q(pv)*R*v**(n-k+1)/(1-vn) for k in range(s,e+1))
    ranges=[(1,n),(1,n//2),(n//2,n),(2,n-1),(1,2),(n-1,n),(5,9),(n//3,2*n//3)]
    got=[f(x["bits"]) for x in ask("CUMPRINC",[[r,n,pv,s,e,0] for s,e in ranges])["rows"]]
    print(r,n,pv)
    for (s,e),g in zip(ranges,got):
        X=ex(s,e)
        # balance scale
        print("  ",s,e,"err ulp %.2f"%(float(Q(g)-X)/math.ulp(g)), "  pv/|res| %.1f"%(pv/abs(g)))
