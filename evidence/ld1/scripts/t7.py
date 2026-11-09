from ask import *
import math
for r,n,pv in [(1.0,20,1-2**-20),(1.0,20,1.0),(3.0,10,1-4**-10),(1.0,53,1.0)]:
    rows=[[r,n,pv,k,k,0] for k in range(1,n+1)]
    a=ask("CUMPRINC",rows)["rows"]
    pm=ask("PMT",[[r,n,pv,0,0]])["rows"][0]["bits"]
    print(r,n,pv,"PMT",pm,f(pm))
    for k,x in zip(range(1,n+1),a):
        P=f(x["bits"]); m=n-k+1
        print("  k",k,x["bits"], "P/(pmt*v^m)-1 = %.3e"%(P/(f(pm)*(1+r)**-m)-1), "P*(1+r)^m = %r"%(P*(1+r)**m))
