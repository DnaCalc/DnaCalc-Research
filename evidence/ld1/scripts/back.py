from ask import *
import sys
A=[39763,44256,39736,39873,0.0575,None,100,2,0]; PR=84.5
def P(ys):
    rows=[A[:5]+[y]+A[6:] for y in ys]
    return [f(x["bits"]) for x in ask("ODDFPRICE",rows)["rows"]]
H=1e-7
def N(y):
    a,b,c=P([y,y+H,y-H]); d=(b-c)/(2*H)
    return y-(a-PR)/d
yK=f("3fb3c65d22ec15b6")
# root by bisection
lo,hi=0.0,0.5
for i in range(80):
    m=(lo+hi)/2
    if P([m])[0]>PR: lo=m
    else: hi=m
root=lo; print("root",repr(root),"errK",yK-root)
y=yK
for k in range(6):
    # find z<y with N(z)=y
    lo,hi=y-0.5,y
    if lo<2e-7: lo=2e-7
    if N(lo)>y: print("no predecessor >=0", N(lo)); break
    for i in range(60):
        m=(lo+hi)/2
        if N(m)<y: lo=m   # N increasing? from left, N(z) increasing toward root
        else: hi=m
    y=(lo+hi)/2
    print(k,repr(y), "err", y-root)
