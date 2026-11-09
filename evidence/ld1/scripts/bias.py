import sys,struct
from fractions import Fraction as Q
f=lambda h: struct.unpack(">d",bytes.fromhex(h))[0]
obs={0.25:0.091,0.5:-0.201,0.625:0.288,0.75:-1.128,0.875:-0.616,1.0:0,1.5:0.889,3.0:0,0.1:-0.481,0.3:0.053,2.0:0.304}
for line in sys.stdin:
    p=line.split(); r=float(p[0]); R=Q(r)
    u1=Q(f(p[2])); ur=Q(f(p[4])); em1=Q(f(p[6]))
    e1=float((u1-1/(1+R))/(1/(1+R)))/2**-52
    e2=float((ur-R/(1+R))/(R/(1+R)))/2**-52
    e3=float((-em1-R/(1+R))/(R/(1+R)))/2**-52
    print("%5s obs %+.3f  u1 %+.3f  RN(u1*r) %+.3f  -em1 %+.3f"%(r,obs[r],e1,e2,e3))
