import json, struct, math
f=lambda h: struct.unpack(">d",bytes.fromhex(h))[0]
d=json.load(open("oz1.json")); T=24.61111111111111
prs=d["prs"]; ys=[f(y) for y in d["ys"]]
def P(Y): return 100*(1+Y/2)**(-T)
def dP(Y): return -T/2*100*(1+Y/2)**(-T-1)
def N(y,pr): return y-(P(y)-pr)/dP(y)
for pr,Y in list(zip(prs,ys))[4:40:3]:
    root=2*((100/pr)**(1/T)-1)
    y=Y; seq=[y]
    for k in range(10):
        lo,hi=-1.99,y
        if N(lo,pr)>y: break
        for i in range(200):
            m=(lo+hi)/2
            if N(m,pr)<y: lo=m
            else: hi=m
        y=(lo+hi)/2; seq.append(y)
    lin=2*(100/pr-1)/T
    print("pr %.2f root %.12f steps %d earliest %.12f  lin %.12f  ratio %.6f"%(pr,root,len(seq)-1,seq[-1],lin,seq[-1]/root), ["%.6f"%s for s in seq[1:]])
