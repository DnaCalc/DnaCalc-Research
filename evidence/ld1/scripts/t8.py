from ask import *
import struct
def nxt(x,i): return struct.unpack(">d",struct.pack(">Q",struct.unpack(">Q",struct.pack(">d",x))[0]+i))[0]
r,n=0.75,40
pvs=[nxt(3.3,i) for i in range(3000)]
pm=[x["bits"] for x in ask("PMT",[[r,n,pv,0,0] for pv in pvs])["rows"]]
for k in [40,20,1]:
    P=[x["bits"] for x in ask("CUMPRINC",[[r,n,pv,k,k,0] for pv in pvs])["rows"]]
    d={}
    conflict=0
    for a,b in zip(pm,P):
        if a in d and d[a]!=b: conflict+=1
        d.setdefault(a,b)
    print("k",k,"distinct pmt",len(set(pm)),"distinct P",len(set(P)),"conflicts",conflict)
    open("pvscan-%d.txt"%k,"w").write("\n".join("%s %s %s"%(h(pv),a,b) for pv,a,b in zip(pvs,pm,P)))
