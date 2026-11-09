from ask import *
import struct,sys
def nxt(x,i): return struct.unpack(">d",struct.pack(">Q",struct.unpack(">Q",struct.pack(">d",x))[0]+i))[0]
r=float(sys.argv[1]); n=int(sys.argv[2]); k=int(sys.argv[3]); tag=sys.argv[4]
pvs=[nxt(3.3,i) for i in range(2000)]
pm=[x["bits"] for x in ask("PMT",[[r,n,pv,0,0] for pv in pvs])["rows"]]
P=[x["bits"] for x in ask("CUMPRINC",[[r,n,pv,k,k,0] for pv in pvs])["rows"]]
open("pvscan-%s.txt"%tag,"w").write("\n".join("%s %s %s"%(h(pv),a,b) for pv,a,b in zip(pvs,pm,P)))
d={};c=0
for a,b in zip(pm,P):
    if a in d and d[a]!=b: c+=1
    d.setdefault(a,b)
print("distinct pmt",len(set(pm)),"P",len(set(P)),"conflicts vs pmt",c)
