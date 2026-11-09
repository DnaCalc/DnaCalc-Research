from ask import *
import struct
def nxt(x,i): return struct.unpack(">d",struct.pack(">Q",struct.unpack(">Q",struct.pack(">d",x))[0]+i))[0]
pvs=[nxt(3.3,i) for i in range(300)]
a=[x["bits"] for x in ask("CUMPRINC",[[0.75,1000,pv,1000,1000,0] for pv in pvs])["rows"]]
b=[x["bits"] for x in ask("CUMPRINC",[[0.75,1000,pv,1000,1000,0] for pv in reversed(pvs)],order="serial")["rows"]][::-1]
c=[ask("CUMPRINC",[[0.75,1000,pv,1000,1000,0]])["rows"][0]["bits"] for pv in pvs[:20]]
print(sum(x==y for x,y in zip(a,b)), len(a), sum(x==y for x,y in zip(a[:20],c)))
old=[l.split()[2] for l in open("pvscan-big1.txt")][:300]
print(sum(x==y for x,y in zip(a,old)))
