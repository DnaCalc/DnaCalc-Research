from ask import *
import struct
def nxt(x,i): return struct.unpack(">d",struct.pack(">Q",struct.unpack(">Q",struct.pack(">d",x))[0]+i))[0]
pvs=[nxt(3.3,i) for i in range(200)]
for (r,n,k) in [(0.75,1000,1000),(0.75,40,20),(0.1,400,300)]:
    a=[f(x["bits"]) for x in ask("CUMPRINC",[[r,n,pv,k,k,0] for pv in pvs])["rows"]]
    b=[f(x["bits"]) for x in ask("CUMPRINC",[[r,n,pv*2,k,k,0] for pv in pvs])["rows"]]
    c=[f(x["bits"]) for x in ask("CUMPRINC",[[r,n,pv*1024,k,k,0] for pv in pvs])["rows"]]
    print(r,n,k, sum(2*x==y for x,y in zip(a,b)), sum(1024*x==y for x,y in zip(a,c)), len(a))
