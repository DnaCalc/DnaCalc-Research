from ask import *
import sys, json
o=json.load(open("oy1.json"))
out=[]
for idx in range(int(sys.argv[1]), int(sys.argv[2])):
    x=o[idx]; a=x["args"]; PR=a[5]
    def P(ys):
        rows=[a[:5]+[y]+a[6:] for y in ys]
        return [f(z["bits"]) for z in ask("ODDFPRICE",rows,priority="interactive")["rows"]]
    H=1e-7
    def N(y):
        p,b,c=P([y,y+H,y-H]); d=(b-c)/(2*H)
        return y-(p-PR)/d
    yK=f(x["y"])
    y=yK; seq=[y]
    try:
        for k in range(8):
            lo,hi=max(y-0.5,2e-7),y
            if N(lo)>y: break
            for i in range(55):
                m=(lo+hi)/2
                if N(m)<y: lo=m
                else: hi=m
            y=(lo+hi)/2; seq.append(y)
    except Exception as ex:
        seq.append(str(ex))
    out.append(dict(args=a,seq=seq))
    print(idx, a, seq[-1], len(seq), flush=True)
json.dump(out,open("back2-%s.json"%sys.argv[1],"w"))
