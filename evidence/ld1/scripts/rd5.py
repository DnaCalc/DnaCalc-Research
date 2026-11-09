exec(open("rd4.py").read().split("# restrict to near-tie rows")[0])
from collections import Counter
C=Counter()
def V(x,shift):
    v=Q(x); neg=shift<0
    for p in chunks(abs(shift)):
        v=E64.rn(v*const(p,neg))
    return v
def Vold(x,shift):
    c=E64.rn(Q(10)**abs(shift)); return E64.rn(Q(x)*c) if shift>0 else E64.rn(Q(x)/c)
for x,s,R in rows:
    v=V(x,s); w=v.numerator//v.denominator; fr=v-w
    vo=Vold(x,s); wo=vo.numerator//vo.denominator; fro=vo-wo
    ex=Q(x)*Q(10)**s; we=ex.numerator//ex.denominator
    if fr==Q(1,2) or fro==Q(1,2) or abs(ex-we-Q(1,2))<Q(2)**(we.bit_length()-62):
        up = (R==we+1)
        C[("S4tie" if fr==Q(1,2) else "S4below" if fr<Q(1,2) else "S4above", "oldtie" if fro==Q(1,2) else "oldbelow" if fro<Q(1,2) else "oldabove", "exact>half" if ex-we>Q(1,2) else "exact<half" if ex-we<Q(1,2) else "exact=half", "excel up" if up else "excel down", s>0)]+=1
for k,v in sorted(C.items(), key=lambda t:-t[1]): print(v,k)
