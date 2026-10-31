"""Certified binary64/Ext80 loads of universal Poisson normalization constants.

No Excel output is an input. Machin's identity and the alternating arctangent
series bracket pi rationally. Decimal's correctly rounded ln/sqrt, each widened
by a full decimal ULP, bracket C=ln(2*pi)/2, S=1-C and G=sqrt(2*pi).
Both endpoints must force exactly the same nearest/even 53- and 64-bit load.
"""
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as Q
from pathlib import Path
import json
import struct
from intervals import B64, E64, log2floor, pow2


def atan(q, n):
    s=Q(0)
    for i in range(n):
        s+=(-1)**i*q**(2*i+1)/(2*i+1)
    t=(-1)**n*q**(2*n+1)/(2*n+1)
    return min(s,s+t),max(s,s+t)


def decimal_bound(q, lower):
    with localcontext() as c:
        c.prec=200
        c.rounding=ROUND_FLOOR if lower else ROUND_CEILING
        return Decimal(q.numerator)/Decimal(q.denominator)


def function_bounds(low, high, name):
    with localcontext() as c:
        c.prec=200
        a=getattr(decimal_bound(low,True),name)()
        b=getattr(decimal_bound(high,False),name)()
        return Q(a)-Q(Decimal(10)**(a.adjusted()-199)),Q(b)+Q(Decimal(10)**(b.adjusted()-199))


def pair(q):return [str(q.numerator),str(q.denominator)]


def generate():
    a,b=atan(Q(1,5),256),atan(Q(1,239),128)
    pi=(16*a[0]-4*b[1],16*a[1]-4*b[0])
    lc=function_bounds(2*pi[0],2*pi[1],'ln')
    cb=(lc[0]/2,lc[1]/2)
    bounds=dict(C=cb,S=(1-cb[1],1-cb[0]),G=function_bounds(2*pi[0],2*pi[1],'sqrt'),TAU=(2*pi[0],2*pi[1]))
    out=dict(provenance='Tier 1 NIST PMF, DLMF 5.11 normalization, Machin identity and own exact interval derivation; Python Decimal correctly rounded elementary functions. Never derived from Excel answers.',
             definitions=dict(C='ln(2*pi)/2',S='1-C: exact Stirling remainder at integer count 1',G='sqrt(2*pi)',TAU='2*pi'),constants={})
    for name,(low,high) in bounds.items():
        assert B64.rn(low)==B64.rn(high)
        assert E64.rn(low)==E64.rn(high)
        b64=B64.rn(low)
        e64=E64.rn(low)
        exponent=log2floor(e64)
        mantissa=int(e64/pow2(exponent-63))
        raw=mantissa.to_bytes(8,'little')+(exponent+16383).to_bytes(2,'little')
        out['constants'][name]=dict(lower=pair(low),upper=pair(high),
            binary64=f'{struct.unpack(">Q",struct.pack(">d",float(b64)))[0]:016x}',
            ext80_le=raw.hex(),rounded_decimal=str(float(b64)))
    return out


if __name__=='__main__':
    import sys
    out=generate()
    with Path(sys.argv[1]).open('x') as f:json.dump(out,f,indent=2)
    print(json.dumps({k:{p:v for p,v in c.items()if p not in ['lower','upper']}for k,c in out['constants'].items()}))
