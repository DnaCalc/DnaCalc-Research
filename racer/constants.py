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


def normal_constants(source_dir):
    """Published literals only; never inspect an Excel output."""
    import re, hashlib
    from decimal import Decimal, localcontext
    source_dir = Path(source_dir)
    source = (source_dir/'cdflib.c').read_text()
    source = source[source.index('struct TupleDD cumnor(double x)'):]
    source = source[:source.index('struct TupleDD cumpoi')]
    tables = {}
    for name, count, body in re.findall(r'double ([abcdpq])\[(\d+)\] = \{([^}]+)\}', source):
        values = [v.strip() for v in body.split(',')]
        assert len(values) == int(count)
        tables['direct_'+name] = values
    assert len(tables) == 6
    source = (source_dir/'calerf.f').read_text()
    for name in 'ABCDPQ':
        match = re.search(r'CD    DATA '+name+r'/([^/]+)/', source)
        assert match
        body = re.sub(r'\nCD   [ 0-9]', '', match[1])
        tables['erfc_'+name.lower()] = [v.strip().replace('D','e') for v in body.split(',')]
    tables['direct_sqrpi'] = ['3.9894228040143267794e-1']
    tables['erfc_sqrpi'] = ['5.6418958354775628695e-1']
    with localcontext() as c:
        c.prec = 100
        root = Decimal(2).sqrt()
        tables['sqrt2'] = [str(root)]
        tables['invsqrt2'] = [str(1/root)]
    def encode(literal):
        rational = Q(literal)
        a, b = B64.rn(rational), E64.rn(rational)
        e = log2floor(b)
        raw = int(b/pow2(e-63)).to_bytes(8, 'little')+(e+16383).to_bytes(2, 'little')
        return dict(decimal=literal, binary64=struct.pack('>d',float(a)).hex(), ext80_le=raw.hex())
    return dict(provenance='Published Cody1969/SPECFUN CALERF and SciPy1.14.1 CDFLIB cumnor tables; round each stated decimal once. sqrt2 and reciprocal derived at100decimal digits, never tuned.',
        sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [source_dir/'cdflib.c',source_dir/'calerf.f']},
        tables={k:[encode(v) for v in values] for k,values in tables.items()})


def hart_constants(source_dir):
    """Fixed numerical table in West2004 Figure2; no program text is ported."""
    out = normal_constants(source_dir)
    tables = dict(
        hart_p=['220.206867912376','221.213596169931','112.079291497871',
                '33.912866078383','6.37396220353165','0.700383064443688','3.52624965998911e-2'],
        hart_q=['440.413735824752','793.826512519948','637.333633378831',
                '296.564248779674','86.7807322029461','16.064177579207','1.75566716318264','8.83883476483184e-2'])
    def encode(literal):
        a,b = B64.rn(Q(literal)),E64.rn(Q(literal))
        e = log2floor(b)
        raw = int(b/pow2(e-63)).to_bytes(8,'little')+(e+16383).to_bytes(2,'little')
        return dict(decimal=literal,binary64=struct.pack('>d',float(a)).hex(),ext80_le=raw.hex())
    out['tables'].update({k:[encode(v) for v in values] for k,values in tables.items()})
    out['hart_provenance'] = dict(source='Graeme West, Better Approximations to Cumulative Normal Functions, December8,2004,Figure2,p3',
        url='https://s2.smu.edu/~aleskovs/emis/sqc2/accuratecumnorm.pdf',
        scope='Only the numerical coefficient table and published body bound7.07106781186547. No VBA program text copied. Shorter published decimal table, not the uninspected full Hart1968 table.',
        body_bound='7.07106781186547')
    return out


if __name__=='__main__':
    import sys
    out=hart_constants(sys.argv[2]) if len(sys.argv)>3 and sys.argv[3]=='hart' else normal_constants(sys.argv[2]) if len(sys.argv)>2 else generate()
    with Path(sys.argv[1]).open('x') as f:json.dump(out,f,indent=2)
    print(json.dumps({k:{p:v for p,v in c.items()if p not in ['lower','upper']}for k,c in out.get('constants',{}).items()}))
