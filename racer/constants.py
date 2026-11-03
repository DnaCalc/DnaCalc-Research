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


def bessel_asymptotic_constants():
    """Fixed decimal quantization of universal constants, never of Excel values."""
    a,b=atan(Q(1,5),128),atan(Q(1,239),64)
    pi=(16*a[0]-4*b[1],16*a[1]-4*b[0])
    def quantize(span,digits,mode):
        scale=10**digits
        def integer(q):
            q*=scale;k=q.numerator//q.denominator
            if mode=='ceil':return k+(q!=k)
            if mode=='nearest_even':return k+(q-k>Q(1,2)or q-k==Q(1,2)and k%2==1)
            return k
        lo,hi=map(integer,span)
        if lo!=hi:raise ValueError('ambiguous universal decimal quantization')
        return published_literal(str(Q(lo,scale)))
    phases=[dict(id=0,rule='correct_binary64',values=[published_literal(str(B64.rn(Q(2*n+1,4)*pi[0])))for n in (0,1)])]
    scales=[dict(id=0,rule='correct_binary64',value=published_literal(str(B64.rn(2/pi[1]))))]
    for digits in (8,9,10):
        for mode in ('floor','nearest_even','ceil'):
            phases.append(dict(id=len(phases),rule=mode,decimal_places=digits,
                values=[quantize((Q(2*n+1,4)*pi[0],Q(2*n+1,4)*pi[1]),digits,mode)for n in (0,1)]))
            scales.append(dict(id=len(scales),rule=mode,decimal_places=digits,value=quantize((2/pi[1],2/pi[0]),digits,mode)))
    return dict(provenance='DLMF10.17.1--4; Machin rational pi bounds. Universal constants only: phase=(2n+1)pi/4, amplitude squared=2/pi. '
                'Fixed floor/nearest-even/ceil decimal quantizers at8,9,10places and correctly rounded binary64; no observed answer enters generation.',
                phases=phases,scales=scales)


def published_literal(literal):
    """Round an explicitly stated rational/decimal once; never inspect targets."""
    a,b=B64.rn(Q(literal)),E64.rn(Q(literal));magnitude=abs(b)
    if not magnitude:raw=bytes(10)
    else:
        e=log2floor(magnitude);significand=int(magnitude/pow2(e-63))
        raw=significand.to_bytes(8,'little')+(e+16383+(0x8000 if b<0 else 0)).to_bytes(2,'little')
    return dict(decimal=literal,binary64=struct.pack('>d',float(a)).hex(),ext80_le=raw.hex())


def as_j01_constants(page369,page370):
    """A&S9.4.1/9.4.4 facts, independently checked against the scanned pages."""
    import hashlib
    values=dict(AS_J0=['1','-2.2499997','1.2656208','-0.3163866','0.0444479','-0.0039444','0.0002100'],
                AS_J1=['0.5','-0.56249985','0.21093573','-0.03954289','0.00443319','-0.00031761','0.00001109'])
    return dict(provenance='A&S Handbook pages369/370 equations9.4.1 and9.4.4; coefficient facts only, no approximation remainder or target-derived value.',
                source_sha256={Path(p).name:hashlib.sha256(Path(p).read_bytes()).hexdigest()for p in (page369,page370)},
                tables={k:[published_literal(x)for x in xs]for k,xs in values.items()},
                constants={'ONE_NINTH':published_literal('1/9')})


def cody_j0_constants(source_path):
    """Cody SPECFUN's stated decimal coefficient facts, never Excel-derived."""
    import hashlib,re
    source_path=Path(source_path);source=source_path.read_text()
    code=''.join(line[6:]for line in source.splitlines()if line.startswith('CD'))
    def table(name):
        match=re.search(r'\bDATA\s+'+name+r'/([^/]+)/',code)
        if match is None:raise ValueError('missing published table '+name)
        return [x.strip().replace('D','e')for x in match[1].split(',')]
    def constant(name):
        match=re.search(r'\b'+name+r'/([^/]+)/',code)
        if match is None:raise ValueError('missing published root literal '+name)
        return match[1].strip().replace('D','e')
    tables={name:table(name)for name in ('PJ0','QJ0')}
    if {k:len(v)for k,v in tables.items()}!={'PJ0':7,'QJ0':5}:raise ValueError('wrong published table length')
    return dict(provenance='Cody SPECFUN CALJY0 small-argument rational coefficient literals and split-root facts only; no program text copied and no Excel-derived coefficient.',
                source_url='https://www.netlib.org/specfun/j0y0',source_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
                tables={name:[published_literal(x)for x in values]for name,values in tables.items()},
                constants={name:published_literal(constant(name))for name in ('XJ0','XJ01','XJ02')})


if __name__=='__main__':
    import sys
    if len(sys.argv)>2 and sys.argv[2]=='as-j01':out=as_j01_constants(*sys.argv[3:])
    else:out=cody_j0_constants(sys.argv[2]) if len(sys.argv)>3 and sys.argv[3]=='cody-j0' else hart_constants(sys.argv[2]) if len(sys.argv)>3 and sys.argv[3]=='hart' else normal_constants(sys.argv[2]) if len(sys.argv)>2 else generate()
    with Path(sys.argv[1]).open('x') as f:json.dump(out,f,indent=2)
    print(json.dumps({k:{p:v for p,v in c.items()if p not in ['lower','upper']}for k,c in out.get('constants',{}).items()}))
