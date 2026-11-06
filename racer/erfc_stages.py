"""Invert the upper ERFC product, and pair the normal wrapper with public ERFC.

Original MIT tooling. No rational coefficients are fitted, and no match count
chooses a numerical program. Scalar receipts include every original readback.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from fingerprint import CachedService,load_service,Native,bits,value,sha,write
from exact import BINARY64 as B,Span,binary64,inverse
from x87 import load64


def run(checkout,output,cache,library,count=256,offline=False):
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=300),cache,offline)
    native=Native(library)
    centers=(.5001,.7,1.,1.11567925,1.1217,1.5,2.5,3.,4.,5.5,8.,16.,26.)
    xs=[f'{int(bits(x),16)-count//2+k:016x}' for x in centers for k in range(count)]
    answers=service.ask('ERFC.PRECISE',[[client.number(x)] for x in xs])['rows']
    square={x:f'{native.publish(native.leaf("mul",load64(int(x,16)),load64(int(x,16)),0x133f)):016x}' for x in xs}
    exponents=sorted({bits(-value(z)) for z in square.values()})
    exponential=dict(zip(exponents,service.ask('EXP',[[client.number(x)] for x in exponents])['rows']))
    counts=Counter();records=[];span_hash=hashlib.sha256()
    for x,answer in zip(xs,answers):
        e=exponential[bits(-value(square[x]))]
        if answer['type']!=e['type'] or answer['type']!='number':
            raise ValueError('normal finite upper ERFC/EXP must be numeric')
        target=binary64(answer['bits']);factor=binary64(e['bits'])
        domain=inverse(Span(target,target),factor,'mul',True,'stored_pc64')
        span_hash.update(json.dumps(domain.json(),sort_keys=True).encode()+b'\n')
        coordinates=B.coordinates(domain)
        n=0 if coordinates is None else coordinates[1]-coordinates[0]+1
        if n>4:
            raise ValueError('unexpectedly broad F preimage')
        row=dict(input=x,square=square[x],exponential=e,erfc=answer,binary64_F_count=n,
            binary64_F_values=[] if coordinates is None else [bits(float(B.value(k))) for k in range(coordinates[0],coordinates[1]+1)])
        if not n:
            row['adjacent_F_bits']=[bits(float(B.value(k))) for k in
                (B.floor_coordinate(domain.low),B.ceiling_coordinate(domain.high))]
        counts['rows']+=1;counts['binary64_F_'+str(n)]+=1
        records.append(row)
    # Independent public wrapper calls use their own staged z, not an assumed
    # exact inverse of multiplying an ERFC argument by sqrt(2).
    ys=[f'{int(bits(y),16)-64+k:016x}' for y in (.9,1.,2.,4.,8.,16.,32.,37.5) for k in range(128)]
    c=load64(0x3fe6a09e667f3bcd)
    zs={y:f'{native.publish(native.leaf("mul",load64(int(y,16)),c,0x133f)):016x}' for y in ys}
    distinct_z=sorted(set(zs.values()))
    public_erfc=dict(zip(distinct_z,service.ask('ERFC.PRECISE',[[client.number(z)] for z in distinct_z])['rows']))
    normals=service.ask('NORM.S.DIST',[[client.number(sign*value(y)),client.boolean(True)] for y in ys for sign in (-1.,1.)])['rows']
    wrappers=[]
    for i,y in enumerate(ys):
        q=public_erfc[zs[y]]
        if q['type']!='number':
            raise ValueError('normal wrapper ERFC control must be numeric')
        half=.5*value(q['bits'])
        if half<2.**-1022:
            half=0.
        for j,sign in enumerate((-1.,1.)):
            prediction=bits(half if sign<0 else 1.-half)
            row=dict(input=bits(sign*value(y)),z=zs[y],public_erfc=q,
                     answer=normals[2*i+j],prediction=prediction)
            counts['normal_rows']+=1
            counts['normal_wrapper_matched']+=row['answer']==dict(type='number',bits=prediction)
            wrappers.append(row)
    result=dict(format='dnacalc-erfc-product-preimages/1',complete=True,counts=counts,records=records,
        normal_wrappers=wrappers,receipts=service.receipts,F_domain_sha256=span_hash.hexdigest(),
        source_bindings={str(p):sha(p) for p in (Path(__file__),Path(library))},
        selection='Thirteen fixed upper-ERFC sites, adjacent double windows; eight fixed normal-wrapper windows with both signs.',
        limits='The F inverse conditions on z=storedPC64(x*x), independently measured public EXP(-z), and a final storedPC64 E*F. Empty binary64 preimages refute a stored-F composition under these premises; they do not identify a rational, coefficient, ABI, or the alternative association E*N/D. Normal wrapper predictions use separately observed ERFC values at staged z. No numerical kernel is selected by counts.')
    write(output,result);print(json.dumps(counts,indent=2),flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('checkout');p.add_argument('output')
    p.add_argument('--cache',required=True);p.add_argument('--library',required=True)
    p.add_argument('--count',type=int,default=256);p.add_argument('--offline',action='store_true')
    a=p.parse_args()
    if not 2<=a.count<=512:
        raise ValueError('small windows only')
    run(a.checkout,a.output,a.cache,a.library,a.count,a.offline)
