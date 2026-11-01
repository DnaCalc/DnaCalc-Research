"""Finite GEOMEAN x87 programs using the existing native ISA interpreter.

The four paths use a product log with division, a reciprocal multiplier,
sum(log2(x)), or sum(ln(x)). All constants are 0, 1 or the exact list length,
plus the x87 ROM logarithm constants. No fitted leaf or selector is admitted.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
import itertools
import json
from pathlib import Path
import platform
import struct
import time

from x87_programs import Builder,Native,WORDS,digest
from bessel import require


def flatten(value):
    return list(itertools.chain.from_iterable(map(flatten,value)))if isinstance(value,list)else[value]


def fold(items,order,combine,identity):
    if order==2:
        if len(items)==0:return identity
        if len(items)==1:return items[0]
        middle=len(items)//2
        return combine(fold(items[:middle],2,combine,identity),fold(items[middle:],2,combine,identity))
    result=identity
    for item in items if order==0 else items[::-1]:result=combine(result,item)
    return result


def coordinates(row):
    if row['function']!='GEOMEAN'or row['result_index']is not None:return None,'surface'
    inputs=flatten(row['args'])
    if not all(isinstance(x,str)and x.startswith('0x')and len(x)==18 for x in inputs):return None,'non_numeric_argument'
    values=[Q(struct.unpack('>d',bytes.fromhex(x[2:]))[0])for x in inputs]
    if not(2<=len(values)<=32)or any(x<=0 for x in values):return None,'length_or_domain'
    # A 22-binade margin bounds every permitted rounded partial product away
    # from binary64 under/overflow, including all three declared reduction trees.
    for order in range(3):
        def multiply(a,b):
            value=a*b
            if not Q(1,2**1000)<=value<=Q(2**1000):raise ValueError('partial product magnitude')
            return value
        try:fold(values,order,multiply,Q(1))
        except ValueError:return None,'partial_product_bound'
    if not row['expected'].startswith('0x'):return None,'non_numeric_expected'
    return [int(x,16)for x in inputs],None


@lru_cache(maxsize=128)
def graph(count,path,order):
    g=Builder(count);one=g.one;zero=g.constant('0000000000000000')
    n=g.constant(struct.pack('>d',float(count)).hex())
    def op(kind,a=0,b=0,site=None):
        result=g.op(kind,a,b,store=False)
        if site is not None:g.nodes[result]['site']=site
        return result
    if path<2:
        reduction=fold(list(range(count)),order,lambda a,b:op('mul',a,b,0),one)
        g.nodes[reduction]['site']=1
        if path==0:
            logarithm=op('fyl2x',one,reduction,2)
            t=op('div',logarithm,n,3)
        else:
            coefficient=op('div',one,n,3)
            t=op('fyl2x',coefficient,reduction,2)
    else:
        coefficient=one if path==2 else op('fldln2')
        logs=[op('fyl2x',coefficient,i,2)for i in range(count)]
        reduction=fold(logs,order,lambda a,b:op('add',a,b,0),zero)
        g.nodes[reduction]['site']=1
        t=op('div',reduction,n,3)
        if path==3:t=op('mul',t,op('fldl2e'),3)
    g.phase=1
    integer=op('frndint',t)
    fraction=op('sub',t,integer,4)
    small=op('f2xm1',fraction,0,5)
    mantissa=op('add',one,small,6)
    root=op('fscale',mantissa,integer,7)
    g.sites=8
    result=g.graph(root,root)
    require(len(result['nodes'])<=128, 'native graph exceeds interpreter bound')
    return result


def physical_mask(mask):
    return mask|2 if mask&1 else mask


def profiles():
    for identity,(path,order,cw,inner,mask)in enumerate(itertools.product(range(4),range(3),WORDS,WORDS,range(256))):
        yield dict(id=identity,path=path,order=order,cw=cw,inner_cw=inner,mask=mask)


def run(input_path,library,destination):
    start=time.monotonic();dest=Path(destination);dest.mkdir()
    original=[json.loads(line)for line in Path(input_path).read_text().splitlines()]
    rows=[];scope=Counter()
    with(dest/'scoped.jsonl').open('x')as good,(dest/'excluded.jsonl').open('x')as bad:
        for index,row in enumerate(original):
            inputs,reason=coordinates(row)
            record=dict(trace_index=index,raw=row)
            if reason:record['reason']=reason;bad.write(json.dumps(record)+'\n');scope[reason]+=1
            else:record['inputs']=inputs;good.write(json.dumps(record)+'\n');rows.append(record);scope['scoped']+=1
    ps=list(profiles());require(len(ps)==12288 and rows, 'incomplete grammar or empty GEOMEAN scope')
    with(dest/'profiles.jsonl').open('x')as out:
        for profile in ps:out.write(json.dumps(profile)+'\n')
    native=Native(library);alive=list(range(len(ps)));executions=0
    source=digest(Path(__file__))
    with(dest/'first-kills.jsonl').open('x')as kills,(dest/'witnesses.jsonl').open('x')as witnesses:
        for order,row in enumerate(rows):
            before=len(alive);survivors=[]
            for identity in alive:
                p=ps[identity];g=graph(len(row['inputs']),p['path'],p['order'])
                result=native.evaluate(g,row['inputs'],p['cw'],physical_mask(p['mask']),inner_cw=p['inner_cw'])
                got=f'0x{native.publish(result[:10]):016x}';executions+=1
                if got==row['raw']['expected']:survivors.append(identity)
                else:
                    root,raw,values=native.evaluate(g,row['inputs'],p['cw'],physical_mask(p['mask']),inner_cw=p['inner_cw'],trace=True)
                    require(root==result, 'traced native replay differs')
                    kills.write(json.dumps(dict(candidate_id=identity,order=order,trace_index=row['trace_index'],
                        root=root[:10].hex(),raw=raw.hex(),values=values.hex(),got=got,expected=row['raw']['expected']))+'\n')
            alive=survivors
            witnesses.write(json.dumps(dict(order=order,trace_index=row['trace_index'],killed=before-len(alive),remaining=len(alive)))+'\n')
    with(dest/'survivors.json').open('x')as out:json.dump(alive,out);out.write('\n')
    require(digest(Path(__file__))==source, 'source changed during run')  # ubs:ignore — public integrity hashes, not secrets
    cpu=next((line.split(':',1)[1].strip()for line in Path('/proc/cpuinfo').read_text().splitlines()if line.startswith('model name')),'unknown')
    summary=dict(complete=True,programs=len(ps),recorded=len(original),scoped=len(rows),scope=scope,
                 executions=executions,survivors=len(alive),seconds=time.monotonic()-start,
                 cpu=cpu,platform=platform.platform(),native_library_sha256=digest(library),source_sha256=source,
                 input_sha256={str(input_path):digest(input_path)},
                 claim='Finite native x87 GEOMEAN instruction paths, coherent reduction stores and independent PC53/PC64 preparation/reconstruction. '
                       'All logical labels retained, including dead/redundant store bits.',
                 limits='Original discovery inputs only. Native microcode is measured on this CPU, not assumed correctly rounded or equivalent to Windows hardware. '
                        'No private Excel program, new Excel reading, fitted constant or selector is inferred.')
    with(dest/'summary.json').open('x')as out:json.dump(summary,out,indent=2);out.write('\n')
    print(json.dumps(summary))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('input');parser.add_argument('library');parser.add_argument('destination')
    args=parser.parse_args();run(args.input,args.library,args.destination)
