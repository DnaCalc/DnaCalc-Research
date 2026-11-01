"""Independent ISA grammar and rational arithmetic checks for GEOMEAN.

Only the earlier independent checker primitives are reused. No search module
or product source is imported. FYL2X/F2XM1/ROM leaves use a separate native
single-instruction replay; they are measurements, not mathematical substitutes.
"""
from collections import Counter
from fractions import Fraction as Q
import itertools
import json
from pathlib import Path
import sys
import time

from exact import binary64,two
from check import records,require,digest
from x87 import Exact,Instructions,published_bits

WORDS=(0x123f,0x133f)


def coordinates(row):
    if row['function']!='GEOMEAN'or row['result_index']is not None:return None,'surface'
    def flatten(v):
        if not isinstance(v,list):return[v]
        return list(itertools.chain.from_iterable(flatten(a)for a in v))
    args=flatten(row['args'])
    if not all(isinstance(x,str)and x.startswith('0x')and len(x)==18 for x in args):return None,'non_numeric_argument'
    values=list(map(binary64,args))
    if not(2<=len(values)<=32)or any(x<=0 for x in values):return None,'length_or_domain'
    def product(xs,order):
        def multiply(a,b):
            z=a*b
            require(two(-1000)<=z<=two(1000),'partial_product_bound')
            return z
        if order==2:
            if not xs:return Q(1)
            if len(xs)==1:return xs[0]
            middle=len(xs)//2
            return multiply(product(xs[:middle],2),product(xs[middle:],2))
        result=Q(1)
        for x in xs if order==0 else xs[::-1]:result=multiply(result,x)
        return result
    try:
        for order in range(3):product(values,order)
    except ValueError:return None,'partial_product_bound'
    if not row['expected'].startswith('0x'):return None,'non_numeric_expected'
    return [int(x,16)for x in args],None


def grammar(count,path,order):
    nodes=[dict(op='input',a=i)for i in range(count)]
    def const(bits):nodes.append(dict(op='const',bits=bits));return len(nodes)-1
    one=const('3ff0000000000000');zero=const('0000000000000000')
    # A count <=32 is an exact binary64 integer; construct its encoding by
    # integer arithmetic rather than the search's host-float conversion.
    exponent=count.bit_length()-1
    n=const(f'{((exponent+1023)<<52)|((count-(1<<exponent))<<(52-exponent)):016x}')
    phase=0
    def append(kind,a=0,b=0,site=None):
        node=dict(op=kind,a=a,b=b)
        if phase:node['phase']=phase
        if site is not None:node['site']=site
        nodes.append(node);return len(nodes)-1
    def reduce(items,kind,identity):
        if order==2:
            if not items:return identity
            if len(items)==1:return items[0]
            middle=len(items)//2
            return append(kind,reduce(items[:middle],kind,identity),reduce(items[middle:],kind,identity),0)
        result=identity
        for item in items if order==0 else items[::-1]:result=append(kind,result,item,0)
        return result
    if path in (0,1):
        reduced=reduce(list(range(count)),'mul',one);nodes[reduced]['site']=1
        if path==0:
            log=append('fyl2x',one,reduced,2);argument=append('div',log,n,3)
        else:
            coefficient=append('div',one,n,3);argument=append('fyl2x',coefficient,reduced,2)
    else:
        coefficient=one if path==2 else append('fldln2')
        logs=[append('fyl2x',coefficient,i,2)for i in range(count)]
        reduced=reduce(logs,'add',zero);nodes[reduced]['site']=1
        argument=append('div',reduced,n,3)
        if path==3:argument=append('mul',argument,append('fldl2e'),3)
    phase=1
    integer=append('frndint',argument)
    fraction=append('sub',argument,integer,4)
    small=append('f2xm1',fraction,0,5)
    mantissa=append('add',one,small,6)
    result=append('fscale',mantissa,integer,7)
    require(len(nodes)<=128,'native interpreter bound')
    return dict(nodes=nodes,outputs=[result,result],sites=8)


def check(input_path,directory,library,output):
    start=time.monotonic();source=Path(__file__);directory=Path(directory)
    hashes={str(p):digest(p)for p in (source,source.with_name('x87.py'),source.with_name('exact.py'),source.with_name('x87_isa.c'))}
    original=list(records(input_path));scoped=[];excluded=[];scope=Counter()
    for index,raw in enumerate(original):
        inputs,reason=coordinates(raw);row=dict(trace_index=index,raw=raw)
        if reason:row['reason']=reason;excluded.append(row);scope[reason]+=1
        else:row['inputs']=inputs;scoped.append(row);scope['scoped']+=1
    require(scoped==list(records(directory/'scoped.jsonl'))and excluded==list(records(directory/'excluded.jsonl')),'changed input/scope')
    ps=[dict(id=i,path=p,order=o,cw=c,inner_cw=w,mask=m)for i,(p,o,c,w,m)in enumerate(itertools.product(range(4),range(3),WORDS,WORDS,range(256)))]
    require(ps==list(records(directory/'profiles.jsonl'))and len(ps)==12288,'incomplete grammar')
    summary=json.loads((directory/'summary.json').read_text());require(summary['complete']is True,'partial run')
    instructions=Instructions(library);exact=Exact(instructions)
    killed=set();counts=Counter();executions=0;per_witness=Counter()
    def replay(identity,order,traced):
        p=ps[identity];row=scoped[order]
        g=grammar(len(row['inputs']),p['path'],p['order'])
        # Storing each reduction result also stores the final reduction node.
        mask=p['mask']|(2 if p['mask']&1 else 0)
        return exact.graph(g,row['inputs'],p['cw'],mask,inner_cw=p['inner_cw'],trace=traced)
    for record in records(directory/'first-kills.jsonl'):
        identity,order=record['candidate_id'],record['order']
        require(identity not in killed and 0<=identity<len(ps)and 0<=order<len(scoped),'duplicate/unknown program or witness')
        require(record['trace_index']==scoped[order]['trace_index']and record['expected']==scoped[order]['raw']['expected'],'changed observation')
        for previous in range(order):
            root=replay(identity,previous,False)[:10];executions+=1
            require(f'0x{published_bits(root):016x}'==scoped[previous]['raw']['expected'],'not first contradiction')
        roots,raw,values=replay(identity,order,True);executions+=1
        root=roots[:10];got=f'0x{published_bits(root):016x}'
        require(root.hex()==record['root']and raw.hex()==record['raw']and values.hex()==record['values'],'full-width instruction trace differs')
        require(got==record['got']and got!=record['expected'],'invalid first kill')
        counts['full_width_instruction_values']+=len(values)//10
        counts['certified_first_kills']+=1;killed.add(identity);per_witness[order]+=1
    alive=json.loads((directory/'survivors.json').read_text())
    require(len(set(alive))==len(alive)and killed.isdisjoint(alive)and killed|set(alive)==set(range(len(ps))),'unsearched program')
    for identity in alive:
        for order,row in enumerate(scoped):
            root=replay(identity,order,False)[:10];executions+=1
            require(f'0x{published_bits(root):016x}'==row['raw']['expected'],'invalid survivor')
    remaining=len(ps);witnesses=list(records(directory/'witnesses.jsonl'))
    require(len(witnesses)==len(scoped),'missing witness accounting')
    for order,(record,row)in enumerate(zip(witnesses,scoped)):
        remaining-=per_witness[order]
        require(record==dict(order=order,trace_index=row['trace_index'],killed=per_witness[order],remaining=remaining),'incorrect witness accounting')
    require(summary['programs']==len(ps)and summary['recorded']==len(original)and summary['scoped']==len(scoped)and
            summary['scope']==scope and summary['survivors']==len(alive)and summary['executions']==executions,'incorrect summary')
    require(all(digest(p)==h for p,h in hashes.items()),'checker changed during run')  # ubs:ignore — public integrity hashes, not secrets
    transcript=Path(str(output)+'.native-leaves.jsonl')
    with transcript.open('x')as out:
        for (kind,a,b,cw),value in sorted(instructions.transcript.items()):
            out.write(json.dumps(dict(kind=kind,a=a.hex(),b=b.hex(),cw=cw,value=value.hex()))+'\n')
    result=dict(complete=True,programs=len(ps),scoped=len(scoped),survivors=len(alive),executions=executions,
                counts=counts,rational_operations=exact.counts,native_leaves=instructions.counts,
                seconds=time.monotonic()-start,source_sha256=hashes,native_replay_sha256=digest(library),
                input_sha256={str(p):digest(p)for p in (Path(input_path),transcript,*sorted(directory.iterdir()))if p.is_file()},
                limits='Declared x87 grammar on the recorded native CPU only. Microcode leaves are replayed instructions, not correctly-rounded mathematical functions. No general Windows/CPU equivalence or complete Excel kernel identification.')
    with open(output,'x')as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(result))


if __name__=='__main__':check(*sys.argv[1:])
