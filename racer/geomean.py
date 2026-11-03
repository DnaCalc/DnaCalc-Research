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


def exact_root(product, count):
    """A correctly rounded root, certified by rational midpoint powers."""
    import math
    # The approximation only starts the search. It never determines an answer
    # without exact comparisons at both rounding boundaries.
    candidate = int.from_bytes(struct.pack('>d', math.pow(float(product), 1.0/count)), 'big')
    def value(bits): return Q(struct.unpack('>d', bits.to_bytes(8, 'big'))[0])
    for _ in range(1024):
        center = value(candidate)
        low = ((value(candidate-1)+center)/2)**count
        high = ((value(candidate+1)+center)/2)**count
        even = candidate % 2 == 0
        if product < low or (product == low and not even): candidate -= 1
        elif product > high or (product == high and not even): candidate += 1
        else: return f'0x{candidate:016x}'
    raise ValueError('initial root approximation failed to bracket; no result emitted')


def exact_root_run(input_path, destination):
    """Eliminate fixed rounded-product/correct-root programs, not native log approximations."""
    from bessel import MODELS, operation, pair
    dest = Path(destination); dest.mkdir()
    original = [json.loads(line) for line in Path(input_path).read_text().splitlines()]
    rows, scope = [], Counter()
    for index, raw in enumerate(original):
        inputs, reason = coordinates(raw)
        if reason: scope[reason] += 1
        else:
            rows.append(dict(trace_index=index, raw=raw, inputs=inputs)); scope['scoped'] += 1
    require(rows, 'empty exact-root scope')
    ps = [dict(id=i, order=order, model=model) for i,(order,model) in
          enumerate(itertools.product(range(3), (*MODELS, 'exact')))]
    def calculate(profile, row):
        xs = [Q(struct.unpack('>d', x.to_bytes(8, 'big'))[0]) for x in row['inputs']]
        def multiply(a,b): return a*b if profile['model']=='exact' else operation('mul', profile['model'], a, b)
        product = fold(xs, profile['order'], multiply, Q(1))
        return exact_root(product, len(xs)), product
    for name, data in [('profiles.json', ps), ('scoped.json', rows)]:
        (dest/name).write_text(json.dumps(data)+'\n')
    alive = list(range(len(ps))); executions = 0
    with (dest/'first-kills.jsonl').open('x') as kills:
        for order, row in enumerate(rows):
            following = []
            for identity in alive:
                got, product = calculate(ps[identity], row); executions += 1
                if got == row['raw']['expected']: following.append(identity)
                else: kills.write(json.dumps(dict(candidate_id=identity, order=order, got=got,
                            product=pair(product), expected=row['raw']['expected']))+'\n')
            alive = following
    (dest/'survivors.json').write_text(json.dumps(alive)+'\n')
    result = dict(complete=True, programs=len(ps), observations=len(rows), recorded=len(original),
        scope=scope, survivors=len(alive), executions=executions,
        input_sha256={str(input_path):digest(input_path)}, source_sha256=digest(Path(__file__)),
        limits='Three product orders, five fixed rounding/store policies plus exact products, followed by a correctly rounded root proved by exact rational powers. This is a specified alternative program family, not a claim that Excel uses a correctly rounded root. All old observations are discovery; no fitting or per-row choice.')
    (dest/'summary.json').write_text(json.dumps(result, indent=2)+'\n'); print(json.dumps(result))


def product_closure(raw):
    """All product trees for two or three modest positive operands."""
    from bessel import MODELS, operation, pair
    from intervals import B64, Interval, number, pow2
    inputs, reason = coordinates(raw)
    require(reason is None and 2 <= len(inputs) <= 3, 'product-closure arity/domain')
    xs = [Q(number(f'0x{x:016x}')) for x in inputs]
    require(all(pow2(-100) <= x <= pow2(100) for x in xs), 'product-closure magnitude')
    table = {1 << i: {x} for i, x in enumerate(xs)}
    for mask in range(1, 1 << len(xs)):
        if mask in table:
            continue
        values = set()
        left = (mask-1) & mask
        while left:
            right = mask ^ left
            for a, b in itertools.product(table[left], table[right]):
                values.add(a*b)
                values.update(operation('mul', m, a, b) for m in MODELS)
            left = (left-1) & mask
        table[mask] = values
    target = Q(number(raw['expected']))
    span = B64.preimage(Interval(target, target))
    # Root is monotone: pull its publication interval back by exact nth powers.
    lo, hi = span.lo**len(xs), span.hi**len(xs)
    roots = table[(1 << len(xs))-1]
    reachable = any((p > lo or p == lo and span.lc) and
                    (p < hi or p == hi and span.hc) for p in roots)
    return dict(raw=raw, products={str(k):[pair(x) for x in sorted(v)] for k,v in table.items()},
                root_span=span.json(), reachable=reachable)


def product_closure_run(input_path, destination):
    from intervals import number, pow2
    started = time.monotonic()
    dest = Path(destination)
    dest.mkdir()
    frozen = {str(p):digest(p) for p in (Path(input_path), Path(__file__),
        Path(__file__).with_name('bessel.py'), Path(__file__).with_name('intervals.py'))}
    counts = Counter()
    with Path(input_path).open() as source, (dest/'closures.jsonl').open('x') as out, (dest/'scope.jsonl').open('x') as scope:
        for index, line in enumerate(source):
            raw = json.loads(line)
            inputs, reason = coordinates(raw)
            if reason is None and not 2 <= len(inputs) <= 3:
                reason = 'arity_outside_2_to_3'
            if reason is None and any(not pow2(-100) <= Q(number(v)) <= pow2(100) for v in flatten(raw['args'])):
                reason = 'outside_modest_magnitudes'
            if reason is None and number(raw['expected']) <= 0:
                reason = 'nonpositive_expected'
            scope.write(json.dumps(dict(trace_index=index, reason=reason))+'\n')
            if reason:
                counts[reason] += 1
                continue
            proof = product_closure(raw)
            out.write(json.dumps(dict(trace_index=index, **proof))+'\n')
            counts['scoped'] += 1
            counts['unreachable'] += not proof['reachable']
            counts['root_product_values'] += len(proof['products'][str((1 << len(inputs))-1)])
    require(counts['scoped'] and frozen == {p:digest(p) for p in frozen}, 'empty or changed product closure')
    result = dict(complete=True, counts=counts, input_sha256=frozen, seconds=time.monotonic()-started,
        limits='All multiplication trees/permutations with an independent choice at every node among five precision/store models or exact rational multiplication, followed by a correctly rounded real nth root. '
               'Only two/three positive operands in[2^-100,2^100]. Exact inverse midpoint powers determine reachability. '
               'No native logarithm/power leaf, complete Excel algorithm, tolerance or product gain is claimed.')
    (dest/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='product-closure':
        product_closure_run(*sys.argv[2:]); raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='exact-root':
        exact_root_run(*sys.argv[2:]); raise SystemExit(0)
    parser=argparse.ArgumentParser();parser.add_argument('input');parser.add_argument('library');parser.add_argument('destination')
    args=parser.parse_args();run(args.input,args.library,args.destination)
