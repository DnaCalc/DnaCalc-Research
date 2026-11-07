"""Separate adjacent constants of a frozen arithmetic shape before asking Excel.

This is a question selector, not a coefficient fitter. It reuses the x87 DAG
interpreter and the audited scalar service. No separator on a finite search means
undetermined within that search, never a uniqueness or equivalence proof. ERFC's
stored positive Horner stages additionally admit exact absorption certificates.
"""
import argparse
import copy
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import random
import re
import sys

import numpy as np

try:
    from .x87_programs import Builder, Native, build_library
    from .fingerprint import bits, value, CachedService, load_service, write, sha
except ImportError:
    from x87_programs import Builder, Native, build_library
    from fingerprint import bits, value, CachedService, load_service, write, sha

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'checker'))
from exact import BINARY64 as B, EXTENDED64 as E, binary64, two, exponent_of


def adjacent(raw, direction):
    x = value(raw)
    if not math.isfinite(x):
        raise ValueError('finite constant required')
    y = math.nextafter(x, math.inf if direction > 0 else -math.inf)
    if not math.isfinite(y):
        raise ValueError('neighbor is nonfinite')
    return bits(y)


def changed_shape(shape, name, raw):
    result = copy.deepcopy(shape['graph'])
    result['nodes'][shape['constants'][name]]['bits'] = raw
    return result


def probe(shape, inputs, native, keep=12):
    """Keep bounded separating inputs for both immediate neighbors of each constant."""
    graph = shape['graph']; mask = shape['mask']; cw = shape['cw']
    if not len(inputs) or keep<1:
        raise ValueError('nonempty search grid and a positive question budget required')
    if not shape['constants']:
        raise ValueError('no named constants in the shape')
    for name,node in shape['constants'].items():
        if not isinstance(node,int)or not 0<=node<len(graph['nodes'])or graph['nodes'][node]['op']!='const':
            raise ValueError('named constant is not a constant node: '+name)
    baseline = native.grid(graph, inputs, cw, mask)
    result = []
    for name, node in shape['constants'].items():
        raw = graph['nodes'][node]['bits']
        alternatives = []
        for direction in (-1, 1):
            neighbor = adjacent(raw, direction)
            other = native.grid(changed_shape(shape, name, neighbor), inputs, cw, mask)
            sites = np.flatnonzero(other[:, 0] != baseline[:, 0])
            alternatives.append(dict(bits=neighbor, direction=direction,
                status='separating_questions' if len(sites) else 'undetermined_within_search',
                separating_inputs=len(sites), examples=[dict(
                    inputs=[f'{int(v):016x}' for v in inputs[i]],
                    baseline=f'{int(baseline[i,0]):016x}',
                    alternative=f'{int(other[i,0]):016x}') for i in sites[:keep]]))
        result.append(dict(name=name, bits=raw, neighbors=alternatives,
            status='separating_questions' if any(a['separating_inputs'] for a in alternatives)
                   else 'undetermined_within_search'))
        print(name, result[-1]['status'], [a['separating_inputs'] for a in alternatives], flush=True)
    return result


def service_questions(shape, constants):
    """Write a reusable numeric scalar-service batch, without predicted answers."""
    selected={}
    columns=shape['argument_inputs']
    width=1+max((n['a']for n in shape['graph']['nodes']if n['op']=='input'),default=-1)
    if not isinstance(columns,list)or not columns or any(type(i)is not int or not 0<=i<width for i in columns):
        raise ValueError('argument inputs must name nonnegative graph input columns')
    for c in constants:
        for alternative in c['neighbors']:
            for r in alternative['examples']:
                args=[dict(type='number',bits=r['inputs'][i])for i in columns]
                selected.setdefault(json.dumps(args,sort_keys=True),args)
    return dict(function=shape['function'],rows=list(selected.values()))


def erfc_shape(checkout):
    source = Path(checkout)/'crates/calc-func/src/ported/normal_dist_common.rs'
    text = source.read_text(); tables = {}
    for name in ('N', 'D'):
        body = re.search(r'const ERFC_TAIL_'+name+r':[^=]+=(.*?);', text, re.S)
        if body is None:
            raise ValueError('missing recovered ERFC table')
        tables[name] = [m.replace('_', '')[2:] for m in re.findall(r'0x[0-9a-f_]+', body[1])]
    g = Builder(2); constants = {}; stages = {}
    def horner(name):
        table = tables[name]
        nodes = [g.constant(raw) for raw in table]
        constants.update({name+str(i): node for i, node in enumerate(nodes)})
        acc = nodes[-1]
        for i in range(len(nodes)-2, -1, -1):
            product = g.op('mul', acc, 0)
            acc = g.op('add', product, nodes[i])
            stages[name+str(i)] = product
        return acc
    n, d = horner('N'), horner('D')
    f = g.op('div', n, d); out = g.op('mul', 1, f)
    # D0 fixes the scale and is not a recovered free coefficient.
    del constants['D0']
    return dict(format='dnacalc-plateau-shape/1', function='ERFC.PRECISE',argument_inputs=[0],
        graph=g.graph(out, f), mask=(1 << g.sites)-1, cw=0x133f,
        constants=constants, tables=tables, stages=stages,
        source_bindings={str(source.resolve()): sha(source)},
        domain=dict(low=bits(math.nextafter(4., math.inf)), high=bits(26.5)),
        leaf='Input1 is independently observed EXP(-storedPC64(x*x)); the search uses host EXP only to propose questions.')


def erfc_inputs(shape, count, seed):
    """Mapping grid, single-bit windows, exponent edges and Horner binade edges."""
    low, high = (value(shape['domain'][k]) for k in ('low', 'high'))
    rng = random.Random(seed)
    xs = {bits(low), bits(high)}
    anchors = [4., 4.25, 5., 5.3, 6., 8., 12.3, 16., 20.3, 26.5]
    tables = {n: list(map(value, t)) for n,t in shape['tables'].items()}
    # Positive Horner product stages are monotone. Locate every power-of-two
    # crossing independently of observations or current match counts.
    for table in tables.values():
        def products(x):
            a = table[-1]; out = []
            for c in table[-2::-1]:
                p = a*x; out.append(p); a = p+c
            return out
        a, b = products(low), products(high)
        for i, (left, right) in enumerate(zip(a,b)):
            for e in range(math.ceil(math.log2(left)), math.floor(math.log2(right))+1):
                lo, hi = low, high
                for _ in range(56):
                    mid = (lo+hi)/2
                    if products(mid)[i] < 2.**e: lo = mid
                    else: hi = mid
                anchors.append((lo+hi)/2)
    for center in anchors:
        raw = int(bits(center),16)
        for offset in range(-64,65):
            x = value(f'{raw+offset:016x}')
            if low <= x <= high: xs.add(bits(x))
    for i in range(count):
        xs.add(bits(low+(high-low)*i/max(1,count-1)))
        xs.add(f'{rng.randrange(int(bits(low),16),int(bits(high),16)+1):016x}')
    inputs = []
    for raw in sorted(xs):
        x = value(raw)
        z = float(E.round(Q(x)*Q(x)))
        inputs.append((int(raw,16), int(bits(math.exp(-z)),16)))
    return np.array(inputs,dtype=np.uint64)


def absorption_certificate(shape, name, neighbor):
    """Prove a neighboring load is absorbed at its first stored Horner sum.

    Positive coefficients make endpoint propagation an enclosure for every x
    in the declared domain. Within fixed product/sum binades, translating a
    stored binary64 product by four coarser ULPs preserves both nearest/even
    rounding parities. Enumerate those residue classes, plus all binade cuts.
    This overapproximates the attainable products; a failure proves nothing.
    """
    prefix, index = name[0], int(name[1:])
    table = [binary64(x) for x in shape['tables'][prefix]]
    domain=[binary64(shape['domain'][k])for k in ('low','high')]
    other=binary64(neighbor)
    if (not 0<=index<len(table)or any(c<=0 for c in table)or other<=0
            or not 0<domain[0]<=domain[1]):
        raise ValueError('absorption proof requires positive coefficients and an ordered positive domain')
    if index == len(table)-1:
        return dict(proved=False, reason='leading coefficient is not an added constant')
    def store(x): return B.round(E.round(x))
    def bound(x):
        acc = table[-1]
        for i in range(len(table)-2, index-1, -1):
            p = store(acc*x)
            if i == index: return p
            acc = store(p+table[i])
    lo, hi = [bound(x)for x in domain]
    if not 0<lo<=hi:
        return dict(proved=False,reason='positive stored product enclosure required')
    constants = (table[index],other); checked = 0
    for e in range(exponent_of(lo), exponent_of(hi)+1):
        h = two(e-52)
        kmin = max(1<<52, math.ceil(lo/h)); kmax = min((1<<53)-1, math.floor(hi/h))
        if kmin > kmax: continue
        cuts = {kmin,kmax+1}
        for c in constants:
            for s in range(exponent_of(kmin*h+c), exponent_of(kmax*h+c)+2):
                k = math.ceil((two(s)-c)/h)
                if kmin < k <= kmax: cuts.add(k)
        cuts = sorted(cuts)
        for left,right in zip(cuts,cuts[1:]):
            # The common period must cover the larger publication lattice,
            # with an extra factor two for its tie parity.
            s = max(exponent_of(left*h+c) for c in constants)
            period = max(4,1 << max(0,s-e+1))
            if period > 4096:
                return dict(proved=False,reason='residue period exceeds budget')
            for k in range(left,min(right,left+period)):
                checked += 1
                a,b = (store(k*h+c) for c in constants)
                if a != b:
                    return dict(proved=False,product=bits(float(k*h)),baseline=bits(float(a)),
                                alternative=bits(float(b)),residue_checks=checked)
    if not checked:
        return dict(proved=False,reason='no product residue checked')
    return dict(proved=True,domain=shape['domain'],product_bounds=[bits(float(lo)),bits(float(hi))],
                residue_checks=checked,reason='neighbor absorbed at first stored sum for every possible binary64 product in enclosure')


def erfc(checkout, output, cache, count=30000, seed=35):
    output=Path(output); output.mkdir(parents=True,exist_ok=False)
    library=output/'native.so'; build_library(library); native=Native(library)
    shape=erfc_shape(checkout); write(output/'shape.json',shape)
    inputs=erfc_inputs(shape,count,seed)
    constants=probe(shape,inputs,native)
    selected={}
    for c in constants:
        for a in c['neighbors']:
            a['absorption']=absorption_certificate(shape,c['name'],a['bits'])
            if a['absorption']['proved']:a['status']='undetermined_proved_absorbed'
            for row in a['examples']: selected[row['inputs'][0]]=None
        if all(a['absorption']['proved'] for a in c['neighbors']):
            c['status']='undetermined_proved_absorbed'
        elif any(a['status'].startswith('undetermined')for a in c['neighbors']):
            c['status']='partially_separated_constant_remains_undetermined'
    # E is a conditional measured leaf. Read it independently before freezing
    # the public ERFC predictions; host EXP never supplies a claimed answer.
    client,_=load_service(checkout); service=CachedService(client.Client(timeout=180),cache)
    xs=list(selected)
    zs=[bits(-float(E.round(Q(value(x))*Q(value(x))))) for x in xs]
    exp=service.ask('EXP',[[client.number(z)] for z in zs])['rows'] if xs else []
    if any(a['type']!='number' for a in exp): raise ValueError('non-numeric EXP leaf')
    measured=np.array([[int(x,16),int(a['bits'],16)] for x,a in zip(xs,exp)],dtype=np.uint64).reshape((-1,2))
    question_constants=probe(shape,measured,native,keep=len(xs)) if xs else []
    # Every emitted question still separates a public result with the actual E.
    questions={r['inputs'][0] for c in question_constants for a in c['neighbors'] for r in a['examples']}
    batch=[dict(input=x,arguments=[client.number(x)]) for x in xs if x in questions]
    predictions=dict(format='dnacalc-plateau-questions/1',function=shape['function'],
        rows=[r['arguments']for r in batch],shape_sha256=sha(output/'shape.json'),
        constants=question_constants,questions=batch,leaf_receipts=service.receipts.copy())
    write(output/'questions.json',predictions)
    answers=service.ask('ERFC.PRECISE',[r['arguments'] for r in batch])['rows'] if batch else []
    observed={r['input']:a for r,a in zip(batch,answers)}
    for c in question_constants:
        for a in c['neighbors']:
            rows=[r for r in a['examples'] if r['inputs'][0] in observed]
            a['baseline_agrees']=sum(observed[r['inputs'][0]]==dict(type='number',bits=r['baseline']) for r in rows)
            a['alternative_agrees']=sum(observed[r['inputs'][0]]==dict(type='number',bits=r['alternative']) for r in rows)
            a['fresh_answers']=len(rows)
    result=dict(format='dnacalc-plateau-erfc/1',complete=True,search_inputs=len(inputs),seed=seed,
        constants=constants,fresh_comparison=question_constants,
        observations=[dict(input=r['input'],answer=a) for r,a in zip(batch,answers)],
        receipts=service.receipts,source_bindings={**shape['source_bindings'],
            str(Path(__file__).resolve()):sha(__file__),str((ROOT/'racer/x87_native.c').resolve()):sha(ROOT/'racer/x87_native.c')},
        limits='One constant moves at a time. Finite no-separator results remain undetermined within this search. Absorption proofs concern only the stated x>4 through26.5 shape and domain. No constant is chosen by match percentage; no universal coefficient uniqueness is claimed. Public predictions condition on independently observed EXP leaves.')
    write(output/'result.json',result)
    print(json.dumps(dict(search_inputs=len(inputs),questions=len(batch),statuses={c['name']:c['status'] for c in constants})),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('erfc');p.add_argument('checkout');p.add_argument('output');p.add_argument('cache')
    p.add_argument('--count',type=int,default=30000);p.add_argument('--seed',type=int,default=35)
    p=sub.add_parser('probe');p.add_argument('shape');p.add_argument('inputs');p.add_argument('library');p.add_argument('output')
    args=parser.parse_args()
    if args.command=='erfc':erfc(args.checkout,args.output,args.cache,args.count,args.seed)
    else:
        if Path(args.output).exists():
            raise FileExistsError('plateau questions need a new output path')
        shape=json.loads(Path(args.shape).read_text())
        inputs=np.array([[int(x,16)for x in row]for row in json.loads(Path(args.inputs).read_text())],dtype=np.uint64)
        constants=probe(shape,inputs,Native(args.library))
        write(args.output,dict(constants=constants,**service_questions(shape,constants),search_inputs=len(inputs),
            shape_sha256=sha(args.shape),inputs_sha256=sha(args.inputs)))


if __name__=='__main__':main()
