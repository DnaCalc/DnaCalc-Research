"""Finite single-predictor least-squares identities, with exact rounding.

Original algebra, MIT. This extends the existing arithmetic research machinery;
it is not an alternate application evaluator. No coefficient or rank tolerance
is selected or fitted. Multiple-predictor and LOGEST programs are outside scope.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
import itertools
import json
from pathlib import Path
import struct
import time

from bessel import MODELS,operation,pair,digest,require
from intervals import B64,number,pow2


def extract(row):
    if row['function']!='LINEST' or row['result_index']not in ([1,1],[1,2]):return None,'surface_or_selected_statistic'
    a=row['args']
    if len(a)<2:return None,'implicit_predictor'
    def vector(arg):
        if not isinstance(arg,list):raise ValueError('not an array')
        if all(isinstance(r,list)for r in arg):
            if len(arg)>1 and any(len(r)!=1 for r in arg):raise ValueError('multiple predictors')
            arg=list(itertools.chain.from_iterable(arg))
        if not all(isinstance(x,str)and x.startswith('0x')for x in arg):raise ValueError('non-numeric vector')
        return [Q(number(x))for x in arg]
    try:y,x=map(vector,a[:2])
    except ValueError:return None,'outside_numeric_single_predictor'
    if len(y)!=len(x)or not(2<=len(x)<=32):return None,'shape_or_length'
    if any(v and not pow2(-100)<=abs(v)<=pow2(100)for v in x+y):return None,'magnitude_bound'
    if len(a)<3:constant=True
    elif isinstance(a[2],bool):constant=a[2]
    elif isinstance(a[2],str)and a[2].startswith('0x'):constant=number(a[2])!=0
    else:return None,'const_coercion_outside_scope'
    if (constant and len(set(x))==1)or(not constant and all(v==0 for v in x)):return None,'not_exact_full_rank'
    if not row['expected'].startswith('0x'):return None,'non_numeric_expected'
    return (x,y,constant),None


def profiles():
    for identity,(algorithm,linear,squared,models)in enumerate(itertools.product(range(3),range(3),range(3),itertools.product(MODELS,repeat=5))):
        yield dict(id=identity,algorithm=algorithm,linear_order=linear,square_order=squared,models=list(models))


def evaluate(profile,values,column):
    x,y,constant=values;m=profile['models'];algorithm=profile['algorithm']
    def op(kind,site,a,b):return operation(kind,m[site],a,b)
    def sub(site,a,b):return op('add',site,a,-b)
    def reduce(items,order,site):
        if order==2:
            if not items:return Q(0)
            if len(items)==1:return items[0]
            middle=len(items)//2
            return op('add',site,reduce(items[:middle],2,site),reduce(items[middle:],2,site))
        total=Q(0)
        for value in items if order==0 else items[::-1]:total=op('add',site,total,value)
        return total
    n=Q(len(x));linear=profile['linear_order'];squared=profile['square_order']
    sx,sy=reduce(x,linear,0),reduce(y,linear,0)
    mx,my=op('div',1,sx,n),op('div',1,sy,n)
    if algorithm==0 and constant:
        dx=[sub(1,v,mx)for v in x];dy=[sub(1,v,my)for v in y]
        xx=[op('mul',2,v,v)for v in dx];xy=[op('mul',2,v,w)for v,w in zip(dx,dy)]
        denominator=reduce(xx,squared,3);numerator=reduce(xy,squared,3)
    else:
        xx=[op('mul',2,v,v)for v in x];xy=[op('mul',2,v,w)for v,w in zip(x,y)]
        denominator=reduce(xx,squared,3);numerator=reduce(xy,squared,3)
        if constant:
            if algorithm==1:
                denominator=sub(4,denominator,op('mul',2,sx,mx))
                numerator=sub(4,numerator,op('mul',2,sx,my))
            else:
                denominator=sub(4,op('mul',4,n,denominator),op('mul',2,sx,sx))
                numerator=sub(4,op('mul',4,n,numerator),op('mul',2,sx,sy))
    if denominator==0:return None,'error:Div0'
    slope=op('div',4,numerator,denominator)
    value=slope if column==1 else sub(4,my,op('mul',2,slope,mx))if constant else Q(0)
    published=B64.rn(value)
    if published is None:return value,'error:Num'
    if published==0:return value,'0x0000000000000000'
    if abs(published)<pow2(-1022):return value,'fault:NotImplemented'
    return value,'0x'+struct.pack('>d',float(published)).hex()


def run(input_path,destination):
    start=time.monotonic();dest=Path(destination);dest.mkdir()
    original=[json.loads(line)for line in Path(input_path).read_text().splitlines()]
    scope=Counter();rows=[];values=[]
    with(dest/'scoped.jsonl').open('x')as good,(dest/'excluded.jsonl').open('x')as bad:
        for index,row in enumerate(original):
            coordinates,reason=extract(row)
            if reason:bad.write(json.dumps(dict(trace_index=index,raw=row,reason=reason))+'\n');scope[reason]+=1
            else:
                rows.append((index,row));values.append(coordinates);scope['scoped']+=1
                good.write(json.dumps(dict(trace_index=index,raw=row))+'\n')
    ps=list(profiles());require(len(ps)==84375 and rows, 'incomplete grammar or empty regression scope')
    with(dest/'profiles.jsonl').open('x')as out:
        for p in ps:out.write(json.dumps(p)+'\n')
    alive=list(range(len(ps)));executions=0;source=digest(Path(__file__))
    with(dest/'first-kills.jsonl').open('x')as kills,(dest/'witnesses.jsonl').open('x')as witnesses:
        for order,((index,row),data)in enumerate(zip(rows,values)):
            before=len(alive);next_alive=[]
            for identity in alive:
                value,got=evaluate(ps[identity],data,row['result_index'][1]);executions+=1
                if got==row['expected']:next_alive.append(identity)
                else:kills.write(json.dumps(dict(candidate_id=identity,order=order,trace_index=index,
                     root=pair(value)if value is not None else None,got=got,expected=row['expected']))+'\n')
            alive=next_alive
            witnesses.write(json.dumps(dict(order=order,trace_index=index,killed=before-len(alive),remaining=len(alive)))+'\n')
    with(dest/'survivors.json').open('x')as out:json.dump(alive,out);out.write('\n')
    require(digest(Path(__file__))==source, 'source changed during run')  # ubs:ignore — public integrity hashes, not secrets
    result=dict(complete=True,programs=len(ps),scoped=len(rows),recorded=len(original),scope=scope,
                executions=executions,survivors=len(alive),seconds=time.monotonic()-start,
                input_sha256={str(input_path):digest(input_path)},source_sha256=source,
                claim='Finite centered/raw-moment single-predictor least-squares coefficient programs with coherent stores and precisions. No coefficient is fitted.',
                limits='LINEST coefficient/intercept publications, exact full-rank numeric vectors and fixed algebra only. '
                       'No QR, multiple predictor, standard-error/statistic, LOGEST, rank-tolerance or complete kernel identification.')
    with(dest/'summary.json').open('x')as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(result))


def rank_coordinates(raw):
    """Read only explicit degrees-of-freedom observations, with exact numeric inputs."""
    if raw['function'] not in ('LINEST','LOGEST') or raw['result_index'] != [4,2]:
        return None, 'not_degrees_of_freedom'
    args = raw['args']
    if len(args)<2 or not isinstance(args[0],list) or not isinstance(args[1],list):
        return None, 'not_explicit_numeric_arrays'
    try:
        y = args[0]
        if all(isinstance(r,list) for r in y):
            if not (len(y)==1 or all(len(r)==1 for r in y)):return None,'response_shape'
            y = list(itertools.chain.from_iterable(y))
        x = args[1]
        if not all(isinstance(r,list) for r in x): x = [[v] for v in x]
        if len(x)!=len(y) or not(2<=len(y)<=32) or not x or not(1<=len(x[0])<=8):
            return None,'matrix_shape'
        if any(len(r)!=len(x[0]) for r in x):return None,'matrix_shape'
        values = [Q(number(v)) for v in y+list(itertools.chain.from_iterable(x))]
        if any(v and not pow2(-100)<=abs(v)<=pow2(100) for v in values):return None,'magnitude_bound'
        columns = [[Q(number(row[c])) for row in x] for c in range(len(x[0]))]
        constant = True if len(args)<3 else args[2] if isinstance(args[2],bool) else number(args[2])!=0
        expected = Q(number(raw['expected']))
        if expected.denominator!=1 or not 0<=expected<=len(y):return None,'not_integer_degrees_of_freedom'
    except (ValueError,TypeError,OverflowError,struct.error):return None,'non_numeric'
    return (columns,constant,int(expected)),None


def rank_profiles():
    # All tolerances are mathematical IEEE roundoff scales, never inferred
    # from a residual or fitted to an Excel rank transition.
    for identity,(threshold,normalizer,zero) in enumerate(itertools.product(
            ('zero','epsilon','epsilon_squared','sqrt_epsilon','gamma_n','gamma_n_squared'),
            ('column','largest','first','absolute'),('keep','drop'))):
        yield dict(id=identity,threshold=threshold,normalizer=normalizer,zero=zero)


def rank_evaluate(profile, data):
    columns, constant, _ = data
    n = len(columns[0]); epsilon=Q(1,2**52);gamma=n*epsilon/(1-n*epsilon)
    tolerances = dict(zero=Q(0),epsilon=epsilon,epsilon_squared=epsilon**2,
                      sqrt_epsilon=Q(1,2**26),gamma_n=gamma,gamma_n_squared=gamma**2)
    values = [[x-sum(c)/n for x in c] if constant else list(c) for c in columns]
    initial = [sum(x*x for x in c) for c in values]
    live = list(range(len(columns)));basis=[];first=None;kept=[];dropped=[]
    while live and len(kept)+int(constant)<n:
        norms = {i:sum(x*x for x in values[i]) for i in live}
        index=max(live,key=lambda i:norms[i]); norm=norms[index];live.remove(index)
        if first is None:first=norm
        scale = {'column':initial[index], 'largest':max(initial), 'first':first, 'absolute':Q(1)}[profile['normalizer']]
        if initial[index]==0:
            if profile['zero']=='keep':kept.append(index)
            else:dropped.append(index)
            continue
        if norm<=tolerances[profile['threshold']]*scale:
            dropped.append(index);continue
        kept.append(index);basis.append(index)
        for other in live:
            coefficient=sum(a*b for a,b in zip(values[index],values[other]))/norm
            values[other]=[b-coefficient*a for a,b in zip(values[index],values[other])]
    dropped.extend(live)
    return n-len(kept)-int(constant), dict(kept=kept, dropped=dropped)


def rank_run(input_path, destination):
    dest=Path(destination);dest.mkdir();rows=[];scope=Counter()
    for index,raw in enumerate(map(json.loads,Path(input_path).read_text().splitlines())):
        data,reason=rank_coordinates(raw)
        if reason:scope[reason]+=1
        else:rows.append((index,raw,data));scope['scoped']+=1
    require(rows,'empty rank-policy scope')
    ps=list(rank_profiles());alive=list(range(len(ps)));executions=0
    (dest/'profiles.json').write_text(json.dumps(ps)+'\n')
    (dest/'scoped.json').write_text(json.dumps([dict(trace_index=i,raw=r) for i,r,_ in rows])+'\n')
    with (dest/'first-kills.jsonl').open('x')as kills:
        for order,(index,raw,data)in enumerate(rows):
            following=[]
            for identity in alive:
                got,trace=rank_evaluate(ps[identity],data);executions+=1
                if got==data[2]:following.append(identity)
                else:kills.write(json.dumps(dict(candidate_id=identity,order=order,trace_index=index,
                     got=got,expected=data[2],trace=trace))+'\n')
            alive=following
    (dest/'survivors.json').write_text(json.dumps(alive)+'\n')
    result=dict(complete=True,programs=len(ps),observations=len(rows),survivors=len(alive),executions=executions,
        scope=scope,source_sha256=digest(Path(__file__)),input_sha256={str(input_path):digest(input_path)},
        limits='Exact rational pivoted Gram-Schmidt geometry with fixed IEEE roundoff-scale comparisons. Degrees of freedom only; no floating QR implementation, coefficient program, fitted tolerance or reusable Excel rank rule is identified. Surviving necessary rank policies cannot be adopted as numeric programs.')
    (dest/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='rank':rank_run(*sys.argv[2:]);raise SystemExit(0)
    parser=argparse.ArgumentParser();parser.add_argument('input');parser.add_argument('destination')
    args=parser.parse_args();run(args.input,args.destination)
