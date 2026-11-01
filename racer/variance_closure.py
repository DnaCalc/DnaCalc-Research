"""Complete value closures for arbitrary-site three-input two-pass variance.

Every addition tree/permutation and every five-model arithmetic annotation is
included. Corrected two-pass bounds deliberately drop the dependency between
sum(delta) and sum(delta^2); an exclusion of that larger set is still valid.
No arithmetic value is inferred from an expected answer.
"""
from collections import Counter
from fractions import Fraction as Q
import itertools
import json
from pathlib import Path
import sys
import time

from bessel import operation,MODELS,pair,digest,require
from intervals import B64,Interval,number,pow2


def encode(values):return [pair(x)for x in sorted(values)]


def sums(leaves):
    """Subset closure includes both ordered children at every binary node."""
    table={1<<i:set(values)for i,values in enumerate(leaves)}
    for mask in range(1,8):
        if mask in table:continue
        result=set();left=(mask-1)&mask
        while left:
            right=mask^left
            if right:
                for a,b,model in itertools.product(table[left],table[right],MODELS):
                    result.add(operation('add',model,a,b))
            left=(left-1)&mask
        table[mask]=result
    return table


def serialized(table):return {str(mask):encode(table[mask])for mask in range(1,8)}


def close(row):
    def flatten(arg):return list(itertools.chain.from_iterable(flatten(x)for x in arg))if isinstance(arg,list)else[arg]
    args=flatten(row['args'])
    require(len(args)==3 and all(isinstance(x,str)and x.startswith('0x')for x in args), 'outside numeric three-input scope')
    xs=[Q(number(x))for x in args]
    require(all(x==0 or pow2(-100)<=abs(x)<=pow2(100)for x in xs), 'outside input magnitude bound')
    require(row['function']in('VAR','VAR.S','VAR.P','VARP','VARA','VARPA'), 'unknown variance function')
    population=row['function']in('VAR.P','VARP','VARPA')
    divisor=Q(3 if population else 2)
    total=sums([{x}for x in xs])
    means={operation('div',model,s,Q(3))for s,model in itertools.product(total[7],MODELS)}
    plain=set();corrected=set();details=[]
    for mean in sorted(means):
        ds=[{operation('add',model,x,-mean)for model in MODELS}for x in xs]
        squares=[{operation('mul',model,d,d)for d,model in itertools.product(values,MODELS)}for values in ds]
        q=sums(squares);t=sums(ds)
        correction={operation('div',md,operation('mul',mm,d,d),Q(3))
                    for d,mm,md in itertools.product(t[7],MODELS,MODELS)}
        outputs={operation('div',model,value,divisor)for value,model in itertools.product(q[7],MODELS)}
        adjusted={operation('div',md,operation('add',ms,value,-c),divisor)
                  for value,c,ms,md in itertools.product(q[7],correction,MODELS,MODELS)}
        details.append(dict(mean=pair(mean),deltas=list(map(encode,ds)),squares=list(map(encode,squares)),
                            square_sum=serialized(q),delta_sum=serialized(t),corrections=encode(correction),
                            plain_roots=encode(outputs),corrected_roots=encode(adjusted)))
        plain.update(outputs);corrected.update(adjusted)
    target=Q(number(row['expected']));span=B64.preimage(Interval(target,target))
    reachable=lambda roots:any(span.contains(value)for value in roots)
    return dict(raw=row,mean_sum=serialized(total),means=encode(means),per_mean=details,
                plain_roots=encode(plain),corrected_roots=encode(corrected),target_span=span.json(),
                plain_reachable=reachable(plain),corrected_reachable=reachable(corrected))


def run(input_path,destination):
    start=time.monotonic();dest=Path(destination);dest.mkdir();source=digest(Path(__file__))
    rows=[json.loads(line)for line in Path(input_path).read_text().splitlines()]
    counts=Counter()
    with(dest/'closures.jsonl').open('x')as out:
        for row in rows:
            proof=close(row);out.write(json.dumps(proof)+'\n');counts['observations']+=1
            counts['plain_unreachable']+=not proof['plain_reachable'];counts['corrected_unreachable']+=not proof['corrected_reachable']
            counts['mean_values']+=len(proof['means']);counts['plain_root_values']+=len(proof['plain_roots'])
            counts['corrected_root_values']+=len(proof['corrected_roots'])
    require(digest(Path(__file__))==source, 'source changed during run')  # ubs:ignore — public integrity hashes, not secrets
    result=dict(complete=True,counts=counts,seconds=time.monotonic()-start,
                source_sha256=source,input_sha256={str(input_path):digest(input_path)},
                claim='All binary sum trees/permutations and arbitrary five-model assignments at every arithmetic site for three-input plain two-pass variance. '
                      'Corrected two-pass uses a conservative dependency-free closure containing every consistent delta/square program.',
                limits='Three-input discovery observations and listed algebra only. Reachable roots need not share a reusable program; conservative corrected roots may not have consistent leaf choices. '
                       'Zero signs are immaterial to squares and division by positive integers. No branch, fit or complete Excel kernel follows.')
    with(dest/'summary.json').open('x')as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(result))


def online_close(row):
    def flatten(v):return list(itertools.chain.from_iterable(map(flatten,v)))if isinstance(v,list)else[v]
    bits=flatten(row['args']);require(len(bits)==3 and all(isinstance(v,str)and v.startswith('0x')for v in bits),'outside three-input online scope')
    xs=list(map(lambda v:Q(number(v)),bits));require(all(v==0 or pow2(-100)<=abs(v)<=pow2(100)for v in xs),'online magnitude scope')
    require(row['function']in('VAR','VAR.S','VAR.P','VARP','VARA','VARPA'),'wrong online variance function')
    divisor=Q(3 if row['function']in('VAR.P','VARP','VARPA')else 2)
    def values(kind,a,b):return{operation(kind,m,a,b)for m in MODELS}
    details=[];roots={'welford':set(),'youngs_cramer':set()}
    for perm in itertools.permutations(range(3)):
        record=dict(indices=list(perm))
        for family in roots:
            states={(xs[perm[0]],Q(0))};stages=[]
            for k,index in enumerate(perm[1:],2):
                x=xs[index];following=set()
                for accumulator,moment in states:
                    if family=='welford':
                        for delta in values('add',x,-accumulator):
                            means={v for q in values('div',delta,Q(k))for v in values('add',accumulator,q)}
                            for mean in means:
                                increments={v for e in values('add',x,-mean)for v in values('mul',delta,e)}
                                for t in increments:
                                    following.update((mean,q)for q in values('add',moment,t))
                    else:
                        for total in values('add',accumulator,x):
                            deviations={v for scaled in values('mul',Q(k),x)for v in values('add',scaled,-total)}
                            increments={v for delta in deviations for square in values('mul',delta,delta)for v in values('div',square,Q(k*(k-1)))}
                            for t in increments:following.update((total,q)for q in values('add',moment,t))
                states=following;stages.append([[pair(a),pair(b)]for a,b in sorted(states)])
            outputs={v for _,q in states for v in values('div',q,divisor)}
            record[family]=dict(stages=stages,roots=encode(outputs));roots[family].update(outputs)
        details.append(record)
    target=Q(number(row['expected']));span=B64.preimage(Interval(target,target))
    return dict(raw=row,permutations=details,roots={k:encode(v)for k,v in roots.items()},target_span=span.json(),reachable={k:any(span.contains(v)for v in values)for k,values in roots.items()})


def online_run(input_path,destination):
    start=time.monotonic();dest=Path(destination);dest.mkdir();source=digest(Path(__file__))
    rows=[json.loads(line)for line in Path(input_path).read_text().splitlines()];require(rows,'empty online observation scope')
    counts=Counter()
    with(dest/'closures.jsonl').open('x')as out:
        for row in rows:
            proof=online_close(row);out.write(json.dumps(proof)+'\n');counts['observations']+=1
            for family in('welford','youngs_cramer'):
                counts[family+'_unreachable']+=not proof['reachable'][family];counts[family+'_roots']+=len(proof['roots'][family])
                counts[family+'_state_pairs']+=sum(len(stage)for p in proof['permutations']for stage in p[family]['stages'])
    require(digest(Path(__file__))==source,'source changed during online closure')  # ubs:ignore — public integrity hashes, not secrets
    summary=dict(complete=True,counts=counts,seconds=time.monotonic()-start,source_sha256=source,input_sha256={str(input_path):digest(input_path)},
                 limits='Three-input Welford and Youngs-Cramer unrolled updates, all6 permutations and independent five-model annotations at every site. Boundary states preserve mean/sum-second-moment correlations. Reachability is not a reusable producer, fitted value or product gain.')
    with(dest/'summary.json').open('x')as out:json.dump(summary,out,indent=2);out.write('\n')
    print(json.dumps(summary))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='online':online_run(*sys.argv[2:])
    else:run(*sys.argv[1:])
