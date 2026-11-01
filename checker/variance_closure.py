"""Independent full-set proof using explicit permutations/trees, not subset DP."""
from collections import Counter
from fractions import Fraction as R
from functools import lru_cache
import itertools
import json
from pathlib import Path
import sys
import time

from exact import BINARY64 as B,Span,binary64,operation,two
from check import records,require,digest

MODELS=('strict','stored_pc64','resident_pc64','stored_pc53','resident_pc53')


@lru_cache(maxsize=262144)
def op(kind,model,a,b):return operation(a,b,kind,model)


def encode(values):return [[str(x.numerator),str(x.denominator)]for x in sorted(values)]


def all_sums(leaves):
    table={}
    for mask in range(1,8):
        indices=[i for i in range(3)if mask&(1<<i)];values=set()
        for perm in itertools.permutations(indices):
            selected=[leaves[i]for i in perm]
            if len(perm)==1:values.update(selected[0])
            elif len(perm)==2:
                for a,b,m in itertools.product(*selected,MODELS):values.add(op('add',m,a,b))
            else:
                for a,b,c,first,last in itertools.product(*selected,MODELS,MODELS):
                    values.add(op('add',last,op('add',first,a,b),c))
                    values.add(op('add',last,a,op('add',first,b,c)))
        table[str(mask)]=values
    return table


def encoded(table):return {mask:encode(values)for mask,values in table.items()}


def verify(row,proof):
    require(proof['raw']==row,'changed input, source or answer')
    def flatten(v):return list(itertools.chain.from_iterable(map(flatten,v)))if isinstance(v,list)else[v]
    bits=flatten(row['args']);require(len(bits)==3,'outside three-input scope')
    xs=list(map(binary64,bits));require(all(x==0 or two(-100)<=abs(x)<=two(100)for x in xs),'magnitude scope')
    require(row['function']in('VAR','VAR.S','VAR.P','VARP','VARA','VARPA'),'wrong function')
    divisor=R(3)if row['function']in('VAR.P','VARP','VARPA')else R(2)
    total=all_sums([{x}for x in xs]);require(encoded(total)==proof['mean_sum'],'incomplete mean-sum closure')
    means={op('div',m,s,R(3))for s,m in itertools.product(total['7'],MODELS)}
    require(encode(means)==proof['means']and len(means)==len(proof['per_mean']),'incomplete shared-mean closure')
    plain=set();corrected=set();nodes=0
    for mean,detail in zip(sorted(means),proof['per_mean']):
        require(encode({mean})[0]==detail['mean'],'wrong mean')
        delta=[{op('add',m,x,-mean)for m in MODELS}for x in xs]
        square=[{op('mul',m,d,d)for d,m in itertools.product(values,MODELS)}for values in delta]
        q,t=all_sums(square),all_sums(delta)
        corrections={op('div',b,op('mul',a,d,d),R(3))for d,a,b in itertools.product(t['7'],MODELS,MODELS)}
        unadjusted={op('div',m,v,divisor)for v,m in itertools.product(q['7'],MODELS)}
        adjusted={op('div',b,op('add',a,v,-c),divisor)for v,c,a,b in itertools.product(q['7'],corrections,MODELS,MODELS)}
        require(list(map(encode,delta))==detail['deltas']and list(map(encode,square))==detail['squares'],'incomplete correlated leaf closure')
        require(encoded(q)==detail['square_sum']and encoded(t)==detail['delta_sum'],'incomplete sum-tree closure')
        require(encode(corrections)==detail['corrections']and encode(unadjusted)==detail['plain_roots']and
                encode(adjusted)==detail['corrected_roots'],'incomplete reconstruction closure')
        plain.update(unadjusted);corrected.update(adjusted)
        nodes+=sum(len(v)for v in q.values())+sum(len(v)for v in t.values())+len(corrections)+len(unadjusted)+len(adjusted)
    require(encode(plain)==proof['plain_roots']and encode(corrected)==proof['corrected_roots'],'incomplete root union')
    target=binary64(row['expected']);span=B.pullback(Span(target,target))
    require(span.json()==proof['target_span'],'incorrect exact target rounding interval')
    # Comparing exact rational roots to the open/closed rounding interval also
    # verifies every final publication, including tie endpoint parity.
    reachable=lambda values:any(span.contains(v)for v in values)
    require(proof['plain_reachable']==reachable(plain)and proof['corrected_reachable']==reachable(corrected),'invalid reachability claim')
    return dict(plain_unreachable=int(not reachable(plain)),corrected_unreachable=int(not reachable(corrected)),
                mean_values=len(means),plain_root_values=len(plain),corrected_root_values=len(corrected)),nodes


def check(input_path,directory,output):
    start=time.monotonic();directory=Path(directory);source=Path(__file__)
    sources={str(p):digest(p)for p in (source,source.with_name('exact.py'),source.with_name('check.py'))}
    rows=list(records(input_path));proofs=list(records(directory/'closures.jsonl'))
    require(len(rows)==len(proofs),'missing original obligation')
    counts=Counter();states=0
    for row,proof in zip(rows,proofs):
        delta,nodes=verify(row,proof);counts.update(delta);counts['observations']+=1;states+=nodes
    summary=json.loads((directory/'summary.json').read_text())
    require(summary['complete']is True and summary['counts']==counts,'incorrect closure summary')
    require(all(digest(p)==h for p,h in sources.items()),'checker changed during proof')  # ubs:ignore — public integrity hashes, not secrets
    receipt=dict(complete=True,counts=counts,full_width_value_states=states,seconds=time.monotonic()-start,
                 arithmetic_cache=op.cache_info()._asdict(),source_sha256=sources,
                 input_sha256={str(p):digest(p)for p in (Path(input_path),*sorted(directory.iterdir()))if p.is_file()},
                 limits='Every specified three-input arithmetic tree/site assignment is covered. Corrected closures are dependency-free supersets. '
                        'No target is used to produce an intermediate, and reachability is not a reusable program or a product gain.')
    with open(output,'x')as out:json.dump(receipt,out,indent=2);out.write('\n')
    print(json.dumps(receipt))


def online_verify(row,proof):
    require(proof['raw']==row,'changed online source observation')
    def flatten(v):return [y for x in v for y in flatten(x)]if isinstance(v,list)else[v]
    bits=flatten(row['args']);require(len(bits)==3 and all(isinstance(v,str)and v.startswith('0x')for v in bits),'online numeric scope')
    xs=list(map(binary64,bits));require(all(v==0 or two(-100)<=abs(v)<=two(100)for v in xs),'online magnitude bound')
    require(row['function']in('VAR','VAR.S','VAR.P','VARP','VARA','VARPA'),'online function scope')
    divisor=R(3)if row['function']in('VAR.P','VARP','VARPA')else R(2)
    def rounded(kind,a,b):return{op(kind,m,a,b)for m in MODELS}
    def advance(family,state,x,k):
        a,q=state;outputs=set()
        if family=='welford':
            differences=rounded('add',x,-a)
            for delta in differences:
                for quotient in rounded('div',delta,R(k)):
                    for mean in rounded('add',a,quotient):
                        for residual in rounded('add',x,-mean):
                            for product in rounded('mul',delta,residual):
                                outputs.update((mean,value)for value in rounded('add',q,product))
        else:
            for total in rounded('add',a,x):
                for multiple in rounded('mul',R(k),x):
                    for residual in rounded('add',multiple,-total):
                        for square in rounded('mul',residual,residual):
                            for contribution in rounded('div',square,R(k*(k-1))):
                                outputs.update((total,value)for value in rounded('add',q,contribution))
        return outputs
    permutations=list(itertools.permutations(range(3)))
    require([p['indices']for p in proof['permutations']]==[list(p)for p in permutations],'incomplete online permutations')
    unions={k:set()for k in('welford','youngs_cramer')};counts=Counter()
    for indices,record in zip(permutations,proof['permutations']):
        for family in unions:
            states={(xs[indices[0]],R(0))};stages=[]
            for k in(2,3):
                following=set()
                for state in states:following.update(advance(family,state,xs[indices[k-1]],k))
                states=following;stages.append([[encode({a})[0],encode({b})[0]]for a,b in sorted(states)])
                counts[family+'_state_pairs']+=len(states)
            require(record[family]['stages']==stages,'incomplete correlated online boundary states')
            values={v for _,q in states for v in rounded('div',q,divisor)}
            require(record[family]['roots']==encode(values),'incomplete online per-order roots')
            unions[family].update(values)
    span=B.pullback(Span(binary64(row['expected']),binary64(row['expected'])))
    require(span.json()==proof['target_span'],'changed online target preimage')
    require(set(proof['roots'])==set(unions)and set(proof['reachable'])==set(unions),'changed online family coverage')
    for family,values in unions.items():
        require(proof['roots'][family]==encode(values),'incomplete online root union')
        reachable=any(span.contains(v)for v in values);require(proof['reachable'][family]==reachable,'incorrect online exclusion')
        counts[family+'_roots']+=len(values);counts[family+'_unreachable']+=not reachable
    counts['observations']+=1
    return counts


def online_check(input_path,directory,output):
    start=time.monotonic();directory=Path(directory);source=Path(__file__)
    hashes={str(p):digest(p)for p in(source,source.with_name('exact.py'),source.with_name('check.py'))}
    rows=list(records(input_path));proofs=list(records(directory/'closures.jsonl'));require(rows and len(rows)==len(proofs),'missing online obligations')
    counts=Counter()
    for row,proof in zip(rows,proofs):counts.update(online_verify(row,proof))
    summary=json.loads((directory/'summary.json').read_text());require(summary['complete']is True and summary['counts']==counts,'incorrect online summary')
    require(all(digest(p)==h for p,h in hashes.items()),'checker changed during online closure')  # ubs:ignore — public integrity hashes, not secrets
    receipt=dict(complete=True,counts=counts,seconds=time.monotonic()-start,source_sha256=hashes,
                 input_sha256={str(p):digest(p)for p in(Path(input_path),*sorted(directory.iterdir()))if p.is_file()},
                 limits='Independent exact boundary-state reproduction for the specified two unrolled online algorithms. All per-site choices and6 index permutations covered; no target-derived state, complete kernel identification or product gain.')
    with open(output,'x')as out:json.dump(receipt,out,indent=2);out.write('\n')
    print(json.dumps(receipt))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='online':online_check(*sys.argv[2:])
    else:check(*sys.argv[1:])
