"""Exact necessary recurrence constraints with unpublished common Ext80 seeds.

Original mathematics, MIT. DLMF 10.6.1 and 10.29.1. No seed, coefficient,
threshold or answer is fitted or selected. Bounds quantify over all Ext80 seeds
whose binary64 publications agree with independently recorded orders 0 and 1.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path
import time

from intervals import B64, E64, E53, Interval, number, pow2

MODELS = ('strict','stored_pc64','resident_pc64','stored_pc53','resident_pc53')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pair(value):
    return [str(value.numerator), str(value.denominator)]


def context(row):
    return json.dumps([row['function'],row['args'][0],row['reading'],row['source']],
                      sort_keys=True,separators=(',',':'))


def groups(rows):
    grouped = defaultdict(dict)
    for i,row in enumerate(rows):
        args = row['args']
        if len(args) != 2 or row['result_index'] is not None:
            continue
        if not all(isinstance(a,str) and a.startswith('0x') for a in args):
            continue
        x, n = map(number,args)
        if x <= 0 or n not in (0,1,2) or not row['expected'].startswith('0x'):
            continue
        y = Q(number(row['expected']))
        if abs(y) < pow2(-1022):
            continue
        key = context(row)
        if int(n) in grouped[key]:
            require(rows[grouped[key][int(n)]]['expected'] == row['expected'], 'conflicting same-context observation')
        else:
            grouped[key][int(n)] = i
    return [dict(context=k, indices=[v[n]for n in range(3)])
            for k,v in sorted(grouped.items()) if set(v)=={0,1,2}]


def seed_bounds(target):
    span = B64.preimage(Interval(target,target))
    lo,hi = E64.ceil(span.lo),E64.floor(span.hi)
    if lo==span.lo and not span.lc:lo=E64.next(lo)
    if hi==span.hi and not span.hc:hi=E64.prev(hi)
    require(lo is not None and hi is not None and lo<=hi, 'empty finite seed bounds')
    return lo,hi


@lru_cache(maxsize=131072)
def operation(kind, model, a, b):
    if model=='strict':
        a,b=B64.rn(a),B64.rn(b)
    exact={'mul':lambda:a*b,'div':lambda:a/b,'add':lambda:a+b}[kind]()
    fmt=B64 if model=='strict'else E53 if model.endswith('pc53')else E64
    value=fmt.rn(exact)
    require(value is not None, 'nonfinite recurrence outside scope')
    if model.startswith('stored'):value=B64.rn(value)
    require(value is not None, 'nonfinite stored recurrence outside scope')
    return value


def recurrence(x, zero, one, function, shape, models):
    m0,m1,m2=models
    if shape==0:
        term=operation('div',m1,operation('mul',m0,Q(2),one),x)
    elif shape==1:
        term=operation('mul',m1,Q(2),operation('div',m0,one,x))
    elif shape==2:
        term=operation('mul',m1,operation('div',m0,Q(2),x),one)
    else:
        term=operation('div',m1,one,operation('div',m0,x,Q(2)))
    if function in ('BESSELJ','BESSELY'):
        return operation('add',m2,term,-zero)
    if function=='BESSELI':
        return operation('add',m2,zero,-term)
    require(function=='BESSELK', 'unknown Bessel consumer')
    return operation('add',m2,zero,term)


def profiles():
    for identity,(shape,models) in enumerate(itertools.product(range(4),itertools.product(MODELS,repeat=3))):
        yield dict(id=identity,shape=shape,models=list(models))


def bounds(rows,group,profile):
    a,b,target=[rows[i]for i in group['indices']]
    function=a['function'];x=Q(number(a['args'][0]))
    zero=seed_bounds(Q(number(a['expected'])))
    one=seed_bounds(Q(number(b['expected'])))
    if function in ('BESSELJ','BESSELY'):
        corners=((zero[1],one[0]),(zero[0],one[1]))
    elif function=='BESSELI':
        corners=((zero[0],one[1]),(zero[1],one[0]))
    else:
        corners=((zero[0],one[0]),(zero[1],one[1]))
    roots=[recurrence(x,z,o,function,profile['shape'],profile['models'])for z,o in corners]
    published=list(map(B64.rn,roots))
    require(roots[0]<=roots[1]and published[0]<=published[1], 'nonmonotone corner orientation')
    expected=Q(number(target['expected']))
    return dict(seed0=list(map(pair,zero)),seed1=list(map(pair,one)),
                corner_seeds=[[pair(z),pair(o)]for z,o in corners],
                root_bounds=list(map(pair,roots)),published_bounds=list(map(pair,published)),
                target=pair(expected),excluded=not(published[0]<=expected<=published[1]))


def run(input_path,destination):
    started=time.monotonic();dest=Path(destination);dest.mkdir()
    rows=[json.loads(line)for line in Path(input_path).read_text().splitlines()]
    gs=groups(rows);ps=list(profiles());counts={};source=digest(Path(__file__))
    for name,obj in [('groups.json',gs),('profiles.json',ps)]:
        with(dest/name).open('x')as out:json.dump(obj,out,indent=2);out.write('\n')
    for function in ('BESSELJ','BESSELY','BESSELI','BESSELK'):
        scoped=[(i,g)for i,g in enumerate(gs)if rows[g['indices'][0]]['function']==function]
        require(scoped, 'empty recurrence cohort')
        survivors=[];executions=0
        with(dest/f'{function}-first-kills.jsonl').open('x')as kills:
            for p in ps:
                for order,(index,g)in enumerate(scoped):
                    result=bounds(rows,g,p);executions+=1
                    if result['excluded']:
                        kills.write(json.dumps(dict(candidate_id=p['id'],order=order,group_index=index,proof=result))+'\n')
                        break
                else:survivors.append(p['id'])
        with(dest/f'{function}-necessary-survivors.json').open('x')as out:json.dump(survivors,out);out.write('\n')
        counts[function]=dict(profiles=len(ps),groups=len(scoped),executions=executions,
                              rejected=len(ps)-len(survivors),necessary_survivors=len(survivors))
    require(digest(Path(__file__))==source, 'source changed during run')  # ubs:ignore — public integrity hashes, not secrets
    result=dict(complete=True,counts=counts,seconds=time.monotonic()-started,
                input_sha256={str(input_path):digest(input_path)},source_sha256=source,
                claim='Necessary constraints on order-two recurrence consumers with arbitrary common Ext80 order-zero/one seeds. '
                      'Monotone corner bounds cover all seed pairs; a surviving bound need not contain any actual seed pair or producer.',
                limits='Same original source/environment and positive input; nonzero normal publications at orders 0,1,2 only. '
                       'No seed generator, higher-order branch, coefficient, cutoff or complete Excel program is identified.')
    with(dest/'summary.json').open('x')as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(result))


def domain_run(input_path, destination):
    """Test fixed published guards as sole explanations, without inventing a seed body."""
    from constants import atan, function_bounds
    dest=Path(destination);dest.mkdir()
    a,b=atan(Q(1,5),64),atan(Q(1,239),32)
    pi=(16*a[0]-4*b[1],16*a[1]-4*b[0])
    maximum=Q(2**1024-2**971)
    logmax=function_bounds(maximum,maximum,'ln')
    def load(span,lattice):
        low,high=map(lattice.rn,span);require(low==high,'ambiguous universal constant load');return low
    constants={name:{kind:pair(load(span,lattice)) for kind,lattice in [('binary64',B64),('extended64',E64)]}
        for name,span in [('pi',pi),('logmax',logmax),('ribesl_exparg',(Q(709),Q(709))),('cody_i0_xmax',(Q('713.986'),Q('713.986')))]}
    rows=[];scope=Counter()
    for index,raw in enumerate(map(json.loads,Path(input_path).read_text().splitlines())):
        if raw['function']not in ('BESSELI','BESSELJ','BESSELY'):
            scope['other_function']+=1;continue
        try:
            require(raw['result_index']is None and len(raw['args'])==2,'shape')
            x,n=map(lambda v:Q(number(v)),raw['args']);require(n.denominator==1 and 0<=n<=100,'order')
            require(raw['expected'].startswith('0x')or raw['expected']=='error:Num','result')
            if raw['function']=='BESSELI':require(abs(x)>=500,'large_argument_scope')
            else:require(abs(x)>=2**24 and n<=4 and (raw['function']!='BESSELY'or x>0),'large_argument_scope')
        except (ValueError,TypeError,OverflowError)as error:
            scope[str(error)]+=1;continue
        rows.append(dict(trace_index=index,raw=raw));scope['scoped']+=1
    require(rows and any(r['raw']['function']=='BESSELI'for r in rows)and any(r['raw']['function']=='BESSELJ'for r in rows),'empty domain group')
    profiles=[]
    for family,forms in [('I',range(3)),('JY',range(5))]:
        models=('strict',)if family=='I'else MODELS
        for form,model,load_kind,inclusive in itertools.product(forms,models,('binary64','extended64'),(False,True)):
            profiles.append(dict(id=len(profiles),family=family,form=form,model=model,load=load_kind,inclusive=inclusive))
    (dest/'constants.json').write_text(json.dumps(constants,indent=2)+'\n')
    (dest/'scoped.json').write_text(json.dumps(rows)+'\n');(dest/'profiles.json').write_text(json.dumps(profiles)+'\n')
    def rejects(p,raw):
        x,n=map(lambda v:Q(number(v)),raw['args']);x=abs(x);n=int(n)
        if p['family']=='I':
            bound=Q(*map(int,constants[('ribesl_exparg','logmax','cody_i0_xmax')[p['form']]][p['load']]))
            return x>=bound if p['inclusive']else x>bound
        pi_value=Q(*map(int,constants['pi'][p['load']]))
        def phase(order,sine):
            multiplier=2*order+(-1 if sine else 1)
            numerator=operation('mul',p['model'],Q(multiplier),pi_value)
            correction=operation('div',p['model'],numerator,Q(4))
            return abs(operation('add',p['model'],x,-correction))
        if p['form']==0:arguments=[x]
        elif p['form']in(1,3):arguments=[phase(n,p['form']==3)]
        else:arguments=[phase(order,p['form']==4)for order in (0,1)]
        return any(v>=2**27 if p['inclusive']else v>2**27 for v in arguments)
    results=[]
    for p in profiles:
        observed=early=sole=0;first_early=first_sole=None
        for row in rows:
            raw=row['raw']
            if (raw['function']=='BESSELI')!=(p['family']=='I'):continue
            rejection=rejects(p,raw);error=raw['expected']=='error:Num';observed+=1
            if rejection and not error:
                early+=1
                if first_early is None:first_early=row['trace_index']
            if rejection!=error:
                sole+=1
                if first_sole is None:first_sole=row['trace_index']
        require(observed,'empty guard profile scope')
        results.append(dict(candidate_id=p['id'],observations=observed,
            rejects_recorded_numeric=early,sole_explanation_mismatches=sole,
            first_early_contradiction=first_early,first_sole_contradiction=first_sole))
    (dest/'results.json').write_text(json.dumps(results,indent=2)+'\n')
    result=dict(complete=True,programs=len(profiles),observations=len(rows),scope=scope,
        sole_explanation_survivors=sum(r['sole_explanation_mismatches']==0 for r in results),
        early_guard_compatible=sum(r['rejects_recorded_numeric']==0 for r in results),
        source_sha256=digest(Path(__file__)),input_sha256={str(input_path):digest(input_path)},
        sources=['https://www.netlib.org/specfun/ribesl','https://www.netlib.org/specfun/i0',
                 'A&S 9.2 asymptotic phases; Machin identity; IEEE754',
                 'DnaCalc real_policy: measured circular trig limit 2^27'],
        limits='A guard rejecting a recorded numeric result is impossible as an early guard in this scope. A guard admitting a recorded NUM fails as its sole explanation, but may still precede a body that fails. No Bessel seed producer, fitted cutoff, new Excel answer or numerical repair is inferred.')
    (dest/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='domains':domain_run(*sys.argv[2:]);raise SystemExit(0)
    parser=argparse.ArgumentParser();parser.add_argument('input');parser.add_argument('destination')
    args=parser.parse_args();run(args.input,args.destination)
