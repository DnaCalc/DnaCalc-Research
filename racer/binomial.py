"""Necessary binomial recurrence constraints; never a fitted initial mass.

The defining ratio P(k+1)/P(k)=(n-k)*p/((k+1)*(1-p)) supplies six fixed
consumer trees. Every initial mass is an arbitrary Ext80 value whose publication
matches the independently observed P(0). Monotonicity propagates exact bounds.
"""
import argparse
from collections import defaultdict,Counter
from fractions import Fraction as Q
import itertools
import json
from pathlib import Path
import sys
import time

from bessel import MODELS,operation,seed_bounds,pair,digest,require
from intervals import B64,number,pow2


def false_flag(arg):
    return arg is False or (isinstance(arg,str)and arg.startswith('0x')and number(arg)==0)


def cohorts(rows):
    grouped=defaultdict(dict);counts=Counter()
    for index,row in enumerate(rows):
        a=row['args'];counts['decoded_rows']+=1
        if len(a)!=4 or row['result_index']is not None or not all(isinstance(x,str)and x.startswith('0x')for x in a[:3]):
            counts['outside_numeric_input_shape']+=1;continue
        if not false_flag(a[3]):counts['cdf_or_non_numeric_flag']+=1;continue
        k,n,p=map(number,a[:3])
        if k!=int(k)or n!=int(n)or not(0<=k<=min(n,64))or not(1<=n<=2**53)or not(0<p<1):
            counts['outside_integer_count_or_probability_scope']+=1;continue
        if not row['expected'].startswith('0x')or abs(Q(number(row['expected'])))<pow2(-1022):
            counts['outside_normal_nonzero_publication']+=1;continue
        key=json.dumps([row['function'],a[1],a[2],row['reading'],row['source']],sort_keys=True,separators=(',',':'))
        k=int(k)
        if k in grouped[key]:
            require(rows[grouped[key][k]]['expected']==row['expected'], 'conflicting same-context reading')
            counts['duplicate_same_count_same_answer']+=1
        else:grouped[key][k]=index
    result=[]
    for key,values in sorted(grouped.items()):
        if 0 not in values or len(values)<2:
            counts['groups_without_initial_or_later_mass']+=1;continue
        result.append(dict(context=key,initial=values[0],targets=[values[k]for k in sorted(values)if k>0]))
    counts['cohorts']=len(result)
    counts['target_observations']=sum(len(g['targets'])for g in result)
    return result,counts


def profiles():
    for identity,(shape,models)in enumerate(itertools.product(range(6),itertools.product(MODELS,repeat=5))):
        yield dict(id=identity,shape=shape,models=list(models))


def update(v,a,b,p,q,profile):
    m=profile['models'][1:];shape=profile['shape']
    def op(kind,x,y,site):return operation(kind,m[site],x,y)
    if shape==0:
        v=op('mul',v,a,0);v=op('div',v,b,1);v=op('mul',v,p,2);return op('div',v,q,3)
    if shape==1:
        v=op('mul',v,p,0);v=op('div',v,q,1);v=op('mul',v,a,2);return op('div',v,b,3)
    if shape==2:
        numerator=op('mul',v,op('mul',a,p,0),1)
        denominator=op('mul',b,q,2)
        return op('div',numerator,denominator,3)
    if shape==3:
        numerator=op('mul',a,p,0);denominator=op('mul',b,q,1)
        ratio=op('div',numerator,denominator,2);return op('mul',v,ratio,3)
    if shape==4:
        v=op('div',v,b,0);v=op('mul',v,a,1);ratio=op('div',p,q,2)
        return op('mul',v,ratio,3)
    ratio=op('div',a,b,0);v=op('mul',v,ratio,1);ratio=op('div',p,q,2)
    return op('mul',v,ratio,3)


def constrain(rows,group,profile):
    initial=rows[group['initial']]
    _,nb,pb=initial['args'][:3]
    n,p=Q(number(nb)),Q(number(pb))
    q=operation('add',profile['models'][0],Q(1),-p)
    require(q>0, 'complement outside positive scope')
    lo,hi=seed_bounds(Q(number(initial['expected'])))
    original=(lo,hi);k=0
    for order,index in enumerate(group['targets']):
        target=rows[index];requested=int(number(target['args'][0]))
        while k<requested:
            lo=update(lo,n-k,Q(k+1),p,q,profile)
            hi=update(hi,n-k,Q(k+1),p,q,profile)
            require(lo<=hi, 'nonmonotone recurrence')
            k+=1
        low,high=B64.rn(lo),B64.rn(hi)
        require(low is not None and high is not None, 'nonfinite publication bounds')
        y=Q(number(target['expected']))
        if not low<=y<=high:
            return dict(target_order=order,target_index=index,count=k,
                        initial_bounds=list(map(pair,original)),complement=pair(q),
                        root_bounds=list(map(pair,(lo,hi))),published_bounds=list(map(pair,(low,high))),target=pair(y))
    return None


def run(input_path,destination):
    start=time.monotonic();dest=Path(destination);dest.mkdir()
    rows=[json.loads(line)for line in Path(input_path).read_text().splitlines()]
    gs,scope=cohorts(rows);ps=list(profiles());require(len(ps)==18750 and gs, 'incomplete grammar or empty cohort')
    for name,data in [('groups.json',gs),('profiles.json',ps),('scope.json',scope)]:
        with(dest/name).open('x')as out:json.dump(data,out,indent=2);out.write('\n')
    survivors=[];executions=0;rejections=Counter()
    source=digest(Path(__file__))
    with(dest/'first-kills.jsonl').open('x')as kills:
        for profile in ps:
            for order,group in enumerate(gs):
                proof=constrain(rows,group,profile);executions+=1
                if proof is not None:
                    kills.write(json.dumps(dict(candidate_id=profile['id'],group_order=order,proof=proof))+'\n')
                    rejections[order]+=1
                    break
            else:survivors.append(profile['id'])
    with(dest/'necessary-survivors.json').open('x')as out:json.dump(survivors,out);out.write('\n')
    require(digest(Path(__file__))==source, 'source changed during run')  # ubs:ignore — public integrity hashes, not secrets
    result=dict(complete=True,profiles=len(ps),cohorts=len(gs),executions=executions,
                rejected=len(ps)-len(survivors),necessary_survivors=len(survivors),seconds=time.monotonic()-start,
                scope=scope,first_kill_groups=dict(rejections),input_sha256={str(input_path):digest(input_path)},source_sha256=source,
                claim='All positive Ext80 initial masses whose binary64 publication equals recorded P(0) are bounded through each fixed recurrence. '
                      'No initial mass or helper producer is selected. Necessary survivors need not share a feasible initial mass across targets.',
                limits='Same original source/environment, exact integer n/k, 0<p<1, k<=64 and normal nonzero PMF publications. '
                       'No CDF, endpoint algorithm, other recurrence tree, additional store/branch or Excel kernel is identified.')
    with(dest/'summary.json').open('x')as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(result))


def cdf_cohorts(rows):
    """Only complete numeric PMF prefixes from the original same context."""
    grouped=defaultdict(lambda:dict(pmf={},cdf={}))
    counts=Counter()
    arities={'BINOM.DIST':4,'POISSON.DIST':3,'HYPGEOM.DIST':5,'NEGBINOM.DIST':4}
    for index,row in enumerate(rows):
        counts['decoded_rows']+=1
        a=row['args'];name=row['function']
        if name not in arities or len(a)!=arities[name]or row['result_index']is not None or not all(isinstance(x,str)and x.startswith('0x')for x in a[:-1]):
            counts['outside_numeric_shape']+=1;continue
        k=number(a[0])
        if not(0<=k<=64)or k!=int(k):counts['outside_count_scope']+=1;continue
        flag=a[-1]
        if isinstance(flag,bool):kind='cdf'if flag else'pmf'
        elif isinstance(flag,str)and flag.startswith('0x'):kind='pmf'if number(flag)==0 else'cdf'
        else:counts['flag_outside_scope']+=1;continue
        if not row['expected'].startswith('0x')or not pow2(-1022)<=Q(number(row['expected']))<=1:
            counts['outside_normal_positive_probability']+=1;continue
        key=json.dumps([name,a[1:-1],row['reading'],row['source']],sort_keys=True,separators=(',',':'))
        bucket=grouped[key][kind];k=int(k)
        if k in bucket:
            require(rows[bucket[k]]['expected']==row['expected'],'conflicting same-context CDF/PMF reading')
            counts['same_answer_duplicates']+=1
        else:bucket[k]=index
    result=[]
    for context,group in sorted(grouped.items()):
        for k,target in sorted(group['cdf'].items()):
            if not all(i in group['pmf']for i in range(k+1)):
                counts['cdf_without_complete_positive_pmf_prefix']+=1;continue
            result.append(dict(context=context,target=target,masses=[group['pmf'][i]for i in range(k+1)]))
    counts['cohorts']=len(result)
    return result,counts


def cdf_profiles():
    for identity,(order,model,clip)in enumerate(itertools.product(range(3),MODELS,(False,True))):
        yield dict(id=identity,order=order,model=model,clip=clip)


def cdf_bounds(rows,group,profile):
    seeds=[seed_bounds(Q(number(rows[i]['expected'])))for i in group['masses']]
    def total(values):
        if profile['order']==2:
            if len(values)==1:return values[0]
            middle=len(values)//2
            return operation('add',profile['model'],total(values[:middle]),total(values[middle:]))
        value=Q(0)
        for x in values if profile['order']==0 else values[::-1]:value=operation('add',profile['model'],value,x)
        return value
    roots=[total([s[side]for s in seeds])for side in (0,1)]
    published=[B64.rn(x)for x in roots]
    if profile['clip']:published=[min(Q(1),x)for x in published]
    target=Q(number(rows[group['target']]['expected']))
    require(roots[0]<=roots[1]and published[0]<=published[1],'nonmonotone CDF sum')
    return dict(seed_bounds=[[pair(lo),pair(hi)]for lo,hi in seeds],root_bounds=list(map(pair,roots)),
                published_bounds=list(map(pair,published)),target=pair(target),
                excluded=not(published[0]<=target<=published[1]))


def cdf_run(input_path,destination):
    start=time.monotonic();dest=Path(destination);dest.mkdir();source=digest(Path(__file__))
    rows=[json.loads(line)for line in Path(input_path).read_text().splitlines()]
    groups,scope=cdf_cohorts(rows);profiles=list(cdf_profiles());require(groups,'empty complete-prefix scope')
    for name,data in [('cdf-groups.json',groups),('cdf-profiles.json',profiles),('cdf-scope.json',scope)]:
        with(dest/name).open('x')as out:json.dump(data,out,indent=2);out.write('\n')
    results={}
    for function in ('BINOM.DIST','POISSON.DIST','HYPGEOM.DIST','NEGBINOM.DIST'):
        cohort=[(i,g)for i,g in enumerate(groups)if rows[g['target']]['function']==function]
        killed=[];alive=[];executions=0
        with(dest/f'{function}-cdf-first-kills.jsonl').open('x')as out:
            for p in profiles:
                for order,(index,g)in enumerate(cohort):
                    proof=cdf_bounds(rows,g,p);executions+=1
                    if proof['excluded']:
                        out.write(json.dumps(dict(candidate_id=p['id'],order=order,group_index=index,proof=proof))+'\n');killed.append(p['id']);break
                else:alive.append(p['id'])
        with(dest/f'{function}-cdf-necessary-survivors.json').open('x')as out:json.dump(alive,out);out.write('\n')
        results[function]=dict(profiles=30,cohorts=len(cohort),executions=executions,rejected=len(killed),necessary_survivors=len(alive))
    require(digest(Path(__file__))==source, 'source changed during CDF run')  # ubs:ignore — public integrity hashes, not secrets
    summary=dict(complete=True,counts=results,scope=scope,seconds=time.monotonic()-start,source_sha256=source,input_sha256={str(input_path):digest(input_path)},
                 claim='Necessary monotone bounds for CDF sums of unpublished Ext80 PMFs constrained by complete independently recorded PMF prefixes.',
                 limits='Original same-context transport and shared PMF-publication assumptions; positive normal masses, k<=64, three sum orders, five models, optional final clipping at1. '
                        'Necessary survivors are not complete programs or consistent latent producers. No CDF kernel or product gain follows.')
    with(dest/'summary.json').open('x')as out:json.dump(summary,out,indent=2);out.write('\n')
    print(json.dumps(summary))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='cdf':
        cdf_run(*sys.argv[2:]);sys.exit(0)
    parser=argparse.ArgumentParser();parser.add_argument('input');parser.add_argument('destination')
    args=parser.parse_args();run(args.input,args.destination)
