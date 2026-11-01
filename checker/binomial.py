"""Independent exact verifier for binomial initial-mass recurrence exclusions."""
from collections import Counter,defaultdict
from fractions import Fraction as R
import itertools
import json
from pathlib import Path
import sys
import time

from exact import BINARY64 as B,EXTENDED64 as E,Span,binary64,operation,two
from check import records,require,digest

MODELS=('strict','stored_pc64','resident_pc64','stored_pc53','resident_pc53')


def pair(value):return [str(value.numerator),str(value.denominator)]


def cohorts(rows):
    grouped=defaultdict(dict);counts=Counter()
    for index,row in enumerate(rows):
        a=row['args'];counts['decoded_rows']+=1
        if len(a)!=4 or row['result_index']is not None or not all(isinstance(x,str)and x.startswith('0x')for x in a[:3]):
            counts['outside_numeric_input_shape']+=1;continue
        flag=a[3]
        if not(flag is False or isinstance(flag,str)and flag.startswith('0x')and binary64(flag)==0):
            counts['cdf_or_non_numeric_flag']+=1;continue
        k,n,p=map(binary64,a[:3])
        if k.denominator!=1 or n.denominator!=1 or not(0<=k<=min(n,64))or not(1<=n<=2**53)or not(0<p<1):
            counts['outside_integer_count_or_probability_scope']+=1;continue
        if not row['expected'].startswith('0x')or abs(binary64(row['expected']))<two(-1022):
            counts['outside_normal_nonzero_publication']+=1;continue
        context=json.dumps([row['function'],a[1],a[2],row['reading'],row['source']],sort_keys=True,separators=(',',':'))
        if int(k)in grouped[context]:
            require(rows[grouped[context][int(k)]]['expected']==row['expected'],'same-context conflict')
            counts['duplicate_same_count_same_answer']+=1
        else:grouped[context][int(k)]=index
    result=[]
    for context,group in sorted(grouped.items()):
        if 0 not in group or len(group)<2:
            counts['groups_without_initial_or_later_mass']+=1;continue
        result.append(dict(context=context,initial=group[0],targets=[group[k]for k in sorted(group)if k>0]))
    counts['cohorts']=len(result);counts['target_observations']=sum(len(g['targets'])for g in result)
    return result,counts


def step(profile,value,numerator,denominator,p,q):
    m=profile['models'][1:];form=profile['shape']
    def op(site,kind,a,b):return operation(a,b,kind,m[site])
    if form==0:
        a=op(0,'mul',value,numerator);b=op(1,'div',a,denominator);c=op(2,'mul',b,p);return op(3,'div',c,q)
    if form==1:
        a=op(0,'mul',value,p);b=op(1,'div',a,q);c=op(2,'mul',b,numerator);return op(3,'div',c,denominator)
    if form==2:
        a=op(0,'mul',numerator,p);b=op(1,'mul',value,a);c=op(2,'mul',denominator,q);return op(3,'div',b,c)
    if form==3:
        a=op(0,'mul',numerator,p);b=op(1,'mul',denominator,q);c=op(2,'div',a,b);return op(3,'mul',value,c)
    if form==4:
        a=op(0,'div',value,denominator);b=op(1,'mul',a,numerator);c=op(2,'div',p,q);return op(3,'mul',b,c)
    a=op(0,'div',numerator,denominator);b=op(1,'mul',value,a);c=op(2,'div',p,q);return op(3,'mul',b,c)


def constraint(rows,g,p):
    initial=rows[g['initial']];n=binary64(initial['args'][1]);prob=binary64(initial['args'][2])
    q=operation(R(1),-prob,'add',p['models'][0]);require(q>0,'nonpositive complement')
    span=B.pullback(Span(binary64(initial['expected']),binary64(initial['expected'])))
    first,last=E.coordinates(span);lower,upper=E.value(first),E.value(last)
    original=(lower,upper);count=0
    for order,index in enumerate(g['targets']):
        target=rows[index];k=int(binary64(target['args'][0]))
        for i in range(count,k):
            lower=step(p,lower,n-i,R(i+1),prob,q);upper=step(p,upper,n-i,R(i+1),prob,q)
            require(lower<=upper,'lost monotonicity')
        count=k;published=[B.round(lower),B.round(upper)];y=binary64(target['expected'])
        if not(published[0]<=y<=published[1]):
            return dict(target_order=order,target_index=index,count=k,initial_bounds=list(map(pair,original)),
                        complement=pair(q),root_bounds=list(map(pair,(lower,upper))),
                        published_bounds=list(map(pair,published)),target=pair(y))
    return None


def check(input_path,directory,output):
    start=time.monotonic();directory=Path(directory);source=Path(__file__)
    hashes={str(p):digest(p)for p in (source,source.with_name('exact.py'),source.with_name('check.py'))}
    rows=list(records(input_path));gs,scope=cohorts(rows)
    require(gs==json.loads((directory/'groups.json').read_text()),'changed cohorts')
    require(scope==json.loads((directory/'scope.json').read_text()),'changed scope accounting')
    ps=[dict(id=i,shape=s,models=list(m))for i,(s,m)in enumerate(itertools.product(range(6),itertools.product(MODELS,repeat=5)))]
    require(ps==json.loads((directory/'profiles.json').read_text())and len(ps)==18750,'incomplete grammar')
    summary=json.loads((directory/'summary.json').read_text());require(summary['complete']is True,'partial run')
    killed=set();executions=0;counts=Counter();per_group=Counter()
    for record in records(directory/'first-kills.jsonl'):
        identity,order=record['candidate_id'],record['group_order']
        require(identity not in killed and 0<=identity<len(ps)and 0<=order<len(gs),'duplicate or unknown profile/group')
        for previous in range(order):
            require(constraint(rows,gs[previous],ps[identity])is None,'not first contradiction');executions+=1
        proof=constraint(rows,gs[order],ps[identity]);executions+=1
        require(proof is not None and proof==record['proof'],'invalid exact initial-mass box certificate')
        counts['certified_first_kills']+=1;per_group[str(order)]+=1;killed.add(identity)
    alive=json.loads((directory/'necessary-survivors.json').read_text())
    require(len(set(alive))==len(alive)and killed.isdisjoint(alive)and killed|set(alive)==set(range(len(ps))),'unsearched profile')
    for identity in alive:
        for g in gs:require(constraint(rows,g,ps[identity])is None,'invalid necessary survivor');executions+=1
    require(summary['profiles']==len(ps)and summary['cohorts']==len(gs)and summary['executions']==executions and
            summary['rejected']==len(killed)and summary['necessary_survivors']==len(alive)and summary['scope']==scope and
            summary['first_kill_groups']==per_group,'incorrect result counts')
    require(all(digest(p)==h for p,h in hashes.items()),'checker changed during run')  # ubs:ignore — public integrity hashes, not secrets
    receipt=dict(complete=True,profiles=len(ps),cohorts=len(gs),rejected=len(killed),necessary_survivors=len(alive),
                 executions=executions,counts=counts,seconds=time.monotonic()-start,source_sha256=hashes,
                 input_sha256={str(p):digest(p)for p in (Path(input_path),*sorted(directory.iterdir()))if p.is_file()},
                 claim='Independent rational reproduction of monotone bounds quantifying over every Ext80 initial mass publishing to observed P(0). No initial mass is selected or fitted.',
                 limits='Specified finite consumers, normal publications and original same-context transport assumptions only. No complete probability producer or kernel is inferred.')
    with open(output,'x')as out:json.dump(receipt,out,indent=2);out.write('\n')
    print(json.dumps(receipt))


def cdf_cohorts(rows):
    grouped=defaultdict(lambda:dict(pmf={},cdf={}));counts=Counter()
    shapes={'BINOM.DIST':4,'POISSON.DIST':3,'HYPGEOM.DIST':5,'NEGBINOM.DIST':4}
    for index,row in enumerate(rows):
        counts['decoded_rows']+=1;a=row['args'];name=row['function']
        if name not in shapes or len(a)!=shapes[name]or row['result_index']is not None or not all(isinstance(v,str)and v.startswith('0x')for v in a[:-1]):
            counts['outside_numeric_shape']+=1;continue
        count=binary64(a[0])
        if not(0<=count<=64)or count.denominator!=1:counts['outside_count_scope']+=1;continue
        flag=a[-1]
        if isinstance(flag,bool):kind='cdf'if flag else'pmf'
        elif isinstance(flag,str)and flag.startswith('0x'):kind='pmf'if binary64(flag)==0 else'cdf'
        else:counts['flag_outside_scope']+=1;continue
        if not row['expected'].startswith('0x')or not two(-1022)<=binary64(row['expected'])<=1:
            counts['outside_normal_positive_probability']+=1;continue
        key=json.dumps([name,a[1:-1],row['reading'],row['source']],sort_keys=True,separators=(',',':'))
        bucket=grouped[key][kind];k=int(count)
        if k in bucket:
            require(rows[bucket[k]]['expected']==row['expected'],'same-context CDF/PMF conflict');counts['same_answer_duplicates']+=1
        else:bucket[k]=index
    result=[]
    for key,group in sorted(grouped.items()):
        for count,target in sorted(group['cdf'].items()):
            if any(k not in group['pmf']for k in range(count+1)):
                counts['cdf_without_complete_positive_pmf_prefix']+=1;continue
            result.append(dict(context=key,target=target,masses=[group['pmf'][k]for k in range(count+1)]))
    counts['cohorts']=len(result);return result,counts


def cdf_bounds(rows,group,profile):
    seeds=[]
    for index in group['masses']:
        value=binary64(rows[index]['expected']);span=B.pullback(Span(value,value))
        lo,hi=E.coordinates(span);seeds.append((E.value(lo),E.value(hi)))
    def add(a,b):return operation(a,b,'add',profile['model'])
    def total(values):
        if profile['order']==2:
            if len(values)==1:return values[0]
            middle=len(values)//2;return add(total(values[:middle]),total(values[middle:]))
        result=R(0)
        for value in values if profile['order']==0 else values[::-1]:result=add(result,value)
        return result
    low=total([lo for lo,_ in seeds]);high=total([hi for _,hi in seeds])
    published=[B.round(low),B.round(high)]
    if profile['clip']:published=[min(R(1),v)for v in published]
    target=binary64(rows[group['target']]['expected'])
    require(low<=high and published[0]<=published[1],'nonmonotone CDF bounds')
    return dict(seed_bounds=[[pair(lo),pair(hi)]for lo,hi in seeds],root_bounds=list(map(pair,(low,high))),
                published_bounds=list(map(pair,published)),target=pair(target),excluded=not(published[0]<=target<=published[1]))


def cdf_check(input_path,directory,output):
    start=time.monotonic();directory=Path(directory);source=Path(__file__)
    hashes={str(p):digest(p)for p in (source,source.with_name('exact.py'),source.with_name('check.py'))}
    rows=list(records(input_path));groups,scope=cdf_cohorts(rows)
    require(groups==json.loads((directory/'cdf-groups.json').read_text())and scope==json.loads((directory/'cdf-scope.json').read_text()),'changed complete-prefix scope')
    profiles=[dict(id=i,order=o,model=m,clip=c)for i,(o,m,c)in enumerate(itertools.product(range(3),MODELS,(False,True)))]
    require(profiles==json.loads((directory/'cdf-profiles.json').read_text())and len(profiles)==30,'changed CDF grammar')
    summary=json.loads((directory/'summary.json').read_text());require(summary['complete']is True,'partial CDF study')
    results={};total_executions=0;exclusions=0
    for name in ('BINOM.DIST','POISSON.DIST','HYPGEOM.DIST','NEGBINOM.DIST'):
        cohort=[(i,g)for i,g in enumerate(groups)if rows[g['target']]['function']==name]
        killed=set();executions=0
        for record in records(directory/f'{name}-cdf-first-kills.jsonl'):
            identity,order=record['candidate_id'],record['order']
            require(identity not in killed and 0<=identity<30 and 0<=order<len(cohort)and record['group_index']==cohort[order][0],'duplicate/unknown CDF proof')
            for previous in range(order):
                require(not cdf_bounds(rows,cohort[previous][1],profiles[identity])['excluded'],'not first CDF contradiction');executions+=1
            proof=cdf_bounds(rows,cohort[order][1],profiles[identity]);executions+=1
            require(proof['excluded']and proof==record['proof'],'invalid complete-prefix bounds');killed.add(identity)
        alive=json.loads((directory/f'{name}-cdf-necessary-survivors.json').read_text())
        require(len(set(alive))==len(alive)and killed.isdisjoint(alive)and killed|set(alive)==set(range(30)),'unsearched CDF profile')
        for identity in alive:
            for _,group in cohort:require(not cdf_bounds(rows,group,profiles[identity])['excluded'],'invalid CDF necessary survivor');executions+=1
        result=dict(profiles=30,cohorts=len(cohort),executions=executions,rejected=len(killed),necessary_survivors=len(alive))
        require(summary['counts'][name]==result,'incorrect CDF result summary');results[name]=result
        total_executions+=executions;exclusions+=len(killed)
    require(all(digest(p)==h for p,h in hashes.items()),'checker changed during CDF run')  # ubs:ignore — public integrity hashes, not secrets
    receipt=dict(complete=True,functions=results,scope=scope,exact_seed_boxes=total_executions,certified_exclusions=exclusions,
                 seconds=time.monotonic()-start,source_sha256=hashes,
                 input_sha256={str(p):digest(p)for p in (Path(input_path),*sorted(directory.iterdir()))if p.is_file()},
                 limits='Necessary bounds under original same-context/shared-PMF assumptions only. No latent probability producer, CDF algorithm or product gain is identified.')
    with open(output,'x')as out:json.dump(receipt,out,indent=2);out.write('\n')
    print(json.dumps(receipt))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='cdf':cdf_check(*sys.argv[2:])
    else:check(*sys.argv[1:])
