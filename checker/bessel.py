"""Independent exact proof of necessary Bessel order-two recurrence bounds."""
from collections import defaultdict,Counter
from fractions import Fraction as R
import itertools
import json
from pathlib import Path
import sys
import time

from exact import BINARY64 as B, EXTENDED64 as E, Span, binary64, operation, two
from check import records, require, digest

MODELS=('strict','stored_pc64','resident_pc64','stored_pc53','resident_pc53')


def decode_pair(pair):
    require(len(pair)==2 and all(isinstance(x,str)for x in pair),'rational certificate encoding')
    return R(*map(int,pair))


def groups(rows):
    result=defaultdict(dict)
    for index,r in enumerate(rows):
        args=r['args']
        if len(args)!=2 or r['result_index']is not None or not all(isinstance(x,str)and x.startswith('0x')for x in args):continue
        x,n=map(binary64,args)
        if x<=0 or n not in (0,1,2)or not r['expected'].startswith('0x'):continue
        if abs(binary64(r['expected']))<two(-1022):continue
        key=json.dumps([r['function'],args[0],r['reading'],r['source']],sort_keys=True,separators=(',',':'))
        if int(n)in result[key]:
            require(rows[result[key][int(n)]]['expected']==r['expected'],'conflicting source/environment')
        else:result[key][int(n)]=index
    return [dict(context=k,indices=[v[n]for n in (0,1,2)])for k,v in sorted(result.items())if set(v)=={0,1,2}]


def seed_bounds(target):
    span=B.pullback(Span(target,target))
    lo,hi=E.coordinates(span)
    return E.value(lo),E.value(hi)


def calculate(profile,function,x,seed0,seed1):
    m=profile['models'];form=profile['shape']
    if form==0:
        t=operation(R(2),seed1,'mul',m[0]);t=operation(t,x,'div',m[1])
    elif form==1:
        t=operation(seed1,x,'div',m[0]);t=operation(R(2),t,'mul',m[1])
    elif form==2:
        t=operation(R(2),x,'div',m[0]);t=operation(t,seed1,'mul',m[1])
    else:
        d=operation(x,R(2),'div',m[0]);t=operation(seed1,d,'div',m[1])
    if function in ('BESSELJ','BESSELY'):return operation(t,-seed0,'add',m[2])
    return operation(seed0,t if function=='BESSELK'else-t,'add',m[2])


def proof(rows,g,p,certificate=None):
    zero,one,target=[rows[i]for i in g['indices']]
    function=zero['function'];x=binary64(zero['args'][0]);y=binary64(target['expected'])
    z=seed_bounds(binary64(zero['expected']));o=seed_bounds(binary64(one['expected']))
    if function in ('BESSELJ','BESSELY'):corners=((z[1],o[0]),(z[0],o[1]))
    elif function=='BESSELI':corners=((z[0],o[1]),(z[1],o[0]))
    else:corners=((z[0],o[0]),(z[1],o[1]))
    roots=[calculate(p,function,x,*pair)for pair in corners]
    published=[B.round(root)for root in roots]
    require(roots[0]<=roots[1]and published[0]<=published[1],'monotone orientation failed')
    excluded=not(published[0]<=y<=published[1])
    if certificate is not None:
        for field,expected in [('seed0',z),('seed1',o),('root_bounds',roots),('published_bounds',published)]:
            require(list(map(decode_pair,certificate[field]))==list(expected),'incorrect '+field)
        require([[decode_pair(a),decode_pair(b)]for a,b in certificate['corner_seeds']]==[list(c)for c in corners],'incorrect extremal seed pairing')
        require(decode_pair(certificate['target'])==y and certificate['excluded']==excluded,'changed target or exclusion')
    return excluded


def check(input_path,directory,output):
    started=time.monotonic();directory=Path(directory)
    source=Path(__file__);hashes={str(p):digest(p)for p in (source,source.with_name('exact.py'),source.with_name('check.py'))}
    rows=list(records(input_path));gs=groups(rows)
    require(gs==json.loads((directory/'groups.json').read_text()),'changed group selection')
    ps=[dict(id=i,shape=s,models=list(m))for i,(s,m)in enumerate(itertools.product(range(4),itertools.product(MODELS,repeat=3)))]
    require(ps==json.loads((directory/'profiles.json').read_text())and len(ps)==500,'changed grammar')
    summary=json.loads((directory/'summary.json').read_text());require(summary['complete']is True,'partial search')
    result={};counts=Counter()
    for function in ('BESSELJ','BESSELY','BESSELI','BESSELK'):
        cohort=[(i,g)for i,g in enumerate(gs)if rows[g['indices'][0]]['function']==function]
        killed=set();executions=0
        for row in records(directory/f'{function}-first-kills.jsonl'):
            identity=row['candidate_id'];order=row['order']
            require(identity not in killed and 0<=identity<len(ps),'duplicate or unknown profile')
            require(0<=order<len(cohort)and row['group_index']==cohort[order][0],'changed first-kill group')
            for previous in range(order):
                require(not proof(rows,cohort[previous][1],ps[identity]),'not first rejection');executions+=1
            require(proof(rows,cohort[order][1],ps[identity],row['proof']),'invalid rejection');executions+=1
            killed.add(identity);counts['independently_checked_seed_box_exclusions']+=1
        alive=json.loads((directory/f'{function}-necessary-survivors.json').read_text())
        require(len(set(alive))==len(alive)and killed.isdisjoint(alive)and killed|set(alive)==set(range(500)),'incomplete profile accounting')
        for identity in alive:
            for _,g in cohort:
                require(not proof(rows,g,ps[identity]),'excluded necessary survivor');executions+=1
        actual=dict(profiles=500,groups=len(cohort),executions=executions,rejected=len(killed),necessary_survivors=len(alive))
        require(summary['counts'][function]==actual,'incorrect summary')
        result[function]=actual
        counts['independently_checked_seed_boxes']+=executions
        counts['full_width_corner_results']+=2*executions
    require(all(digest(p)==h for p,h in hashes.items()),'checker changed during run')  # ubs:ignore — public integrity hashes, not secrets
    receipt=dict(complete=True,counts=counts,functions=result,seconds=time.monotonic()-started,
                 source_sha256=hashes,input_sha256={str(p):digest(p)for p in (Path(input_path),*sorted(directory.iterdir()))if p.is_file()},
                 claim='Exact monotone corner proof over all Ext80 seed pairs publishing to recorded orders zero and one. Necessary constraints only; no seed pair is fitted or selected.',
                 limits='Positive x, normal nonzero publications, fixed specified consumers and same original source/environment. A necessary survivor is not an executable complete program.')
    with open(output,'x')as out:json.dump(receipt,out,indent=2);out.write('\n')
    print(json.dumps(receipt))


def check_domains(input_path,directory,output):
    from reference import log_bounds,LIB
    directory=Path(directory)
    frozen={str(p):digest(p)for p in (Path(input_path),*directory.iterdir(),Path(__file__),
            Path(__file__).with_name('reference.py'),Path(__file__).with_name('exact.py'))if p.is_file()}
    def atan_bounds(denominator,terms):
        total=sum((R((-1)**i,(2*i+1)*denominator**(2*i+1))for i in range(terms)),R(0))
        next_term=R((-1)**terms,(2*terms+1)*denominator**(2*terms+1))
        return min(total,total+next_term),max(total,total+next_term)
    a,b=atan_bounds(5,80),atan_bounds(239,40)
    pi=(16*a[0]-4*b[1],16*a[1]-4*b[0])
    logmax=log_bounds(0x7fefffffffffffff)
    constants={}
    for name,span in [('pi',pi),('logmax',logmax),('ribesl_exparg',(R(709),R(709))),('cody_i0_xmax',(R('713.986'),R('713.986')))]:
        constants[name]={}
        for kind,lattice in [('binary64',B),('extended64',E)]:
            lo,hi=map(lattice.round,span);require(lo==hi,'ambiguous independent constant load')
            constants[name][kind]=[str(lo.numerator),str(lo.denominator)]
    require(json.loads((directory/'constants.json').read_text())==constants,'incorrect domain constants')
    rows=[];scope=Counter()
    for index,raw in enumerate(records(input_path)):
        if raw['function']not in ('BESSELI','BESSELJ','BESSELY'):scope['other_function']+=1;continue
        try:
            require(raw['result_index']is None and len(raw['args'])==2,'shape')
            x,n=map(binary64,raw['args']);require(n.denominator==1 and 0<=n<=100,'order')
            require(raw['expected'].startswith('0x')or raw['expected']=='error:Num','result')
            require(abs(x)>=500 if raw['function']=='BESSELI'else abs(x)>=two(24)and n<=4 and(raw['function']!='BESSELY'or x>0),'large_argument_scope')
        except(ValueError,TypeError,OverflowError)as error:scope[str(error)]+=1;continue
        rows.append(dict(trace_index=index,raw=raw));scope['scoped']+=1
    require(rows and any(r['raw']['function']=='BESSELI'for r in rows)and any(r['raw']['function']=='BESSELJ'for r in rows),'empty domain group')
    require(json.loads((directory/'scoped.json').read_text())==rows,'changed domain scope')
    ps=[]
    for family,forms in [('I',range(3)),('JY',range(5))]:
        for form,model,load,inclusive in itertools.product(forms,('strict',)if family=='I'else MODELS,('binary64','extended64'),(False,True)):
            ps.append(dict(id=len(ps),family=family,form=form,model=model,load=load,inclusive=inclusive))
    require(json.loads((directory/'profiles.json').read_text())==ps,'changed domain grammar')
    def rejects(p,raw):
        x,n=map(binary64,raw['args']);x=abs(x);n=int(n)
        if p['family']=='I':
            value=decode_pair(constants[('ribesl_exparg','logmax','cody_i0_xmax')[p['form']]][p['load']])
            return x>=value if p['inclusive']else x>value
        pi_value=decode_pair(constants['pi'][p['load']])
        orders=[n]if p['form']in(1,3)else[0,1]
        angles=[]
        for order in orders:
            multiple=2*order+(-1 if p['form']in(3,4)else 1)
            product=operation(R(multiple),pi_value,'mul',p['model'])
            offset=operation(product,R(4),'div',p['model'])
            angles.append(abs(operation(x,-offset,'add',p['model'])))
        if p['form']==0:angles=[x]
        return any(v>=two(27)if p['inclusive']else v>two(27)for v in angles)
    expected=[]
    for p in ps:
        hits=[];early=[];observed=0
        for row in rows:
            raw=row['raw']
            if(raw['function']=='BESSELI')!=(p['family']=='I'):continue
            refusal=rejects(p,raw);numeric=raw['expected']!='error:Num';observed+=1
            if refusal and numeric:early.append(row['trace_index'])
            if refusal==numeric:hits.append(row['trace_index'])
        require(observed,'empty domain profile')
        expected.append(dict(candidate_id=p['id'],observations=observed,rejects_recorded_numeric=len(early),
            sole_explanation_mismatches=len(hits),first_early_contradiction=early[0]if early else None,
            first_sole_contradiction=hits[0]if hits else None))
    require(json.loads((directory/'results.json').read_text())==expected,'incorrect domain contradiction')
    summary=json.loads((directory/'summary.json').read_text())
    require(summary['complete']and summary['programs']==len(ps)and summary['observations']==len(rows)and summary['scope']==scope and
            summary['sole_explanation_survivors']==sum(r['sole_explanation_mismatches']==0 for r in expected)and
            summary['early_guard_compatible']==sum(r['rejects_recorded_numeric']==0 for r in expected),'wrong domain summary')
    require(all(digest(p)==h for p,h in frozen.items()),'domain evidence/source changed')
    result=dict(complete=True,programs=len(ps),observations=len(rows),sole_explanation_survivors=summary['sole_explanation_survivors'],
        early_guard_compatible=summary['early_guard_compatible'],input_sha256=frozen,mpfr_version=LIB.mpfr_get_version().decode(),
        limits='Independent rational arithmetic, Machin brackets and directed MPFR logmax loads. Early-guard impossibilities and incomplete sole explanations are counted separately. No numeric Bessel body, fitted cutoff or domain repair is identified.')
    with open(output,'x')as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(result))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='domains':check_domains(*sys.argv[2:])
    else:check(*sys.argv[1:])
