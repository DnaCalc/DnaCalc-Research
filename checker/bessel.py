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


if __name__=='__main__':check(*sys.argv[1:])
