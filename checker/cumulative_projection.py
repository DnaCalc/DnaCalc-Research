"""Independent exact checks of cumulative growth and unknown-principal/payment projections."""
from collections import defaultdict, Counter
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import sys
import time

from check import records, require, digest, graph_trace, tagged_value, write_result
from cumulative import coordinates, node, MODELS
from exact import BINARY64 as B, EXTENDED64 as E, Span, binary64, inverse, operation

SOURCE=digest(Path(__file__))
CUMULATIVE_SOURCE=digest(Path(__file__).with_name('cumulative.py'))


def group(row):
    return json.dumps([row['args'][0],row['args'][1],row['source'],row['reading']],sort_keys=True,separators=(',',':'))


def growth_shape(power,model,degree):
    one=node('const','strict',0x3ff0000000000000)
    q=node('add',model,one,node('arg','strict',0))
    value=one
    if power==0:
        for _ in range(degree):value=node('mul',model,value,q)
    else:
        base=q
        while degree:
            if degree%2:value=node('mul',model,value,base)
            degree//=2
            if degree:base=node('mul',model,base,base)
    return value


def verify_growth(rows,path):
    available=defaultdict(list)
    for i,row in enumerate(rows):available[group(row)].append(i)
    helpers={};nodes_checked=0
    for h in records(path):
        raw=rows[h['row_index']];key=group(raw)
        require(json.loads(h['group'])==json.loads(key),'wrong growth observation group')
        args,n,*_=coordinates(raw)
        require(h['row_index']==available[key][0]and h['periods']==n,'growth representative/bound differs')
        require(list(map(binary64,h['input_bits']))==args,'growth inputs differ')
        model=h['growth_model'];power=h['power']
        require(model in MODELS and power in range(2),'undeclared growth profile')
        graph=h['graph'];values=graph_trace(graph,args)
        require(len(values)==len(h['values']),'growth value count differs')
        shapes=[]
        for j,(instruction,actual,value)in enumerate(zip(graph['nodes'],h['values'],values)):
            kind,operands=next(iter(instruction['op'].items()));m=instruction['model']
            if kind=='arg':
                require(operands==0 and m=='strict','undeclared growth argument')
                shape=node(kind,m,operands)
            elif kind=='const':
                require(operands['bits']in (0,0x3ff0000000000000,0xbff0000000000000)and m=='strict','undeclared growth constant')
                shape=node(kind,m,operands['bits'])
            else:
                require(kind in ('add','mul')and all(0<=k<j for k in operands)and m==model,'undeclared growth operation')
                shape=node(kind,m,*(shapes[k]for k in operands))
            require(tagged_value(actual)==value,'growth arithmetic differs')
            tag='ext80-le:'if kind not in ('arg','const')and m.startswith('resident_')else'binary64:'
            require(actual.startswith(tag),'growth intermediate representation differs')
            shapes.append(shape);nodes_checked+=1
        require(len(h['outputs'])==len(h['factors'])==n+1,'incomplete power outputs')
        factors=[]
        for k,(index,actual)in enumerate(zip(h['outputs'],h['factors'])):
            require(0<=index<len(values)and shapes[index]==growth_shape(power,model,k),'growth expression differs')
            require(actual==h['values'][index],'growth factor publication changed')
            factors.append(values[index])
        hk=(key,power,model);require(hk not in helpers,'duplicate growth profile/group')
        helpers[hk]=factors
    require(set(helpers)==set(itertools.product(available,range(2),MODELS)),'incomplete growth grammar/group coverage')
    return helpers,nodes_checked


def check(trace,growth,directory,output,companions=None):
    started=time.monotonic();directory=Path(directory)
    rows=list(records(trace));require(len(rows)==540,'discovery rows missing')
    for row in rows:coordinates(row)
    helpers,nodes_checked=verify_growth(rows,growth)
    payment_map={}
    if companions:
        for p in records(companions):
            a=p['args'];require(p['function']=='PMT'and p['result_index']is None and len(a)==5 and binary64(a[3])==0,'wrong PMT companion surface')
            k=tuple(a[:3]+[a[4],json.dumps(p['reading'],sort_keys=True)])
            require(k not in payment_map,'duplicate companion')
            require('G6-cumprinc' in p['source'],'companion batch provenance missing')
            payment_map[k]=p
    groups=defaultdict(list)
    for i,r in enumerate(rows):
        a=r['args']
        if (companions or binary64(a[5])==0)and a[3]==a[4]:
            groups[group(r)+':'+a[2]+':'+a[5]].append(i)
    profiles=json.loads((directory/'profiles.json').read_text())
    declared=[dict(id=i,power=p,growth_model=g,consumer_model=c,
        family='shared_payment'if companions else'shared_initial_principal')
        for i,(p,g,c)in enumerate(itertools.product(range(2),MODELS,MODELS))]
    require(profiles==declared,'projection grammar differs')
    actual=list(records(directory/'constraints.jsonl'));position=0;survivors=[];counts=Counter()
    for profile in profiles:
        killed=False
        for loan,indices in sorted(groups.items()):
            require(position<len(actual),'missing group certificate')
            record=actual[position];position+=1
            require(record['profile_id']==profile['id']and record['loan']==loan,'group/profile order differs')
            payment=None;current=None
            if companions:
                raw=rows[indices[0]];a=raw['args']
                payment=payment_map[tuple(a[:3]+[a[5],json.dumps(raw['reading'],sort_keys=True)])]
                y=binary64(payment['expected']);current=B.pullback(Span(y,y))
                require(record['payment_observation']==payment,'PMT binding differs')
            require(record['payment_observation']==payment,'unexpected payment constraint')
            used=[]
            for i in indices:
                r=rows[i];a=r['args'];period=int(binary64(a[3]));y=binary64(r['expected'])
                degree=int(binary64(a[1]))-period+1 if companions else period-1
                factor=helpers[group(r),profile['power'],profile['growth_model']][degree]
                published=B.pullback(Span(y,y))
                required=(published if companions and binary64(a[5])==1 and period==1 else
                    inverse(published,factor,'div'if companions else'mul',True,profile['consumer_model']))
                current=required if not used and not companions else current.meet(required)if current else None
                item=record['preimages'][len(used)]
                require(item['row_index']==i and item['period']==period,'singleton binding differs')
                require(item['required']==(required.json()if required else None),'inverse rounding endpoints differ')
                require(item['intersection']==(current.json()if current else None),'intersection endpoints differ')
                used.append(i);counts['singleton_constraints_checked']+=1
                if current is None or E.coordinates(current)is None:break
            require(len(used)==len(record['preimages']),'incomplete/redundant certificate constraints')
            possible=current is not None and E.coordinates(current)is not None
            reason='nonempty'if possible else'empty_real_intersection'if current is None else'empty_ext80_lattice'
            require(record['feasible']==possible and record['reason']==reason,'feasibility conclusion differs')
            if not possible:
                killed=True;counts[reason]+=1;break
        if not killed:survivors.append(profile['id'])
    require(position==len(actual),'extra group certificates')
    require(survivors==json.loads((directory/'survivors.json').read_text()),'survivor IDs differ')
    require(digest(Path(__file__))==SOURCE and digest(Path(__file__).with_name('cumulative.py'))==CUMULATIVE_SOURCE,'checker changed during execution')  # ubs:ignore — public source hashes, no secret comparison
    write_result(output,dict(complete=True,profiles=50,survivors=survivors,original_available_rows=540,
        scoped_singletons=sum(map(len,groups.values())),loans=len(groups),growth_nodes_checked=nodes_checked,
        **counts,seconds=time.monotonic()-started,projection_checker_sha256=SOURCE,
        inputs={str(p):digest(p)for p in ([Path(trace),Path(growth),directory/'profiles.json',directory/'constraints.jsonl']+([Path(companions)]if companions else[]))},
        claim='Independent exact growth replay and rounding-preimage intersections, preserving endpoints and original typed loan/environment identities. All declared profiles checked; no coefficient selected or fitted.',
        limitations='Only the declared singleton consumers and conditional growth programs. A shared latent is an existential constraint, never an implemented component. Companion source/build/bitness are preserved; historical CPU, transport and settings are unavailable. No full-function identification or fresh-validation claim.'))


if __name__=='__main__':check(*sys.argv[1:])
