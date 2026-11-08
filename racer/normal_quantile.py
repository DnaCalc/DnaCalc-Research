"""AS241 expression orders from its published coefficient tables.

Tier 2 reference: CPython 3.13 Modules/_statisticsmodule.c (PSF licence).
The ISA graphs and searches are original MIT research, not a port of its code.
"""
import itertools
import random
import numpy as np
from .ets import *
from .x87_programs import Builder,Native
from .ets_proof import freeze

# Reuse the existing inverse-program registry instead of keeping a second table.
import re
_TABLE_SOURCE=Path(__file__).with_name('src')/'normal.rs'
_table_text=_TABLE_SOURCE.read_text().split('const AS241:[[f64;8];6]=[',1)[1].split('];',1)[0]
_tables=[[float(v.strip().replace('_','')) for v in row.split(',') if v.strip()]
         for row in re.findall(r'\[([^\[\]]+)\]',_table_text)]
if len(_tables)!=6 or any(len(row)!=8 for row in _tables):raise ValueError('AS241 registry shape changed')
CENTRAL_N,CENTRAL_D,TAIL_N,TAIL_D,FAR_N,FAR_D=_tables

def tail(p):
    r=math.sqrt(-math.log(p if p<.5 else 1-p))
    nn,dd=(TAIL_N,TAIL_D)if r<=5 else(FAR_N,FAR_D);r-=1.6 if r<=5 else 5.
    n=nn[0];d=dd[0]
    for c in nn[1:]:n=n*r+c
    for c in dd[1:]:d=d*r+c
    return math.copysign(n/d,p-.5)

def tail_proof():
    name='step6-normal-tail';rng=random.Random(396777);points=[]
    for i in range(160):
        p=rng.uniform(.926,.999999)if i<80 else math.exp(rng.uniform(-24.5,-3.))
        points.append(p)
    pred=[bits(tail(p))for p in points];freeze(name,dict(points=points,predictions=pred,extra_source=Path(__file__).read_text(),method='Published AS241 moderate tail, binary64 operations, system log and square root.'))
    c,s=service();answers=s.ask_rows('NORM.S.INV',[[c.number(p)]for p in points]);fail=[dict(index=i,prediction=p,answer=a)for i,(p,a)in enumerate(zip(pred,answers))if p!=a.get('bits')];write(OUT/(name+'-answers.json'),dict(points=points,answers=answers,receipts=s.receipts));write(OUT/(name+'-proof.json'),dict(failures=fail,selections=len(points)));print('AS241 tail refutations',len(fail),fail[:5])

def upper_mapping():
    centers=[.936,.937,.93749,.9375,.93751,.938,.94,.95,.975,.99,.9999,math.nextafter(1.,0.)];points=[]
    for p in centers:
        raw=int(bits(p),16);points.extend(value(f'{k:016x}')for k in range(raw-32,min(raw+32,0x3ff0000000000000)))
    pred=[bits(tail(p))for p in points];freeze('step6-normal-upper-mapping',dict(points=points,predictions=pred,extra_source=Path(__file__).read_text(),method='Single-bit regions across the candidate upper-tail transition and exponent edges.'))
    c,s=service();ans=s.ask_rows('NORM.S.INV',[[c.number(p)]for p in points],chunk_size=1000);write(OUT/'step6-normal-upper-mapping-answers.json',dict(points=points,answers=ans,receipts=s.receipts));print([(p0,sum(pred[j]==ans[j].get('bits')for j,p in enumerate(points)if int(bits(p0),16)-32<=int(bits(p),16)<min(int(bits(p0),16)+32,0x3ff0000000000000)))for p0 in centers])
TAGS=['q','square','r','nproduct','nsum','dproduct','dsum','product','quotient']
def central_graph():
    g=Builder(1)
    def op(kind,a,b,tag):
        k=g.op(kind,a,b,store=False);g.nodes[k].update(site=TAGS.index(tag),tag=tag);return k
    def c(v):return g.constant(bits(v))
    q=op('sub',0,c(.5),'q');r=op('sub',c(.180625),op('mul',q,q,'square'),'r')
    def horner(coefficients,prefix):
        v=c(coefficients[0])
        for a in coefficients[1:]:v=op('add',op('mul',v,r,prefix+'product'),c(a),prefix+'sum')
        return v
    n=horner(CENTRAL_N,'n');d=horner(CENTRAL_D,'d')
    before=op('div',op('mul',n,q,'product'),d,'quotient');after=op('mul',q,op('div',n,d,'quotient'),'product');g.sites=len(TAGS);return g.graph(before,after)

def discovery():
    d=read(OUT/'step6-quantile-coupling-normal.json');rows=[];expected=[]
    for conf,a in zip(d['confidences'],d['answers']):
        p=(1+conf)/2
        if p>.925:continue
        rows.append([int(bits(p),16)]);expected.append(int(a['bits'],16))
    native=Native(OUT/'x87-ets-loop.so');g=central_graph();profiles=[]
    for cw,mask in itertools.product([0x123f,0x133f],range(1<<len(TAGS))):
        out=native.grid(g,np.array(rows,dtype='uint64'),cw,mask)
        for output in [0,1]:
            if np.array_equal(out[:,output],expected):profiles.append(dict(cw=cw,mask=mask,output=output,stores=[t for j,t in enumerate(TAGS)if mask&(1<<j)]))
    write(OUT/'step6-normal-central-profiles.json',dict(profiles=profiles,producer_sha256=sha(__file__)));print('central profiles',len(profiles),profiles[:4])

def proof():
    name='step6-normal-central';rng=random.Random(396621);points=[rng.uniform(.076,.924)for _ in range(128)];profiles=read(OUT/'step6-normal-central-profiles.json')['profiles'];native=Native(OUT/'x87-ets-loop.so');g=central_graph();rows=np.array([[int(bits(p),16)]for p in points],dtype='uint64');pred={}
    for i,p in enumerate(profiles):pred[str(i)]=[f'{int(v):016x}'for v in native.grid(g,rows,p['cw'],p['mask'])[:,p['output']]]
    freeze(name,dict(points=points,predictions=pred,profiles=profiles,extra_source=Path(__file__).read_text()))
    c,s=service();answers=s.ask_rows('NORM.S.INV',[[c.number(p)]for p in points]);actual=[a.get('bits')for a in answers];survive=[i for i,v in pred.items()if v==actual];write(OUT/(name+'-answers.json'),dict(points=points,answers=answers,receipts=s.receipts));write(OUT/(name+'-proof.json'),dict(survivors=survive,counts={i:sum(a==b for a,b in zip(v,actual))for i,v in pred.items()}));print('fresh central programs',survive)


def store_refutations(mapping, output, library):
    """Replay the original fixed AS241 store family on any fresh mapping.

    Every rejected profile keeps its first exact contradiction. Selecting the
    shortest/simple profile is allowed only among complete survivors.
    """
    d=read(mapping)
    selected=[i for i,p in enumerate(d['points'])if .075<p<.925]
    if not selected:raise ValueError('no central observations')
    rows=np.array([[int(bits(d['points'][i]),16)]for i in selected],dtype='uint64')
    actual=np.array([int(d['answers'][i]['bits'],16)for i in selected],dtype='uint64')
    native=Native(library);g=central_graph();refutations=[];survivors=[]
    for cw,mask in itertools.product([0x123f,0x133f],range(1<<len(TAGS))):
        out=native.grid(g,rows,cw,mask)
        for branch in [0,1]:
            bad=np.flatnonzero(out[:,branch]!=actual)
            profile=dict(cw=cw,mask=mask,output=branch)
            if len(bad):
                j=int(bad[0]);refutations.append(dict(**profile,case=selected[j],p=d['points'][selected[j]],prediction=f'{int(out[j,branch]):016x}',answer=f'{int(actual[j]):016x}'))
            else:survivors.append(profile)
    write(output,dict(observations=len(rows),profiles=len(refutations)+len(survivors),survivors=survivors,
                     first_contradictions=refutations,source_sha256=sha(__file__),mapping_sha256=sha(mapping)))
    return survivors

if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='proof':proof()
    elif len(sys.argv)>1 and sys.argv[1]=='tail':tail_proof()
    elif len(sys.argv)>1 and sys.argv[1]=='upper':upper_mapping()
    else:discovery()
