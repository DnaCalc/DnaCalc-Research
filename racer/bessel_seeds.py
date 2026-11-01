"""Conservative all-site rounding boxes for the published Cody J0 rational.

Each coefficient independently admits its stated decimal's RN53 and RN64 load.
Every arithmetic node independently admits all five existing precision/store
models. Interval hulls discard dependencies, so an exclusion is conservative.
"""
from collections import Counter
from fractions import Fraction as Q
import itertools
import json
from pathlib import Path
import sys
import time

from bessel import MODELS,operation,require,pair,digest
from intervals import B64,E64,Interval,number


class Boxes:
    def __init__(self,tables,x):
        self.tables=tables;self.nodes=[];self.input=self.append('input',(x,x),absolute=True)
    def append(self,kind,bounds,**fields):
        self.nodes.append(dict(op=kind,bounds=list(map(pair,bounds)),**fields));return len(self.nodes)-1
    def literal(self,value):return self.append('literal',(Q(value),Q(value)),decimal=str(value))
    def load(self,name,index=None):
        raw=self.tables['constants'][name]if index is None else self.tables['tables'][name][index]
        q=Q(raw['decimal']);values=(B64.rn(q),E64.rn(q))
        return self.append('coefficient',(min(values),max(values)),name=name,index=index)
    def bounds(self,index):return [Q(*map(int,p))for p in self.nodes[index]['bounds']]
    def op(self,kind,a,b):
        xs,ys=self.bounds(a),self.bounds(b);values=[]
        if kind=='div':require(ys[0]*ys[1]>0,'rational denominator crosses zero')
        for x,y,model in itertools.product(xs,ys,MODELS):
            values.append(operation('add'if kind=='sub'else kind,model,x,-y if kind=='sub'else y))
        return self.append(kind,(min(values),max(values)),a=a,b=b)
    def horner(self,z,coefficients):
        result=coefficients[0]
        for coefficient in coefficients[1:]:result=self.op('add',self.op('mul',result,z),coefficient)
        return result


def bound(raw,tables,factor,association):
    x=abs(Q(number(raw['args'][0])));b=Boxes(tables,x);argument=b.input
    z=b.op('mul',argument,argument)
    numerator=b.horner(z,[b.load('PJ0',i)for i in (4,5,6,0,1,2,3)])
    denominator=b.horner(z,[b.literal(1)]+[b.load('QJ0',i)for i in (4,0,1,2,3)])
    root=b.load('XJ0')
    if factor==0:
        split=b.op('div',b.load('XJ01'),b.literal(256))
        left=b.op('sub',b.op('sub',argument,split),b.load('XJ02'))
        product=b.op('mul',left,b.op('add',argument,root))
    elif factor==1:
        product=b.op('mul',b.op('sub',argument,root),b.op('add',argument,root))
    else:product=b.op('sub',z,b.op('mul',root,root))
    if association==0:result=b.op('div',b.op('mul',product,numerator),denominator)
    else:result=b.op('mul',product,b.op('div',numerator,denominator))
    lo,hi=b.bounds(result);published=(B64.rn(lo),B64.rn(hi));target=Q(number(raw['expected']))
    return dict(factor=factor,association=association,nodes=b.nodes,output=result,
                published_bounds=list(map(pair,published)),target=pair(target),
                excluded=not(published[0]<=target<=published[1]))


def run(input_path,tables_path,destination):
    start=time.monotonic();dest=Path(destination);dest.mkdir();source=digest(Path(__file__))
    original=[json.loads(line)for line in Path(input_path).read_text().splitlines()]
    tables=json.loads(Path(tables_path).read_text());rows=[];excluded=[]
    for i,row in enumerate(original):
        a=row['args']
        if row['function']=='BESSELJ'and row['result_index']is None and len(a)==2 and all(isinstance(x,str)and x.startswith('0x')for x in a)and number(a[1])==0 and Q(1,4)<=abs(Q(number(a[0])))<=4 and row['expected'].startswith('0x'):
            rows.append(dict(trace_index=i,raw=row))
        else:excluded.append(i)
    require(rows,'empty Cody J0 scope');counts=Counter();alive=set(range(6))
    with(dest/'bounds.jsonl').open('x')as out:
        for row in rows:
            proofs=[]
            for factor,association in itertools.product(range(3),range(2)):
                proof=bound(row['raw'],tables,factor,association);proofs.append(proof)
                counts['coefficient_load_nodes']+=sum(n['op']=='coefficient'for n in proof['nodes'])
                counts['arithmetic_box_nodes']+=sum(n['op']in ('add','sub','mul','div')for n in proof['nodes'])
                counts['excluded_boxes']+=proof['excluded'];counts['boxes']+=1
                if proof['excluded']:alive.discard(2*factor+association)
            out.write(json.dumps(dict(**row,proofs=proofs))+'\n')
    require(digest(Path(__file__))==source,'source changed during seed study')  # ubs:ignore — public integrity hashes, not secrets
    summary=dict(complete=True,recorded=len(original),scoped=len(rows),outside_scope=len(excluded),counts=counts,
                 necessary_survivors=sorted(alive),seconds=time.monotonic()-start,source_sha256=source,
                 input_sha256={str(p):digest(p)for p in (input_path,tables_path)},
                 claim='Conservative enclosure of the six declared Cody J0 small-argument rational factor/combination variants with independent RN53/RN64 coefficient loads and arbitrary five-model arithmetic sites.',
                 limits='Absolute x in[1/4,4], order0, stated coefficient decimals and algebra only; dependency-free boxes may contain unrealizable combinations. No other Bessel program or Excel cutoff is excluded or fitted.')
    with(dest/'summary.json').open('x')as out:json.dump(summary,out,indent=2);out.write('\n')
    print(json.dumps(summary))


def as_bound(raw,tables,zform,form):
    """Fixed A&S polynomial, with every operation rounded independently."""
    x=Q(number(raw['args'][0]));order=int(number(raw['args'][1]));b=Boxes(tables,x);a=b.input
    if zform==0:
        scaled=b.op('div',a,b.literal(3));z=b.op('mul',scaled,scaled)
    elif zform==1:z=b.op('div',b.op('mul',a,a),b.literal(9))
    else:z=b.op('mul',b.op('mul',a,a),b.load('ONE_NINTH'))
    cs=[b.load('AS_J'+str(order),i)for i in range(7)]
    if form==0:result=b.horner(z,cs[::-1])
    elif form==1:
        z2=b.op('mul',z,z);z4=b.op('mul',z2,z2)
        low=b.op('add',b.op('add',cs[0],b.op('mul',cs[1],z)),b.op('mul',b.op('add',cs[2],b.op('mul',cs[3],z)),z2))
        high=b.op('add',b.op('add',cs[4],b.op('mul',cs[5],z)),b.op('mul',cs[6],z2))
        result=b.op('add',low,b.op('mul',high,z4))
    else:
        scheme,addition=divmod(form-2,3);powers={0:b.literal(1),1:z}
        def power(k):
            if k not in powers:
                if scheme==0 or k%2:powers[k]=b.op('mul',power(k-1),z)
                else:
                    half=power(k//2);powers[k]=b.op('mul',half,half)
            return powers[k]
        terms=[b.op('mul',c,power(k))for k,c in enumerate(cs)]
        if addition==1:terms.reverse()
        if addition==2:
            while len(terms)>1:terms=[b.op('add',terms[i],terms[i+1])if i+1<len(terms)else terms[i]for i in range(0,len(terms),2)]
            result=terms[0]
        else:
            result=terms[0]
            for term in terms[1:]:result=b.op('add',result,term)
    if order:result=b.op('mul',a,result)
    lo,hi=b.bounds(result);published=(B64.rn(lo),B64.rn(hi));target=Q(number(raw['expected']))
    return dict(zform=zform,form=form,nodes=b.nodes,output=result,published_bounds=list(map(pair,published)),target=pair(target),excluded=not(published[0]<=target<=published[1]))


def as_run(input_path,tables_path,destination):
    start=time.monotonic();dest=Path(destination);dest.mkdir();source=digest(Path(__file__))
    original=[json.loads(line)for line in Path(input_path).read_text().splitlines()];tables=json.loads(Path(tables_path).read_text())
    scoped=[];outside=[]
    for index,row in enumerate(original):
        a=row['args'];ok=row['function']=='BESSELJ'and row['result_index']is None and len(a)==2 and all(isinstance(v,str)and v.startswith('0x')for v in a)and row['expected'].startswith('0x')
        if ok:ok=number(a[1])in(0,1)and Q(1,4)<=Q(number(a[0]))<=3 and abs(Q(number(row['expected'])))>=Q(1,1<<1022)
        if ok:scoped.append(dict(trace_index=index,raw=row))
        else:outside.append(index)
    counts=Counter();observations=Counter(int(number(x['raw']['args'][1]))for x in scoped)
    require(set(observations)=={0,1},'empty A&S order scope');alive={order:set(range(24))for order in (0,1)}
    with(dest/'scope.json').open('x')as out:json.dump(dict(scoped=[x['trace_index']for x in scoped],outside=outside),out)
    with(dest/'bounds.jsonl').open('x')as out:
        for row in scoped:
            proofs=[];order=int(number(row['raw']['args'][1]))
            for zform,form in itertools.product(range(3),range(8)):
                proof=as_bound(row['raw'],tables,zform,form);proofs.append(proof)
                counts['coefficient_load_nodes']+=sum(n['op']=='coefficient'for n in proof['nodes'])
                counts['arithmetic_box_nodes']+=sum(n['op']in('add','sub','mul','div')for n in proof['nodes'])
                counts['boxes']+=1;counts['excluded_boxes']+=proof['excluded']
                if proof['excluded']:alive[order].discard(8*zform+form)
            out.write(json.dumps(dict(**row,proofs=proofs))+'\n')
    require(digest(Path(__file__))==source,'source changed during A&S study')  # ubs:ignore — public integrity hashes, not secrets
    summary=dict(complete=True,recorded=len(original),scoped=len(scoped),outside_scope=len(outside),observations={str(k):v for k,v in observations.items()},counts=counts,
                 necessary_survivors={str(k):sorted(v)for k,v in alive.items()},seconds=time.monotonic()-start,source_sha256=source,input_sha256={str(p):digest(p)for p in(input_path,tables_path)},
                 limits='Twenty-four labelled variants per order for A&S9.4.1/9.4.4 on positive x in[1/4,3]. Independent RN53/RN64 loads and arbitrary five-model arithmetic sites. Dependency-free hulls are conservative supersets; no remainder, fit, seed producer, cutoff or complete Bessel program is inferred.')
    with(dest/'summary.json').open('x')as out:json.dump(summary,out,indent=2);out.write('\n')
    print(json.dumps(summary))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='as-j01':as_run(*sys.argv[2:])
    else:run(*sys.argv[1:])
