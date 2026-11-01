"""Independent constant-load, rational-box and symbolic-shape verification."""
from collections import Counter
from fractions import Fraction as R
import itertools
import json
from pathlib import Path
import re
import sys
import time

from exact import BINARY64 as B,EXTENDED64 as E,binary64,operation
from check import records,require,digest
from x87 import bits64,encode

MODELS=('strict','stored_pc64','resident_pc64','stored_pc53','resident_pc53')


def value(pair):return R(*map(int,pair))


def constants(path,source_path):
    data=json.loads(Path(path).read_text());source=Path(source_path)
    require(data['source_sha256']==digest(source),'changed published source')  # ubs:ignore — public integrity hashes, not secrets
    # Read only double-precision DATA facts from the primary source. No code
    # from that source is executed or imported by the verifier.
    code=''.join(line[6:]for line in source.read_text().splitlines()if line[:2]=='CD')
    facts={}
    for name,count in [('PJ0',7),('QJ0',5)]:
        match=re.search(r'\bDATA\s+'+name+r'/([^/]+)/',code)
        require(match is not None,'missing primary-source DATA fact')
        xs=[R(x.strip().replace('D','e'))for x in match[1].split(',')]
        require(len(xs)==count and len(data['tables'][name])==count,'wrong table length')
        for index,(x,entry)in enumerate(zip(xs,data['tables'][name])):facts[name,index]=(x,entry)
    for name in ('XJ0','XJ01','XJ02'):
        match=re.search(r'\b'+name+r'/([^/]+)/',code)
        require(match is not None,'missing primary-source split-root fact')
        facts[name,None]=(R(match[1].strip().replace('D','e')),data['constants'][name])
    result={}
    for key,(x,entry)in facts.items():
        require(R(entry['decimal'])==x,'decimal differs from published coefficient')
        a,b=B.round(x),E.round(x)
        require(entry['binary64']==f'{bits64(a):016x}'and entry['ext80_le']==encode(b).hex(),'incorrect literal rounding/load')
        result[key]=(min(a,b),max(a,b))
    return result


def declared_shape(factor,association):
    x=('input',True);z=('mul',x,x)
    def c(name,index=None):return('coefficient',name,index)
    def h(coeffs):
        result=coeffs[0]
        for coefficient in coeffs[1:]:result=('add',('mul',result,z),coefficient)
        return result
    numerator=h([c('PJ0',i)for i in (4,5,6,0,1,2,3)])
    denominator=h([('literal',R(1))]+[c('QJ0',i)for i in (4,0,1,2,3)])
    root=c('XJ0')
    if factor==0:
        split=('div',c('XJ01'),('literal',R(256)))
        left=('sub',('sub',x,split),c('XJ02'))
        product=('mul',left,('add',x,root))
    elif factor==1:product=('mul',('sub',x,root),('add',x,root))
    else:product=('sub',z,('mul',root,root))
    return('div',('mul',product,numerator),denominator)if association==0 else('mul',product,('div',numerator,denominator))


def verify(raw,proof,loads,expected_shape=None,literals=(1,256)):
    if expected_shape is None:
        factor,association=proof['factor'],proof['association']
        require(factor in range(3)and association in range(2),'undeclared algebra variant')
        expected_shape=declared_shape(factor,association)
    shapes=[];bounds=[];dependencies=[];counts=Counter()
    for index,node in enumerate(proof['nodes']):
        kind=node['op']
        if kind=='input':
            x=abs(binary64(raw['args'][0]));span=(x,x);shape=('input',True)
            require(node==dict(op='input',bounds=node['bounds'],absolute=True),'changed input operation')
            used=[]
        elif kind=='literal':
            x=R(node['decimal']);require(x in literals,'unlisted literal')
            span=(x,x);shape=('literal',x);used=[]
            require(set(node)=={'op','bounds','decimal'},'undeclared literal attribute')
        elif kind=='coefficient':
            key=(node['name'],node['index']);require(key in loads,'unpublished coefficient')
            span=loads[key];shape=('coefficient',*key);used=[];counts['coefficient_load_nodes']+=1
            require(set(node)=={'op','bounds','name','index'},'undeclared coefficient attribute')
        else:
            require(kind in ('add','sub','mul','div')and set(node)=={'op','bounds','a','b'},'unlisted arithmetic operation')
            a,b=node['a'],node['b'];require(0<=a<index and 0<=b<index,'forward graph edge')
            if kind=='div':require(bounds[b][0]*bounds[b][1]>0,'denominator crosses zero')
            endpoints=[operation(x,-y if kind=='sub'else y,'add'if kind=='sub'else kind,model)
                       for x,y,model in itertools.product(bounds[a],bounds[b],MODELS)]
            span=(min(endpoints),max(endpoints));shape=(kind,shapes[a],shapes[b]);used=[a,b]
            counts['arithmetic_box_nodes']+=1
        require(tuple(map(value,node['bounds']))==span,'incorrect full-width rounding box')
        bounds.append(span);shapes.append(shape);dependencies.append(used)
    output=proof['output'];require(0<=output<len(shapes),'missing output')
    require(shapes[output]==expected_shape,'not the independently declared expression')
    reachable=set();todo=[output]
    while todo:
        n=todo.pop()
        if n not in reachable:reachable.add(n);todo.extend(dependencies[n])
    require(len(reachable)==len(shapes),'unused evidence nodes')
    published=tuple(B.round(x)for x in bounds[output]);target=binary64(raw['expected'])
    require(tuple(map(value,proof['published_bounds']))==published and value(proof['target'])==target,'changed publication or answer')
    excluded=not(published[0]<=target<=published[1]);require(proof['excluded']==excluded,'invalid exclusion')
    counts['excluded_boxes']+=excluded;counts['boxes']+=1
    return counts


def check(input_path,tables_path,source_path,directory,output):
    start=time.monotonic();directory=Path(directory);source=Path(__file__)
    hashes={str(p):digest(p)for p in (source,source.with_name('exact.py'),source.with_name('check.py'),source.with_name('x87.py'))}
    rows=list(records(input_path));loads=constants(tables_path,source_path);scope=[]
    for index,row in enumerate(rows):
        a=row['args']
        if row['function']=='BESSELJ'and row['result_index']is None and len(a)==2 and all(isinstance(x,str)and x.startswith('0x')for x in a)and binary64(a[1])==0 and R(1,4)<=abs(binary64(a[0]))<=4 and row['expected'].startswith('0x'):
            scope.append(dict(trace_index=index,raw=row))
    proofs=list(records(directory/'bounds.jsonl'));require(len(proofs)==len(scope)and scope,'incomplete scoped observations')
    counts=Counter();alive=set(range(6))
    for expected,record in zip(scope,proofs):
        require(record['trace_index']==expected['trace_index']and record['raw']==expected['raw'],'changed source row')
        require([(p['factor'],p['association'])for p in record['proofs']]==list(itertools.product(range(3),range(2))),'incomplete algebra variants')
        for proof in record['proofs']:
            counts.update(verify(record['raw'],proof,loads))
            if proof['excluded']:alive.discard(2*proof['factor']+proof['association'])
    summary=json.loads((directory/'summary.json').read_text())
    require(summary['complete']is True and summary['recorded']==len(rows)and summary['scoped']==len(scope)and
            summary['outside_scope']==len(rows)-len(scope)and summary['counts']==counts and summary['necessary_survivors']==sorted(alive),'incorrect seed-family summary')
    require(all(digest(p)==h for p,h in hashes.items()),'checker changed during seed proof')  # ubs:ignore — public integrity hashes, not secrets
    result=dict(complete=True,scoped=len(scope),counts=counts,necessary_survivors=sorted(alive),seconds=time.monotonic()-start,
                source_sha256=hashes,input_sha256={str(p):digest(p)for p in (Path(input_path),Path(tables_path),Path(source_path),*sorted(directory.iterdir()))if p.is_file()},
                limits='Conservative all-site five-model boxes for the six declared Cody J0 rational variants and independently rounded published literals only. No fitted coefficient, other Bessel program or Excel cutoff is inferred.')
    with open(output,'x')as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(result))


def as_constants(path,page369,page370):
    data=json.loads(Path(path).read_text())
    require(data['source_sha256']=={Path(p).name:digest(p)for p in(page369,page370)},'changed A&S primary scans')  # ubs:ignore — public integrity hashes, not secrets
    # Separate transcription checked directly against the primary scans by the
    # read-only reviewer. Only coefficient facts enter the verifier.
    facts={('AS_J0',i):R(v)for i,v in enumerate(('1','-2.2499997','1.2656208','-.3163866','.0444479','-.0039444','.0002100'))}
    facts.update({('AS_J1',i):R(v)for i,v in enumerate(('.5','-.56249985','.21093573','-.03954289','.00443319','-.00031761','.00001109'))})
    facts['ONE_NINTH',None]=R(1,9)
    require(set(data['tables'])=={'AS_J0','AS_J1'}and all(len(v)==7 for v in data['tables'].values())and set(data['constants'])=={'ONE_NINTH'},'changed A&S table coverage')
    result={}
    for(name,index),x in facts.items():
        entry=data['constants'][name]if index is None else data['tables'][name][index]
        require(R(entry['decimal'])==x,'changed published A&S coefficient')
        a,b=B.round(x),E.round(x)
        require(entry['binary64']==f'{bits64(a):016x}'and entry['ext80_le']==encode(b).hex(),'incorrect A&S coefficient load')
        result[name,index]=(min(a,b),max(a,b))
    return result


def as_shape(order,zform,form):
    x=('input',True)
    def op(k,a,b):return(k,a,b)
    if zform==0:
        t=op('div',x,('literal',R(3)));z=op('mul',t,t)
    elif zform==1:z=op('div',op('mul',x,x),('literal',R(9)))
    else:z=op('mul',op('mul',x,x),('coefficient','ONE_NINTH',None))
    cs=[('coefficient','AS_J'+str(order),i)for i in range(7)]
    if form==0:
        y=cs[6]
        for i in range(5,-1,-1):y=op('add',op('mul',y,z),cs[i])
    elif form==1:
        square=op('mul',z,z)
        groups=[op('add',cs[i],op('mul',cs[i+1],z))for i in(0,2,4)]+[cs[6]]
        a=op('add',groups[0],op('mul',groups[1],square));b=op('add',groups[2],op('mul',groups[3],square))
        y=op('add',a,op('mul',b,op('mul',square,square)))
    else:
        power_kind,addition=divmod(form-2,3)
        def power(k):
            if k==0:return('literal',R(1))
            if k==1:return z
            if power_kind==1 and k%2==0:
                half=power(k//2);return op('mul',half,half)
            return op('mul',power(k-1),z)
        terms=[op('mul',c,power(i))for i,c in enumerate(cs)]
        if addition==1:terms=terms[::-1]
        def tree(xs):
            if len(xs)==1:return xs[0]
            half=(len(xs)+1)//2;return op('add',tree(xs[:half]),tree(xs[half:]))
        if addition==2:y=tree(terms)
        else:
            y=terms[0]
            for t in terms[1:]:y=op('add',y,t)
    return op('mul',x,y)if order else y


def as_check(input_path,tables_path,page369,page370,directory,output):
    start=time.monotonic();directory=Path(directory);source=Path(__file__)
    hashes={str(p):digest(p)for p in(source,source.with_name('exact.py'),source.with_name('check.py'),source.with_name('x87.py'))}
    rows=list(records(input_path));loads=as_constants(tables_path,page369,page370);scope=[]
    for index,row in enumerate(rows):
        a=row['args'];ok=row['function']=='BESSELJ'and row['result_index']is None and len(a)==2 and all(isinstance(v,str)and v.startswith('0x')for v in a)and row['expected'].startswith('0x')
        if ok:ok=binary64(a[1])in(0,1)and R(1,4)<=binary64(a[0])<=3 and abs(binary64(row['expected']))>=R(1,1<<1022)
        if ok:scope.append(dict(trace_index=index,raw=row))
    require({int(binary64(x['raw']['args'][1]))for x in scope}=={0,1},'empty A&S order scope')
    selected={x['trace_index']for x in scope}
    require(json.loads((directory/'scope.json').read_text())==dict(scoped=sorted(selected),outside=[i for i in range(len(rows))if i not in selected]),'changed A&S selection')
    proofs=list(records(directory/'bounds.jsonl'));require(len(proofs)==len(scope),'incomplete A&S rows')
    counts=Counter();observations=Counter();alive={order:set(range(24))for order in(0,1)}
    for expected,record in zip(scope,proofs):
        require(record['trace_index']==expected['trace_index']and record['raw']==expected['raw'],'changed A&S source observation')
        order=int(binary64(record['raw']['args'][1]));observations[str(order)]+=1
        require([(p['zform'],p['form'])for p in record['proofs']]==list(itertools.product(range(3),range(8))),'incomplete A&S grammar')
        for proof in record['proofs']:
            shape=as_shape(order,proof['zform'],proof['form']);counts.update(verify(record['raw'],proof,loads,shape,(1,3,9)))
            if proof['excluded']:alive[order].discard(8*proof['zform']+proof['form'])
    survivors={str(k):sorted(v)for k,v in alive.items()};summary=json.loads((directory/'summary.json').read_text())
    require(summary['complete']is True and summary['recorded']==len(rows)and summary['scoped']==len(scope)and summary['outside_scope']==len(rows)-len(scope)and summary['counts']==counts and summary['observations']==observations and summary['necessary_survivors']==survivors,'incorrect A&S summary')
    require(all(digest(p)==h for p,h in hashes.items()),'checker changed during A&S proof')  # ubs:ignore — public integrity hashes, not secrets
    result=dict(complete=True,scoped=len(scope),observations=observations,counts=counts,necessary_survivors=survivors,seconds=time.monotonic()-start,source_sha256=hashes,
                input_sha256={str(p):digest(p)for p in(Path(input_path),Path(tables_path),Path(page369),Path(page370),*sorted(directory.iterdir()))if p.is_file()},
                limits='A&S9.4.1/9.4.4 fixed polynomial graphs only, positive x in[1/4,3]. Conservative all-site/load enclosures; no analytic remainder, fitted constant, cutoff or complete Bessel kernel is inferred.')
    with open(output,'x')as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(result))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='as-j01':as_check(*sys.argv[2:])
    else:check(*sys.argv[1:])
