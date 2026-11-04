"""Independent rational verification of finite LINEST coefficient programs."""
from collections import Counter
from fractions import Fraction as R
from functools import lru_cache
import itertools
import json
from pathlib import Path
import sys
import time

from exact import BINARY64 as B,binary64,operation,two
from check import records,require,digest

MODELS=('strict','stored_pc64','resident_pc64','stored_pc53','resident_pc53')


@lru_cache(maxsize=131072)
def arithmetic(kind,model,a,b):return operation(a,b,kind,model)


def extract(row):
    if row['function']!='LINEST'or row['result_index']not in ([1,1],[1,2]):return None,'surface_or_selected_statistic'
    a=row['args']
    if len(a)<2:return None,'implicit_predictor'
    def vector(arg):
        if not isinstance(arg,list):raise ValueError('not vector')
        if all(isinstance(r,list)for r in arg):
            if len(arg)>1 and any(len(r)!=1 for r in arg):raise ValueError('not vector')
            arg=list(itertools.chain.from_iterable(arg))
        if not all(isinstance(x,str)and x.startswith('0x')for x in arg):raise ValueError('not numeric')
        return list(map(binary64,arg))
    try:y,x=map(vector,a[:2])
    except ValueError:return None,'outside_numeric_single_predictor'
    if len(y)!=len(x)or not(2<=len(x)<=32):return None,'shape_or_length'
    if any(v!=0 and not two(-100)<=abs(v)<=two(100)for v in x+y):return None,'magnitude_bound'
    if len(a)<3:constant=True
    elif isinstance(a[2],bool):constant=a[2]
    elif isinstance(a[2],str)and a[2].startswith('0x'):constant=binary64(a[2])!=0
    else:return None,'const_coercion_outside_scope'
    if (constant and len(set(x))==1)or(not constant and all(v==0 for v in x)):return None,'not_exact_full_rank'
    if not row['expected'].startswith('0x'):return None,'non_numeric_expected'
    binary64(row['expected'])
    return (x,y,constant),None


def evaluate(p,data,column):
    xs,ys,constant=data;models=p['models']
    def op(site,kind,a,b):return arithmetic(kind,models[site],a,b)
    def minus(site,a,b):return op(site,'add',a,-b)
    def total(items,order,site):
        if order==2:
            if len(items)==0:return R(0)
            if len(items)==1:return items[0]
            middle=len(items)//2
            return op(site,'add',total(items[:middle],2,site),total(items[middle:],2,site))
        result=R(0)
        for item in items if order==0 else items[::-1]:result=op(site,'add',result,item)
        return result
    size=R(len(xs));linear=p['linear_order'];squared=p['square_order']
    sx,sy=total(xs,linear,0),total(ys,linear,0)
    mx,my=op(1,'div',sx,size),op(1,'div',sy,size)
    if p['algorithm']==0 and constant:
        a=[minus(1,x,mx)for x in xs];b=[minus(1,y,my)for y in ys]
        xx=[op(2,'mul',x,x)for x in a];xy=[op(2,'mul',x,y)for x,y in zip(a,b)]
        denominator,numerator=total(xx,squared,3),total(xy,squared,3)
    else:
        xx=[op(2,'mul',x,x)for x in xs];xy=[op(2,'mul',x,y)for x,y in zip(xs,ys)]
        denominator,numerator=total(xx,squared,3),total(xy,squared,3)
        if constant:
            if p['algorithm']==1:
                denominator=minus(4,denominator,op(2,'mul',sx,mx))
                numerator=minus(4,numerator,op(2,'mul',sx,my))
            else:
                denominator=minus(4,op(4,'mul',size,denominator),op(2,'mul',sx,sx))
                numerator=minus(4,op(4,'mul',size,numerator),op(2,'mul',sx,sy))
    if denominator==0:return None,'error:Div0'
    slope=op(4,'div',numerator,denominator)
    value=slope if column==1 else minus(4,my,op(2,'mul',slope,mx))if constant else R(0)
    published=B.round(value)
    if published is None:return value,'error:Num'
    if published==0:return value,'0x0000000000000000'
    if abs(published)<two(-1022):return value,'fault:NotImplemented'
    return value,f'0x{B.floor_coordinate(abs(published))|((1<<63)if published<0 else 0):016x}'


def check(input_path,directory,output):
    start=time.monotonic();directory=Path(directory);source=Path(__file__)
    hashes={str(p):digest(p)for p in (source,source.with_name('exact.py'),source.with_name('check.py'))}
    rows=list(records(input_path));scoped=[];excluded=[];coordinates=[];scope=Counter()
    for index,raw in enumerate(rows):
        data,reason=extract(raw)
        if reason:excluded.append(dict(trace_index=index,raw=raw,reason=reason));scope[reason]+=1
        else:scoped.append(dict(trace_index=index,raw=raw));coordinates.append(data);scope['scoped']+=1
    require(scoped==list(records(directory/'scoped.jsonl'))and excluded==list(records(directory/'excluded.jsonl')),'changed inputs or scope')
    ps=[dict(id=i,algorithm=a,linear_order=l,square_order=s,models=list(m))for i,(a,l,s,m)in enumerate(itertools.product(range(3),range(3),range(3),itertools.product(MODELS,repeat=5)))]
    require(ps==list(records(directory/'profiles.jsonl'))and len(ps)==84375,'changed program grammar')
    summary=json.loads((directory/'summary.json').read_text());require(summary['complete']is True,'partial run')
    killed=set();executions=0;counts=Counter();per_witness=Counter()
    for record in records(directory/'first-kills.jsonl'):
        identity,order=record['candidate_id'],record['order']
        require(identity not in killed and 0<=identity<len(ps)and 0<=order<len(scoped),'duplicate or unknown program/witness')
        require(record['trace_index']==scoped[order]['trace_index']and record['expected']==scoped[order]['raw']['expected'],'changed identity or answer')
        for earlier in range(order+1):
            row=scoped[earlier]['raw'];root,got=evaluate(ps[identity],coordinates[earlier],row['result_index'][1]);executions+=1
            if earlier<order:require(got==row['expected'],'not the first contradiction')
            else:
                expected_root=[str(root.numerator),str(root.denominator)]if root is not None else None
                require(record['root']==expected_root and record['got']==got and got!=row['expected'],'invalid full-width first kill')
        killed.add(identity);per_witness[order]+=1;counts['certified_first_kills']+=1
    alive=json.loads((directory/'survivors.json').read_text())
    require(len(set(alive))==len(alive)and killed.isdisjoint(alive)and killed|set(alive)==set(range(len(ps))),'unsearched program')
    for identity in alive:
        for row,data in zip(scoped,coordinates):
            _,got=evaluate(ps[identity],data,row['raw']['result_index'][1]);executions+=1
            require(got==row['raw']['expected'],'invalid survivor')
    remaining=len(ps);witnesses=list(records(directory/'witnesses.jsonl'))
    require(len(witnesses)==len(scoped),'missing witness accounting')
    for order,(witness,row)in enumerate(zip(witnesses,scoped)):
        remaining-=per_witness[order]
        require(witness==dict(order=order,trace_index=row['trace_index'],killed=per_witness[order],remaining=remaining),'wrong elimination counts')
    require(summary['programs']==len(ps)and summary['recorded']==len(rows)and summary['scoped']==len(scoped)and
            summary['scope']==scope and summary['executions']==executions and summary['survivors']==len(alive),'incorrect summary')
    require(all(digest(p)==h for p,h in hashes.items()),'checker changed during run')  # ubs:ignore — public integrity hashes, not secrets
    receipt=dict(complete=True,programs=len(ps),scoped=len(scoped),survivors=len(alive),executions=executions,
                 counts=counts,arithmetic_cache=arithmetic.cache_info()._asdict(),seconds=time.monotonic()-start,
                 source_sha256=hashes,input_sha256={str(p):digest(p)for p in (Path(input_path),*sorted(directory.iterdir()))if p.is_file()},
                 limits='Specified single-predictor coefficient algebra only. No QR, multiple predictors, LOGEST, numerical rank threshold, fitted coefficient or complete Excel program is inferred.')
    with open(output,'x')as out:json.dump(receipt,out,indent=2);out.write('\n')
    print(json.dumps(receipt))


def rank_data(raw):
    if raw['function'] not in ('LINEST','LOGEST') or raw['result_index'] != [4,2]:
        return None,'not_degrees_of_freedom'
    args=raw['args']
    if len(args)<2 or not isinstance(args[0],list) or not isinstance(args[1],list):return None,'not_explicit_numeric_arrays'
    try:
        y=args[0]
        if all(isinstance(r,list)for r in y):
            if not(len(y)==1 or all(len(r)==1 for r in y)):return None,'response_shape'
            y=[v for r in y for v in r]
        x=args[1]
        if not all(isinstance(r,list)for r in x):x=[[v]for v in x]
        if len(x)!=len(y)or not 2<=len(y)<=32 or not x or not 1<=len(x[0])<=8 or any(len(r)!=len(x[0])for r in x):return None,'matrix_shape'
        values=list(map(binary64,y+[v for r in x for v in r]))
        if any(v and not two(-100)<=abs(v)<=two(100)for v in values):return None,'magnitude_bound'
        columns=[[binary64(r[c])for r in x]for c in range(len(x[0]))]
        constant=True if len(args)<3 else args[2] if isinstance(args[2],bool)else binary64(args[2])!=0
        target=binary64(raw['expected'])
        if target.denominator!=1 or not 0<=target<=len(y):return None,'not_integer_degrees_of_freedom'
    except (ValueError,TypeError,OverflowError):return None,'non_numeric'
    return (columns,constant,int(target)),None


def rank_schur(profile,data):
    columns,constant,_=data;n=len(columns[0]);p=len(columns)
    # Gram matrix, with exact centering, followed by Schur complements.
    means=[sum(c)/n if constant else R(0)for c in columns]
    matrix=[[sum((a-means[i])*(b-means[j])for a,b in zip(columns[i],columns[j]))
             for j in range(p)]for i in range(p)]
    diagonal=[matrix[i][i]for i in range(p)];remaining=list(range(p));kept=[];dropped=[];first=None
    unit=R(1,2**52);gamma=n*unit/(1-n*unit)
    tolerance={'zero':R(0),'epsilon':unit,'epsilon_squared':unit*unit,
               'sqrt_epsilon':R(1,2**26),'gamma_n':gamma,'gamma_n_squared':gamma*gamma}[profile['threshold']]
    while remaining and len(kept)+int(constant)<n:
        index=max(remaining,key=lambda k:matrix[k][k]);remaining.remove(index)
        pivot=matrix[index][index]
        if first is None:first=pivot
        scales={'column':diagonal[index],'largest':max(diagonal),'first':first,'absolute':R(1)}
        if diagonal[index]==0:
            (kept if profile['zero']=='keep'else dropped).append(index)
        elif pivot<=tolerance*scales[profile['normalizer']]:dropped.append(index)
        else:
            kept.append(index)
            for i in remaining:
                for j in remaining:matrix[i][j]-=matrix[i][index]*matrix[index][j]/pivot
    dropped.extend(remaining)
    return n-len(kept)-int(constant),dict(kept=kept,dropped=dropped)


def check_rank(input_path,directory,output):
    directory=Path(directory);frozen={str(p):digest(p)for p in (Path(input_path),*directory.iterdir(),Path(__file__),Path(__file__).with_name('exact.py'))if p.is_file()}
    rows=[];scope=Counter()
    for index,raw in enumerate(records(input_path)):
        data,reason=rank_data(raw)
        if reason:scope[reason]+=1
        else:rows.append((index,raw,data));scope['scoped']+=1
    require(rows,'empty rank-policy scope')
    ps=[dict(id=i,threshold=t,normalizer=n,zero=z)for i,(t,n,z)in enumerate(itertools.product(
        ('zero','epsilon','epsilon_squared','sqrt_epsilon','gamma_n','gamma_n_squared'),
        ('column','largest','first','absolute'),('keep','drop')))]
    require(json.loads((directory/'profiles.json').read_text())==ps,'changed rank grammar')
    require(json.loads((directory/'scoped.json').read_text())==[dict(trace_index=i,raw=r)for i,r,_ in rows],'changed rank scope')
    seen=set();executions=0
    for kill in records(directory/'first-kills.jsonl'):
        identity,at=kill['candidate_id'],kill['order']
        require(identity not in seen and 0<=identity<len(ps) and 0<=at<len(rows),'bad rank kill identity')
        for prior in range(at+1):
            index,raw,data=rows[prior];got,trace=rank_schur(ps[identity],data);executions+=1
            if prior<at:require(got==data[2],'not first rank contradiction')
            else:require(got==kill['got']!=data[2]==kill['expected']and kill['trace_index']==index and trace==kill['trace'],'invalid rank kill')
        seen.add(identity)
    alive=json.loads((directory/'survivors.json').read_text())
    require(len(set(alive))==len(alive) and seen.isdisjoint(alive) and seen|set(alive)==set(range(len(ps))),'unsearched rank policy')
    for identity in alive:
        for _,_,data in rows:
            require(rank_schur(ps[identity],data)[0]==data[2],'invalid rank survivor');executions+=1
    summary=json.loads((directory/'summary.json').read_text())
    require(summary['complete']and summary['programs']==len(ps)and summary['observations']==len(rows)and summary['scope']==scope and summary['executions']==executions and summary['survivors']==len(alive),'wrong rank summary')
    require(all(digest(p)==h for p,h in frozen.items()),'rank evidence/source changed')
    result=dict(complete=True,programs=len(ps),observations=len(rows),survivors=len(alive),executions=executions,input_sha256=frozen,
                limits='Independent rational Gram-matrix Schur-complement replay of fixed rank-policy hypotheses. Degrees of freedom only, not a floating QR or complete coefficient program; no fitted threshold or adoption follows.')
    with open(output,'x')as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(result))


def rounded_rank_profiles():
    """Physical precision/centering/normalization programs; no fitted tolerance."""
    axes=itertools.product(('binary32','binary64','extended64','stored_pc64','decimal15'),
        ('raw_moments','centered'),('none','power2','largest'),('largest','original'),
        ('keep','drop'),('schur','ldl'))
    return [dict(id=i,precision=p,centering=c,normalization=n,pivot=v,zero=z,factor=f)
            for i,(p,c,n,v,z,f) in enumerate(axes)]


def decimal15(value):
    """Exact nearest/even decimal quantization, with no floating logarithm."""
    if value==0:return R(0)
    sign=-1 if value<0 else 1
    value=abs(value)
    exponent=len(str(value.numerator))-len(str(value.denominator))
    power=lambda e:R(10**e) if e>=0 else R(1,10**-e)
    if value<power(exponent):exponent-=1
    quantum=power(exponent-14)
    scaled=value/quantum
    whole,remainder=divmod(scaled.numerator,scaled.denominator)
    if 2*remainder>scaled.denominator or 2*remainder==scaled.denominator and whole%2:whole+=1
    return sign*whole*quantum


def rounded_rank(profile,data):
    from exact import Lattice,EXTENDED64,exponent_of
    columns,constant,_=data;n=len(columns[0]);p=len(columns)
    lattice={'binary32':Lattice(24,-126,127),'binary64':B,'extended64':EXTENDED64}
    def store(value):
        mode=profile['precision']
        if mode=='decimal15':return decimal15(value)
        if mode=='stored_pc64':
            first=EXTENDED64.round(value)
            rounded=None if first is None else B.round(first)
        else:rounded=lattice[mode].round(value)
        if rounded is None:raise ArithmeticError('finite lattice overflow')
        return rounded
    add=lambda a,b:store(a+b)
    mul=lambda a,b:store(a*b)
    div=lambda a,b:store(a/b)
    def total(values):
        result=R(0)
        for value in values:result=add(result,value)
        return result
    values=[]
    for column in columns:
        largest=max(map(abs,column))
        scale=(two(exponent_of(largest)) if profile['normalization']=='power2' and largest else
               largest if profile['normalization']=='largest' and largest else R(1))
        values.append([store(x/scale)for x in column])
    sums=[total(c)for c in values]
    if constant and profile['centering']=='centered':
        means=[div(s,R(n))for s in sums]
        values=[[add(x,-mean)for x in c]for c,mean in zip(values,means)]
    matrix=[[total(mul(a,b)for a,b in zip(values[i],values[j]))for j in range(p)]for i in range(p)]
    if constant and profile['centering']=='raw_moments':
        matrix=[[add(matrix[i][j],-div(mul(sums[i],sums[j]),R(n)))for j in range(p)]for i in range(p)]
    diagonal=[matrix[i][i]for i in range(p)];remaining=list(range(p));kept=[];dropped=[]
    lower=[];pivots=[];pivot_values=[]
    while remaining and len(kept)+int(constant)<n:
        index=max(remaining,key=lambda i:matrix[i][i]) if profile['pivot']=='largest' else remaining[0]
        remaining.remove(index)
        row=[]
        if profile['factor']=='ldl':
            for k in range(len(kept)):
                value=matrix[index][kept[k]]
                for j in range(k):value=add(value,-mul(mul(row[j],lower[k][j]),pivots[j]))
                row.append(div(value,pivots[k]))
            pivot=matrix[index][index]
            for j in range(len(kept)):pivot=add(pivot,-mul(mul(row[j],row[j]),pivots[j]))
        else:pivot=matrix[index][index]
        pivot_values.append([index,str(pivot)])
        if diagonal[index]==0:
            # Preserve the separately observed zero-column count convention.
            if profile['zero']=='keep':kept.append(index);lower.append([R(0)]*len(lower));pivots.append(R(1))
            else:dropped.append(index)
        elif pivot<=0:dropped.append(index)
        else:
            kept.append(index);lower.append(row);pivots.append(pivot)
            if profile['factor']=='schur':
                for i in remaining:
                    for j in remaining:matrix[i][j]=add(matrix[i][j],-div(mul(matrix[i][index],matrix[index][j]),pivot))
    return n-len(kept)-int(constant),dict(kept=kept,dropped=dropped+remaining,pivots=pivot_values)


def run_rounded_rank(input_path,output):
    """Necessary rank constraints only; every label is killed or fully replayed."""
    import gzip
    output=Path(output);output.mkdir()
    source=Path(input_path);raw=source.read_bytes()
    if source.suffix=='.gz':raw=gzip.decompress(raw)
    loaded=json.loads(raw)
    rows=[];scope=Counter()
    for entry in loaded:
        original=entry['raw'];data,reason=rank_data(original)
        if reason:scope[reason]+=1
        else:rows.append((entry['trace_index'],original,data));scope['scoped']+=1
    require(rows,'empty rounded-rank scope')
    ps=rounded_rank_profiles();alive=list(range(len(ps)));executions=0;kills=[]
    (output/'profiles.json').write_text(json.dumps(ps,indent=2)+'\n')
    for order,(index,original,data)in enumerate(rows):
        following=[]
        for identity in alive:
            try:got,trace=rounded_rank(ps[identity],data)
            except (ArithmeticError,ValueError)as error:got=None;trace=dict(refusal=str(error))
            executions+=1
            if got==data[2]:following.append(identity)
            else:kills.append(dict(candidate_id=identity,order=order,trace_index=index,
                got=got,expected=data[2],trace=trace,original=original))
        alive=following
    (output/'first-kills.jsonl').write_text(''.join(json.dumps(r)+'\n'for r in kills))
    result=dict(complete=True,program_labels=len(ps),observations=len(rows),scope=dict(scope),
        survivors=alive,executions=executions,input_sha256={str(source):digest(source)},
        source_sha256={str(p):digest(p)for p in (Path(__file__),Path(__file__).with_name('exact.py'))},
        sources='IEEE754 binary32/binary64/x87 precision; published decimal15 precision as an explicitly hypothetical internal arithmetic; centered/raw Gram identities and Schur/LDL factorizations. Original MIT rational implementation, no external code copied.',
        limits='Degrees-of-freedom constraints only, not a coefficient/QR program. No fitted tolerance: pivots are compared only with zero. '
            'Physical precision, centering, scaling, pivot and zero-column labels are fixed before evaluation. '
            'Decimal arithmetic is a hypothesis, not an observed Excel intermediate. Scope is the unchanged original rank corpus; all existing answers are discovery. '
            'First contradictions exclude only these finite programs; no complete rank algorithm family or numerical leaf is ruled out or adopted.')
    (output/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='rank-rounded':run_rounded_rank(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='rank':check_rank(*sys.argv[2:])
    else:check(*sys.argv[1:])
