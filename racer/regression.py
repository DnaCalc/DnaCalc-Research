"""Finite regression arithmetic programs, with exact rounding.

Original algebra, MIT. This extends the existing arithmetic research machinery;
it is not an alternate application evaluator. No coefficient or rank tolerance
is selected or fitted. Statistical programs explicitly condition on original
QR, residual and mean operands collected through the registered-route control.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
import itertools
import json
from pathlib import Path
import struct
import time

try:
    from .bessel import MODELS,operation,pair,digest,require
    from .intervals import B64,number,pow2
except ImportError:
    from bessel import MODELS,operation,pair,digest,require
    from intervals import B64,number,pow2


def extract(row):
    if row['function']!='LINEST' or row['result_index']not in ([1,1],[1,2]):return None,'surface_or_selected_statistic'
    a=row['args']
    if len(a)<2:return None,'implicit_predictor'
    def vector(arg):
        if not isinstance(arg,list):raise ValueError('not an array')
        if all(isinstance(r,list)for r in arg):
            if len(arg)>1 and any(len(r)!=1 for r in arg):raise ValueError('multiple predictors')
            arg=list(itertools.chain.from_iterable(arg))
        if not all(isinstance(x,str)and x.startswith('0x')for x in arg):raise ValueError('non-numeric vector')
        return [Q(number(x))for x in arg]
    try:y,x=map(vector,a[:2])
    except ValueError:return None,'outside_numeric_single_predictor'
    if len(y)!=len(x)or not(2<=len(x)<=32):return None,'shape_or_length'
    if any(v and not pow2(-100)<=abs(v)<=pow2(100)for v in x+y):return None,'magnitude_bound'
    if len(a)<3:constant=True
    elif isinstance(a[2],bool):constant=a[2]
    elif isinstance(a[2],str)and a[2].startswith('0x'):constant=number(a[2])!=0
    else:return None,'const_coercion_outside_scope'
    if (constant and len(set(x))==1)or(not constant and all(v==0 for v in x)):return None,'not_exact_full_rank'
    if not row['expected'].startswith('0x'):return None,'non_numeric_expected'
    return (x,y,constant),None


def profiles():
    for identity,(algorithm,linear,squared,models)in enumerate(itertools.product(range(3),range(3),range(3),itertools.product(MODELS,repeat=5))):
        yield dict(id=identity,algorithm=algorithm,linear_order=linear,square_order=squared,models=list(models))


def evaluate(profile,values,column):
    x,y,constant=values;m=profile['models'];algorithm=profile['algorithm']
    def op(kind,site,a,b):return operation(kind,m[site],a,b)
    def sub(site,a,b):return op('add',site,a,-b)
    def reduce(items,order,site):
        if order==2:
            if not items:return Q(0)
            if len(items)==1:return items[0]
            middle=len(items)//2
            return op('add',site,reduce(items[:middle],2,site),reduce(items[middle:],2,site))
        total=Q(0)
        for value in items if order==0 else items[::-1]:total=op('add',site,total,value)
        return total
    n=Q(len(x));linear=profile['linear_order'];squared=profile['square_order']
    sx,sy=reduce(x,linear,0),reduce(y,linear,0)
    mx,my=op('div',1,sx,n),op('div',1,sy,n)
    if algorithm==0 and constant:
        dx=[sub(1,v,mx)for v in x];dy=[sub(1,v,my)for v in y]
        xx=[op('mul',2,v,v)for v in dx];xy=[op('mul',2,v,w)for v,w in zip(dx,dy)]
        denominator=reduce(xx,squared,3);numerator=reduce(xy,squared,3)
    else:
        xx=[op('mul',2,v,v)for v in x];xy=[op('mul',2,v,w)for v,w in zip(x,y)]
        denominator=reduce(xx,squared,3);numerator=reduce(xy,squared,3)
        if constant:
            if algorithm==1:
                denominator=sub(4,denominator,op('mul',2,sx,mx))
                numerator=sub(4,numerator,op('mul',2,sx,my))
            else:
                denominator=sub(4,op('mul',4,n,denominator),op('mul',2,sx,sx))
                numerator=sub(4,op('mul',4,n,numerator),op('mul',2,sx,sy))
    if denominator==0:return None,'error:Div0'
    slope=op('div',4,numerator,denominator)
    value=slope if column==1 else sub(4,my,op('mul',2,slope,mx))if constant else Q(0)
    published=B64.rn(value)
    if published is None:return value,'error:Num'
    if published==0:return value,'0x0000000000000000'
    if abs(published)<pow2(-1022):return value,'fault:NotImplemented'
    return value,'0x'+struct.pack('>d',float(published)).hex()


def run(input_path,destination):
    start=time.monotonic();dest=Path(destination);dest.mkdir()
    original=[json.loads(line)for line in Path(input_path).read_text().splitlines()]
    scope=Counter();rows=[];values=[]
    with(dest/'scoped.jsonl').open('x')as good,(dest/'excluded.jsonl').open('x')as bad:
        for index,row in enumerate(original):
            coordinates,reason=extract(row)
            if reason:bad.write(json.dumps(dict(trace_index=index,raw=row,reason=reason))+'\n');scope[reason]+=1
            else:
                rows.append((index,row));values.append(coordinates);scope['scoped']+=1
                good.write(json.dumps(dict(trace_index=index,raw=row))+'\n')
    ps=list(profiles());require(len(ps)==84375 and rows, 'incomplete grammar or empty regression scope')
    with(dest/'profiles.jsonl').open('x')as out:
        for p in ps:out.write(json.dumps(p)+'\n')
    alive=list(range(len(ps)));executions=0;source=digest(Path(__file__))
    with(dest/'first-kills.jsonl').open('x')as kills,(dest/'witnesses.jsonl').open('x')as witnesses:
        for order,((index,row),data)in enumerate(zip(rows,values)):
            before=len(alive);next_alive=[]
            for identity in alive:
                value,got=evaluate(ps[identity],data,row['result_index'][1]);executions+=1
                if got==row['expected']:next_alive.append(identity)
                else:kills.write(json.dumps(dict(candidate_id=identity,order=order,trace_index=index,
                     root=pair(value)if value is not None else None,got=got,expected=row['expected']))+'\n')
            alive=next_alive
            witnesses.write(json.dumps(dict(order=order,trace_index=index,killed=before-len(alive),remaining=len(alive)))+'\n')
    with(dest/'survivors.json').open('x')as out:json.dump(alive,out);out.write('\n')
    require(digest(Path(__file__))==source, 'source changed during run')  # ubs:ignore — public integrity hashes, not secrets
    result=dict(complete=True,programs=len(ps),scoped=len(rows),recorded=len(original),scope=scope,
                executions=executions,survivors=len(alive),seconds=time.monotonic()-start,
                input_sha256={str(input_path):digest(input_path)},source_sha256=source,
                claim='Finite centered/raw-moment single-predictor least-squares coefficient programs with coherent stores and precisions. No coefficient is fitted.',
                limits='LINEST coefficient/intercept publications, exact full-rank numeric vectors and fixed algebra only. '
                       'No QR, multiple predictor, standard-error/statistic, LOGEST, rank-tolerance or complete kernel identification.')
    with(dest/'summary.json').open('x')as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(result))


def rank_coordinates(raw):
    """Read only explicit degrees-of-freedom observations, with exact numeric inputs."""
    if raw['function'] not in ('LINEST','LOGEST') or raw['result_index'] != [4,2]:
        return None, 'not_degrees_of_freedom'
    args = raw['args']
    if len(args)<2 or not isinstance(args[0],list) or not isinstance(args[1],list):
        return None, 'not_explicit_numeric_arrays'
    try:
        y = args[0]
        if all(isinstance(r,list) for r in y):
            if not (len(y)==1 or all(len(r)==1 for r in y)):return None,'response_shape'
            y = list(itertools.chain.from_iterable(y))
        x = args[1]
        if not all(isinstance(r,list) for r in x): x = [[v] for v in x]
        if len(x)!=len(y) or not(2<=len(y)<=32) or not x or not(1<=len(x[0])<=8):
            return None,'matrix_shape'
        if any(len(r)!=len(x[0]) for r in x):return None,'matrix_shape'
        values = [Q(number(v)) for v in y+list(itertools.chain.from_iterable(x))]
        if any(v and not pow2(-100)<=abs(v)<=pow2(100) for v in values):return None,'magnitude_bound'
        columns = [[Q(number(row[c])) for row in x] for c in range(len(x[0]))]
        constant = True if len(args)<3 else args[2] if isinstance(args[2],bool) else number(args[2])!=0
        expected = Q(number(raw['expected']))
        if expected.denominator!=1 or not 0<=expected<=len(y):return None,'not_integer_degrees_of_freedom'
    except (ValueError,TypeError,OverflowError,struct.error):return None,'non_numeric'
    return (columns,constant,int(expected)),None


def statistic_profiles():
    for identity,(model,root,form)in enumerate(itertools.product(MODELS,('binary64','coherent'),range(6))):
        yield dict(id=identity,model=model,root=root,form=form)


def statistic_value(profile,leaf,covariance_override=None):
    from intervals import E64,E53
    ss,covariance,df=(Q(number(leaf[k]))for k in ('ss','covariance','df'))
    if covariance_override is not None:
        covariance=covariance_override
    model=profile['model']
    def op(kind,a,b):return operation(kind,model,a,b)
    def root(x):
        if profile['root']=='binary64' or model=='strict':
            result=B64.rn_sqrt(B64.rn(x))
        else:
            fmt=E53 if model.endswith('pc53')else E64
            result=fmt.rn_sqrt(x)
            if model.startswith('stored'):result=B64.rn(result)
        require(result is not None,'nonfinite square root outside the finite study')
        return result
    form=profile['form']
    if form==0:value=root(op('div',op('mul',ss,covariance),df))
    elif form==1:value=root(op('mul',op('div',ss,df),covariance))
    elif form==2:value=root(op('mul',ss,op('div',covariance,df)))
    elif form==3:value=op('mul',root(op('div',ss,df)),root(covariance))
    elif form==4:value=op('mul',root(ss),root(op('div',covariance,df)))
    else:value=op('div',root(op('mul',ss,covariance)),root(df))
    value=B64.rn(value)
    if value is None:return 'error:Num'
    return '0x'+struct.pack('>d',float(value)).hex()


def statistic_run(input_path,destination):
    """Fixed exact standard-error consumers with explicitly conditional QR leaves."""
    from intervals import E64,E53
    import math
    import random
    start=time.monotonic();dest=Path(destination);dest.mkdir()
    source=Path(__file__);intervals=source.with_name('intervals.py')
    hashes={str(p):digest(p)for p in (source,intervals,Path(input_path))}
    # Binary64 controls independently use hardware sqrt; exact square and tie
    # controls exercise both extended precisions without a floating oracle.
    rng=random.Random(21003);controls=0
    for _ in range(1000):
        bits=rng.randrange(1,0x7ff0000000000000)
        x=struct.unpack('>d',bits.to_bytes(8,'big'))[0]
        require(float(B64.rn_sqrt(Q(x)))==math.sqrt(x),'binary64 root control')
        controls+=1
    for fmt in (B64,E64,E53):
        for exponent in (-100,0,100):
            for offset in range(12):
                lo=pow2(exponent)+offset*pow2(exponent-fmt.p+1)
                hi=lo+pow2(exponent-fmt.p+1);mid=(lo+hi)/2
                for squared,wanted in ((lo*lo,lo),(hi*hi,hi),(mid*mid,lo if offset%2==0 else hi)):
                    require(fmt.rn_sqrt(squared)==wanted,'exact square/root midpoint control');controls+=1
    rows=[json.loads(line)for line in Path(input_path).read_text().splitlines()]
    ps=list(statistic_profiles());counts={};first={};scope=Counter();baseline_controls=0
    with(dest/'predictions.jsonl').open('x')as out:
        for index,record in enumerate(rows):
            raw,leaf=record['raw'],record['leaf'];role=str(leaf['role'])
            values=[Q(number(leaf[k]))for k in ('ss','covariance','df')]
            if any(v and not pow2(-500)<=abs(v)<=pow2(500)for v in values):
                scope['outside_research_magnitude_bound']+=1;continue
            require(values[0]>=0 and values[1]>=0 and values[2]>0,'invalid finite leaf domain')
            baseline=dict(model='strict',root='binary64',form=0 if role=='0'else 1)
            require(statistic_value(baseline,leaf)==raw['actual']==leaf['value'],'original standard-error arithmetic changed')
            baseline_controls+=1;scope['scoped']+=1
            stats=counts.setdefault(role,[[0]*5 for _ in ps]);kills=first.setdefault(role,[None]*len(ps))
            was=raw['actual']==raw['expected'];predictions=[]
            for p in ps:
                prediction=statistic_value(p,leaf);predictions.append(prediction);now=prediction==raw['expected']
                for i,value in enumerate((1,now,was,now and not was,was and not now)):stats[p['id']][i]+=value
                if not now and kills[p['id']]is None:kills[p['id']]=dict(input_index=index,trace_index=record['trace_index'],id=raw['id'],prediction=prediction,expected=raw['expected'])
            out.write(json.dumps(dict(input_index=index,id=raw['id'],role=role,predictions=predictions,expected=raw['expected'],baseline=raw['actual']),separators=(',',':'))+'\n')
    require(scope['scoped']>0 and sum(scope.values())==len(rows),'incomplete statistic scope')
    require(all(digest(p)==h for p,h in hashes.items()),'statistic source/input changed')  # ubs:ignore — public evidence-integrity hashes, not authentication secrets
    summary=dict(complete=True,input_rows=len(rows),scope=scope,programs_per_role=len(ps),roles={'0':'predictor','1':'constant'},
        profiles=ps,counts_all_exact_baseline_gains_losses=counts,first_contradictions=first,
        exact_root_controls=controls,baseline_controls=baseline_controls,input_sha256=hashes,
        seconds=time.monotonic()-start,
        source='Defining standard-error identity sqrt(residual_sum_squares*covariance/degrees_of_freedom); IEEE754 and Intel x87 models; exact integer midpoint square roots.',
        limits='Six fixed algebraic expressions, five arithmetic modes and binary64/coherent square roots, separately counted by mathematical coefficient role. Current QR/covariance operands are conditional binary64 leaves, not identified Excel producers. No fitting, per-row selector, tolerance or complete regression-kernel claim.')
    (dest/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items()if k not in ('profiles','first_contradictions','input_sha256')}))


def factor_profiles():
    for identity,(shape,order,model,store)in enumerate(itertools.product(range(4),range(3),MODELS,(False,True))):
        yield dict(id=identity,shape=shape,order=order,model=model,store_factor=store)


def factor_value(profile,leaf):
    factors=leaf['factors'];w=list(map(lambda x:Q(number(x)),factors['weights']))
    inverse=[[Q(number(x))for x in row]for row in factors['inverse']]
    require(len(inverse)==len(w)and all(len(row)==len(w)for row in inverse),'invalid triangular-factor shape')
    model=profile['model'];order=profile['order'];shape=profile['shape'];size=len(w)
    def op(kind,a,b):return operation(kind,model,a,b)
    def total(values):
        values=list(values)
        if order==2 and len(values)>1:
            middle=len(values)//2;return op('add',total(values[:middle]),total(values[middle:]))
        value=Q(0)
        for x in values if order!=1 else values[::-1]:value=op('add',value,x)
        return value
    if shape==0:
        def covariance(a,b):return total(op('mul',inverse[a][k],inverse[b][k])for k in range(size))
        value=total(op('mul',op('mul',w[a],w[b]),covariance(a,b))for a in range(size)for b in range(size))
    else:
        start=0 if shape==1 else 1
        weighted=[total(op('mul',w[a],inverse[a][k])for a in range(start,size))for k in range(start,size)]
        value=total(op('mul',x,x)for x in weighted)
        if shape>=2:
            base=op('mul',inverse[0][0],inverse[0][0])if shape==2 else op('div',Q(1),Q(factors['observations']))
            value=op('add',base,value)
    return B64.rn(value)if profile['store_factor']else value


def factor_run(input_path,destination):
    """Quadratic form versus weighted inverse-R norms, with frozen QR factors."""
    start=time.monotonic();dest=Path(destination);dest.mkdir()
    rows=[json.loads(line)for line in Path(input_path).read_text().splitlines()]
    scoped=[r for r in rows if r['leaf']['role']==1]
    require(scoped,'empty intercept standard-error scope')
    producers=list(factor_profiles());consumers=list(statistic_profiles())
    profiles=[dict(id=i,producer=p,consumer=c)for i,(p,c)in enumerate(itertools.product(producers,consumers))]
    counts=[[0]*5 for _ in profiles];first=[None]*len(profiles);controls=0
    for at,record in enumerate(scoped):
        raw,leaf=record['raw'],record['leaf']
        baseline=dict(shape=0,order=0,model='strict',store_factor=False)
        require(factor_value(baseline,leaf)==Q(number(leaf['covariance'])),'original covariance-factor arithmetic changed')
        controls+=1;values=[factor_value(p,leaf)for p in producers];was=raw['actual']==raw['expected']
        for p in profiles:
            covariance=values[p['producer']['id']]
            require(covariance is not None and covariance>=0,'outside finite positive covariance-factor study')
            prediction=statistic_value(p['consumer'],leaf,covariance);now=prediction==raw['expected'];tally=counts[p['id']]
            for i,v in enumerate((1,now,was,now and not was,was and not now)):tally[i]+=v
            if not now and first[p['id']]is None:first[p['id']]=dict(input_index=at,id=raw['id'],prediction=prediction,expected=raw['expected'])
    summary=dict(complete=True,input_rows=len(rows),scoped=len(scoped),excluded_other_roles=len(rows)-len(scoped),
        producers=len(producers),consumers=len(consumers),profiles=len(profiles),baseline_factor_controls=controls,
        counts_all_exact_baseline_gains_losses=counts,first_contradictions=first,
        input_sha256={str(p):digest(p)for p in (Path(input_path),Path(__file__),Path(__file__).with_name('intervals.py'))},seconds=time.monotonic()-start,
        source='Linear least-squares covariance identity w*R^-1*(R^-1)^T*w^T = sum((w*R^-1)_k^2); centered intercept identity 1/n plus centered predictor term; IEEE754 exact algebra and square-root rounding.',
        limits='Conditional original binary64 inverse-R and weight operands, with120 fixed factor producers and60 standard-error consumers. Four mathematical forms, three summation orders, five arithmetic modes and a final store choice. No fitted coefficient, cutoff, per-row selector or claim of identifying a full QR/regression algorithm.')
    (dest/'profiles.json').write_text(json.dumps(profiles)+'\n')
    (dest/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items()if k not in ('counts_all_exact_baseline_gains_losses','first_contradictions','input_sha256')}))


def intercept_solve_factor(profile, leaf):
    """Solve R' z = w directly, without first constructing (X'X)^-1."""
    matrix=[[Q(number(x))for x in row]for row in leaf['triangular']['r']]
    weights=[Q(number(x))for x in leaf['factors']['weights']]
    n=Q(leaf['factors']['observations']);model=profile['model'];order=profile['order']
    def op(kind,a,b):return operation(kind,model,a,b)
    def store(x):return B64.rn(x)if profile['local_store']else x
    def total(terms):
        terms=list(terms)
        if order=='balanced'and len(terms)>1:
            k=len(terms)//2;return op('add',total(terms[:k]),total(terms[k:]))
        result=Q(0)
        for t in terms[::-1]if order=='reverse'else terms:result=op('add',result,t)
        return result
    start=int(profile['shape']!='full');z=[Q(0)]*len(weights)
    for j in range(start,len(weights)):
        indices=list(range(start,j))
        if order=='reverse':indices.reverse()
        if profile['solve']=='sum_then_subtract':
            value=op('add',weights[j],-total(op('mul',matrix[i][j],z[i])for i in range(start,j)))
            z[j]=store(op('div',value,matrix[j][j]))
        elif profile['solve']=='running_difference':
            value=weights[j]
            for i in indices:value=op('add',value,-op('mul',matrix[i][j],z[i]))
            z[j]=store(op('div',value,matrix[j][j]))
        else:
            value=op('div',weights[j],matrix[j][j])
            for i in indices:value=op('add',value,-op('mul',op('div',matrix[i][j],matrix[j][j]),z[i]))
            z[j]=store(value)
    factor=total(op('mul',v,v)for v in z[start:])
    if start:
        base=op('div',1,n)if profile['shape']=='centered_1/n'else op('mul',Q(number(leaf['factors']['inverse'][0][0])),Q(number(leaf['factors']['inverse'][0][0])))
        factor=op('add',base,factor)
    return B64.rn(factor)if profile['store_factor']else factor


def intercept_solve_run(input_path,destination):
    """First contradictions for direct and centered inverse-R mean solves."""
    rows=[r for r in map(json.loads,Path(input_path).read_text().splitlines())if r['leaf']['role']==1]
    require(rows and all(r['leaf']['triangular']for r in rows),'no triangular intercept operands')
    dest=Path(destination);dest.mkdir(exist_ok=False);profiles=[];survivors=[]
    descriptions=itertools.product(('full','centered_1/n','centered_inverse00'),
        ('sum_then_subtract','running_difference','divide_first'),('forward','reverse','balanced'),
        MODELS,(True,False),(True,False))
    consumers=list(statistic_profiles())
    for shape,solve,order,model,local_store,store_factor in descriptions:
        producer=dict(shape=shape,solve=solve,order=order,model=model,local_store=local_store,store_factor=store_factor)
        active=set(range(len(consumers)));first=[None]*len(consumers);checked=[0]*len(consumers)
        for index,row in enumerate(rows):
            factor=intercept_solve_factor(producer,row['leaf'])
            for i in sorted(active):
                prediction=statistic_value(consumers[i],row['leaf'],factor);checked[i]+=1
                if prediction!=row['raw']['expected']:
                    first[i]=dict(input_index=index,id=row['raw']['id'],expected=row['raw']['expected'],prediction=prediction)
                    active.remove(i)
            if not active:break
        for i,consumer in enumerate(consumers):
            identity=len(profiles)
            if first[i]is None:survivors.append(identity)
            profiles.append(dict(id=identity,producer=producer,consumer=consumer,checked=checked[i],first_contradiction=first[i]))
    result=dict(complete=True,rows=len(rows),programs=profiles,survivors=survivors,
        source_sha256={str(p):digest(p)for p in(Path(__file__),Path(input_path),Path(__file__).with_name('intervals.py'))},
        limits='Direct triangular solve of R transpose against the intercept weights, or the centered predictor submatrix with1/n or inverseR00 squared. Original QR and residual operands remain conditional. No new Excel answers, constant fit or complete regression claim. A first exact contradiction refutes the whole program.')
    (dest/'summary.json').write_text(json.dumps(result,separators=(',',':'))+'\n')
    print(json.dumps(dict(programs=len(profiles),survivors=survivors,rows=len(rows))))


def intercept_raw_moment_run(input_path,destination):
    """Single-retained-predictor intercept errors from raw second moments."""
    all_rows=[r for r in map(json.loads,Path(input_path).read_text().splitlines())if r['leaf']['role']==1]
    rows=[]
    for row in all_rows:
        leaf=row['leaf'];tri=leaf['triangular'];n=leaf['factors']['observations']
        if len(tri['order'])!=2 or tri['order'][0]is not None:continue
        index=tri['order'][1];raw=row['raw']['args'][1]
        if not isinstance(raw,list):continue
        if raw and not isinstance(raw[0],list):raw=[raw]
        if len(raw)==n:xs=[Q(number(r[index]))for r in raw]
        elif len(raw[index])==n:xs=list(map(lambda v:Q(number(v)),raw[index]))
        else:continue
        rows.append((row,xs))
    require(rows,'no raw one-predictor statistics');dest=Path(destination);dest.mkdir(exist_ok=False)
    programs=[];survivors=[]
    for norm,average,outer,order,form in itertools.product(MODELS,MODELS,MODELS,('forward','reverse','balanced'),range(7)):
        first=None;checked=0
        def op(k,m,a,b):return operation(k,m,Q(a),Q(b))
        def total(xs):
            if order=='balanced'and len(xs)>1:
                k=len(xs)//2;return op('add',norm,total(xs[:k]),total(xs[k:]))
            value=Q(0)
            for x in xs[::-1]if order=='reverse'else xs:value=op('add',norm,value,x)
            return value
        def root(x):return B64.rn_sqrt(B64.rn(x))
        for row,xs in rows:
            leaf=row['leaf'];n=Q(len(xs));ss,df=map(lambda k:Q(number(leaf[k])),('ss','df'))
            inv=Q(number(leaf['factors']['inverse'][1][1]));r=Q(number(leaf['triangular']['r'][1][1]));mean=-Q(number(leaf['factors']['weights'][1]))
            rawsum=total([op('mul',norm,x,x)for x in xs]);rawmean=op('div',average,rawsum,n)
            covariance=B64.rn(inv*inv);sem=root(B64.rn(B64.rn(ss*covariance)/df))
            mse=op('div',outer,ss,df);sey=root(B64.rn(ss/df));sxx=op('mul',outer,r,r)
            centered=op('add',outer,op('div',average,1,n),op('div',outer,op('mul',outer,mean,mean),sxx))
            if form==0:result=op('mul',outer,sem,root(rawmean))
            elif form==1:result=root(op('mul',outer,op('mul',outer,sem,sem),rawmean))
            elif form==2:result=root(op('mul',outer,mse,op('div',outer,rawsum,op('mul',outer,n,sxx))))
            elif form==3:result=root(op('mul',outer,mse,centered))
            elif form==4:
                ratio=op('div',outer,op('mul',outer,n,op('mul',outer,mean,mean)),sxx)
                result=root(op('mul',outer,op('div',outer,mse,n),op('add',outer,1,ratio)))
            elif form==5:result=op('mul',outer,sey,root(centered))
            else:result=root(op('mul',outer,op('mul',outer,sey,sey),centered))
            prediction='0x'+struct.pack('>d',float(B64.rn(result))).hex();checked+=1
            if prediction!=row['raw']['expected']:
                first=dict(id=row['raw']['id'],expected=row['raw']['expected'],prediction=prediction);break
        identity=len(programs)
        if first is None:survivors.append(identity)
        programs.append(dict(id=identity,norm=norm,average=average,outer=outer,order=order,form=form,checked=checked,first_contradiction=first))
    result=dict(complete=True,original_intercept_rows=len(all_rows),single_retained_predictor_rows=len(rows),programs=programs,survivors=survivors,
        input_sha256={str(p):digest(p)for p in(Path(input_path),Path(__file__))},
        limits='Raw second-moment identity se(intercept)=se(slope)*sqrt(mean(x^2)), and centered one-predictor forms using the original R diagonal, residual and published standard errors. No multi-predictor extension, fit or fresh-validation claim.')
    (dest/'summary.json').write_text(json.dumps(result,separators=(',',':'))+'\n');print(json.dumps(dict(programs=len(programs),rows=len(rows),survivors=survivors)))


def rank_profiles():
    # All tolerances are mathematical IEEE roundoff scales, never inferred
    # from a residual or fitted to an Excel rank transition.
    for identity,(threshold,normalizer,zero) in enumerate(itertools.product(
            ('zero','epsilon','epsilon_squared','sqrt_epsilon','gamma_n','gamma_n_squared'),
            ('column','largest','first','absolute'),('keep','drop'))):
        yield dict(id=identity,threshold=threshold,normalizer=normalizer,zero=zero)


def rank_evaluate(profile, data):
    columns, constant, _ = data
    n = len(columns[0]); epsilon=Q(1,2**52);gamma=n*epsilon/(1-n*epsilon)
    tolerances = dict(zero=Q(0),epsilon=epsilon,epsilon_squared=epsilon**2,
                      sqrt_epsilon=Q(1,2**26),gamma_n=gamma,gamma_n_squared=gamma**2)
    values = [[x-sum(c)/n for x in c] if constant else list(c) for c in columns]
    initial = [sum(x*x for x in c) for c in values]
    live = list(range(len(columns)));basis=[];first=None;kept=[];dropped=[]
    while live and len(kept)+int(constant)<n:
        norms = {i:sum(x*x for x in values[i]) for i in live}
        index=max(live,key=lambda i:norms[i]); norm=norms[index];live.remove(index)
        if first is None:first=norm
        scale = {'column':initial[index], 'largest':max(initial), 'first':first, 'absolute':Q(1)}[profile['normalizer']]
        if initial[index]==0:
            if profile['zero']=='keep':kept.append(index)
            else:dropped.append(index)
            continue
        if norm<=tolerances[profile['threshold']]*scale:
            dropped.append(index);continue
        kept.append(index);basis.append(index)
        for other in live:
            coefficient=sum(a*b for a,b in zip(values[index],values[other]))/norm
            values[other]=[b-coefficient*a for a,b in zip(values[index],values[other])]
    dropped.extend(live)
    return n-len(kept)-int(constant), dict(kept=kept, dropped=dropped)


def rank_run(input_path, destination):
    dest=Path(destination);dest.mkdir();rows=[];scope=Counter()
    for index,raw in enumerate(map(json.loads,Path(input_path).read_text().splitlines())):
        data,reason=rank_coordinates(raw)
        if reason:scope[reason]+=1
        else:rows.append((index,raw,data));scope['scoped']+=1
    require(rows,'empty rank-policy scope')
    ps=list(rank_profiles());alive=list(range(len(ps)));executions=0
    (dest/'profiles.json').write_text(json.dumps(ps)+'\n')
    (dest/'scoped.json').write_text(json.dumps([dict(trace_index=i,raw=r) for i,r,_ in rows])+'\n')
    with (dest/'first-kills.jsonl').open('x')as kills:
        for order,(index,raw,data)in enumerate(rows):
            following=[]
            for identity in alive:
                got,trace=rank_evaluate(ps[identity],data);executions+=1
                if got==data[2]:following.append(identity)
                else:kills.write(json.dumps(dict(candidate_id=identity,order=order,trace_index=index,
                     got=got,expected=data[2],trace=trace))+'\n')
            alive=following
    (dest/'survivors.json').write_text(json.dumps(alive)+'\n')
    result=dict(complete=True,programs=len(ps),observations=len(rows),survivors=len(alive),executions=executions,
        scope=scope,source_sha256=digest(Path(__file__)),input_sha256={str(input_path):digest(input_path)},
        limits='Exact rational pivoted Gram-Schmidt geometry with fixed IEEE roundoff-scale comparisons. Degrees of freedom only; no floating QR implementation, coefficient program, fitted tolerance or reusable Excel rank rule is identified. Surviving necessary rank policies cannot be adopted as numeric programs.')
    (dest/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


def triangular_inverse(matrix, algorithm, model, local_store):
    """Upper/non-unit inverse: original back substitution or LAPACK DTRTI2.

    LAPACK3.12.1, 6ec7f2bc4ecf4c4a93496aa2fa519575bc0e39ca,
    SRC/dtrti2.f and BLAS/SRC/{dtrmv,dscal}.f, modified BSD.
    Copyright Universities of Tennessee, California Berkeley, Colorado Denver.
    See notices/LAPACK.txt. Array writes always store binary64.
    """
    size = len(matrix)
    require(size > 0 and all(len(r) == size for r in matrix), 'invalid triangular shape')
    require(all(matrix[i][i] != 0 for i in range(size)), 'zero triangular diagonal')
    def op(kind, a, b): return operation(kind, model, a, b)
    def store(x):
        value = B64.rn(x)
        require(value is not None, 'nonfinite inverse array element')
        return value
    if algorithm == 'back_substitution':
        inverse = [[Q(0) for _ in range(size)] for _ in range(size)]
        for j in range(size):
            inverse[j][j] = store(op('div', Q(1), matrix[j][j]))
            for i in range(j-1, -1, -1):
                total = Q(0)
                for k in range(i+1, j+1):
                    total = op('add', total, op('mul', matrix[i][k], inverse[k][j]))
                    if local_store: total = store(total)
                inverse[i][j] = store(op('div', -total, matrix[i][i]))
        return inverse
    require(algorithm == 'lapack_dtrti2', 'unknown triangular program')
    inverse = [list(row) for row in matrix]
    for j in range(size):
        inverse[j][j] = store(op('div', Q(1), inverse[j][j]))
        scale = -inverse[j][j]
        # DTRMV Upper/NoTranspose/NonUnit, unit stride, in source order.
        for k in range(j):
            if inverse[k][j] != 0:
                value = inverse[k][j]
                for i in range(k):
                    inverse[i][j] = store(op('add', inverse[i][j], op('mul', value, inverse[i][k])))
                inverse[k][j] = store(op('mul', inverse[k][j], inverse[k][k]))
        for i in range(j):
            inverse[i][j] = store(op('mul', scale, inverse[i][j]))
    return inverse


def inverse_profiles():
    for identity, (algorithm, model, local_store) in enumerate(itertools.product(
            ('back_substitution', 'lapack_dtrti2'), MODELS, (False, True))):
        yield dict(id=identity, algorithm=algorithm, model=model, local_store=local_store)


def inverse_covariance(profile, leaf, inverse):
    if leaf['role'] == 1:
        copy = dict(leaf, factors=dict(leaf['factors'], inverse=[
            ['0x'+struct.pack('>d', float(x)).hex() for x in row] for row in inverse]))
        return factor_value(profile, copy)
    original = leaf['triangular']
    predictor = original['predictors']-leaf['column']
    position = original['order'].index(predictor)
    values = inverse[position]
    model = profile['model']
    def total(xs):
        xs = list(xs)
        if profile['order'] == 2 and len(xs) > 1:
            middle = len(xs)//2
            return operation('add', model, total(xs[:middle]), total(xs[middle:]))
        result = Q(0)
        for x in xs if profile['order'] != 1 else xs[::-1]:
            result = operation('add', model, result, x)
        return result
    result = total(operation('mul', model, x, x) for x in values)
    return B64.rn(result) if profile['store_factor'] else result


def inverse_run(input_path, destination, mode='original_consumers'):
    """Whole inverse/covariance producers with explicitly conditional QR leaves."""
    require(mode in ('original_consumers', 'all_consumers'), 'unknown inverse consumer scope')
    start = time.monotonic(); dest = Path(destination); dest.mkdir()
    rows = [json.loads(line) for line in Path(input_path).read_text().splitlines()]
    require(rows and all(r['leaf'].get('triangular') for r in rows), 'missing original R observations')
    inverses = list(inverse_profiles()); factors = list(factor_profiles())
    consumers = list(statistic_profiles())
    grammar = dict(inverse_programs=inverses, factor_programs=factors,
                   consumer_programs=consumers if mode == 'all_consumers' else 'original role-specific strict consumer')
    (dest/'grammar.json').write_text(json.dumps(grammar, indent=2)+'\n')
    controls = 0; scope = Counter(r['leaf']['role'] for r in rows)
    counts = {role: [[0]*5 for _ in range(len(inverses)*len(factors)*(len(consumers) if mode == 'all_consumers' else 1))]
              for role in scope}
    first = {role: [None for _ in values] for role, values in counts.items()}
    for at, record in enumerate(rows):
        raw, leaf = record['raw'], record['leaf']; role = leaf['role']
        original = leaf['triangular']; matrix = [[Q(number(x)) for x in row] for row in original['r']]
        observed = [[Q(number(x)) for x in row] for row in original['inverse']]
        baseline = triangular_inverse(matrix, 'back_substitution', 'strict', False)
        require(baseline == observed, 'original inverse arithmetic changed')
        baseline_factor = dict(shape=0, order=0, model='strict', store_factor=False)
        require(inverse_covariance(baseline_factor, leaf, baseline) == Q(number(leaf['covariance'])),
                'original inverse/covariance route changed')
        baseline_consumer = dict(model='strict', root='binary64', form=role)
        require(statistic_value(baseline_consumer, leaf) == raw['actual'], 'original statistic route changed')
        controls += 1; was = raw['actual'] == raw['expected']
        predictions = {}
        for ip in inverses:
            try: inverse = triangular_inverse(matrix, ip['algorithm'], ip['model'], ip['local_store'])
            except (ValueError, ZeroDivisionError, OverflowError): inverse = None
            for fp in factors:
                try: covariance = inverse_covariance(fp, leaf, inverse) if inverse is not None else None
                except (ValueError, ZeroDivisionError, OverflowError): covariance = None
                used = consumers if mode == 'all_consumers' else [baseline_consumer]
                for ci, cp in enumerate(used):
                    identity = (ip['id']*len(factors)+fp['id'])*len(used)+ci
                    key = (covariance, cp['model'], cp['root'], cp['form'])
                    if key not in predictions:
                        try:
                            require(covariance is not None and covariance >= 0, 'nonfinite covariance')
                            prediction = statistic_value(cp, leaf, covariance)
                        except (ValueError, ZeroDivisionError, OverflowError): prediction = 'error:Num'
                        predictions[key] = prediction
                    prediction = predictions[key]; now = prediction == raw['expected']
                    for index, value in enumerate((1, now, was, now and not was, was and not now)):
                        counts[role][identity][index] += value
                    if not now and first[role][identity] is None:
                        first[role][identity] = [at, prediction]
    rankings = {role: sorted(range(len(values)), key=lambda i: (values[i][4], -values[i][1], i))
                for role, values in counts.items()}
    result = dict(complete=True, input_rows=len(rows), role_counts=dict(scope), original_route_controls=controls,
        counts_all_exact_baseline_gains_losses=counts, first_contradictions=first,
        zero_loss_improvements={role: [i for i,c in enumerate(values) if c[3] and not c[4]] for role,values in counts.items()},
        full_survivors={role: [i for i,c in enumerate(values) if c[0] == c[1]] for role,values in counts.items()},
        best={role: [dict(id=i, counts=counts[role][i]) for i in ids[:10]] for role,ids in rankings.items()},
        inputs_sha256={str(p): digest(p) for p in (Path(input_path), Path(__file__), Path(__file__).with_name('intervals.py'))},
        seconds=time.monotonic()-start, mode=mode,
        source='LAPACK3.12.1 6ec7f2bc4ecf4c4a93496aa2fa519575bc0e39ca DTRTI2/DTRMV/DSCAL, modified BSD; least-squares covariance identities and IEEE arithmetic.',
        limits='Original observed binary64 QR R, means, residual sums of squares and degrees of freedom remain conditional. Two inverse algorithms, five precisions/local stores,120 existing covariance labels and declared consumers. Array writes store binary64. No coefficient, threshold, rank tolerance or row-specific program fitted. A survivor still requires full application/bank checks.')
    (dest/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('counts_all_exact_baseline_gains_losses','first_contradictions','inputs_sha256')}))


def delivered_rank_constraints(delivery, destination):
    """Bind delivered rank readings to their exact reference-cell operands.

    Report necessary intervals for simple relative-residual comparisons.
    An interval is an observed constraint, never a fitted implementation.
    """
    import re
    import sys
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checker'))
    from questions_r13 import delivery_files,check_workbook_readings
    from corpus_cases import outcome
    delivery=Path(delivery);dest=Path(destination);dest.mkdir(exist_ok=False)
    files=delivery_files(delivery)
    question_paths=[p for p in(delivery/'packet').rglob('workbook-questions.json')
                    if json.loads(p.read_text())['sheet'].endswith('_regression')]
    require(len(question_paths)==1,'expected one delivered regression question book')
    question_path=question_paths[0]
    questions=json.loads(question_path.read_text())
    answers=[json.loads(line)for line in(delivery/'answers.jsonl').read_text().splitlines()
             if json.loads(line).get('observation',{}).get('sheet')==questions['sheet']]
    checked=check_workbook_readings(delivery,question_path,answers)
    values={a['question']['row']:a['observation']['excel_raw']for a in answers}
    predictions={p['id']:p for p in map(json.loads,(question_path.parent/'application-predictions.jsonl').read_text().splitlines())}
    rows=[];observations=[];groups={};unselected=Counter()
    expression=r'INDEX\((LINEST|LOGEST)\((D\d+:D\d+),CHOOSE\(\{1,2\},(D\d+:D\d+),(D\d+:D\d+)\),(TRUE|FALSE),TRUE\),([1-5]),([1-3])\)'
    def vector(reference):
        first,last=map(lambda x:int(x[1:]),reference.split(':'))
        return [['0x'+values[i]['number']]for i in range(first,last+1)]
    for answer in answers:
        q=answer['question']
        if q['role']!='question':continue
        match=re.fullmatch(expression,q['formula'])
        if match is None:
            unselected[q.get('function','unidentified')]+=1
            continue
        name,y,x1,x2,constant,r,c=match.groups();constant=constant=='TRUE';r,c=int(r),int(c)
        a,b=vector(x1),vector(x2);matrix=[a[i]+b[i]for i in range(len(a))];response=vector(y)
        require(len(a)==len(b)==len(response)==4,'rank reference shape changed')
        expected=answer['observation']['excel_raw']
        raw=dict(function=name,id=q['id'],args=[response,matrix,constant,True],result_index=[r,c],
                 expected=outcome(expected),actual=predictions[q['id']]['calc_func_prediction'],
                 original_question=q,source=str(delivery/'answers.jsonl'),original_observation=answer['observation'])
        rows.append(raw)
        if (r,c)!=(4,2):continue
        df=Q(number(raw['expected']));retained=4-int(df)-int(constant)
        require(df.denominator==1 and retained in (1,2),'unexpected observed rank')
        columns=[[Q(number(v[0]))for v in a],[Q(number(v[0]))for v in b]]
        if constant:columns=[[v-sum(col)/4 for v in col]for col in columns]
        xx,yy=[sum(v*v for v in col)for col in columns];xy=sum(x*y for x,y in zip(*columns))
        require(0<xx<yy and xy*xy<xx*yy,'not the expected exact full-rank pivot geometry')
        residual=xx-xy*xy/yy
        metrics=dict(column=residual/xx,largest=residual/yy,absolute=residual)
        reverse=q.get('reverse',False)
        observation=dict(id=q['id'],function=name,constant=constant,reverse=reverse,
                         changed_operand=q.get('changed_operand',q.get('input_bits')),df=int(df),retained_predictors=retained,
                         necessary_threshold_coordinates={k:[str(v.numerator),str(v.denominator)]for k,v in metrics.items()})
        observations.append(observation)
        group=(name,constant,reverse);bounds=groups.setdefault(group,{k:dict(drop=[],keep=[])for k in metrics})
        for k,v in metrics.items():bounds[k]['drop'if retained==1 else'keep'].append((v,q['id']))
    intervals=[]
    for group,bounds in groups.items():
        for normalizer,values_ in bounds.items():
            low=max(values_['drop'],default=(Q(0),None));high=min(values_['keep'],default=(None,None))
            intervals.append(dict(function=group[0],constant=group[1],reverse=group[2],normalizer=normalizer,
                lower_closed=[str(low[0].numerator),str(low[0].denominator)],lower_witness=low[1],
                upper_open=None if high[0]is None else[str(high[0].numerator),str(high[0].denominator)],upper_witness=high[1],
                consistent=high[0]is None or low[0]<high[0],drop_observations=len(values_['drop']),keep_observations=len(values_['keep'])))
    (dest/'inputs.jsonl').write_text(''.join(json.dumps(r)+'\n'for r in rows))
    result=dict(complete=True,delivery=files,workbook=checked,rank_observations=len(observations),selected_calls=len(rows),
        observations=observations,necessary_intervals=intervals,unselected_questions=unselected,
        limits='Exact four-row, two-predictor residual geometry and observed degrees of freedom only. Each threshold interval is a necessary constraint on one predeclared comparison; it is not a fitted tolerance or a QR/coefficient implementation. Row reversal and const modes stay separate. All original numeric answers, predictions and echoes remain unchanged.')
    (dest/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(selected_calls=len(rows),rank_observations=len(observations),consistent_intervals=sum(i['consistent']for i in intervals),intervals=len(intervals))))


def service_design(mode='rank'):
    """One reusable design for numeric rank, reflector and residual questions."""
    import math
    cases = []
    if mode == 'rank':
        for constant, reverse, swap, scales, exponent in itertools.product(
                (False, True), (False, True), (False, True),
                ((0, 0), (-20, -20), (20, 20), (-10, 10), (10, -10)),
                (None, -30, -26, -24, -22, -20, -18, -16, -14, -12, -10, -8, -4, 0)):
            delta = 0. if exponent is None else 2.**exponent
            xs = [[math.ldexp(x, scales[0]), math.ldexp(2*x+(delta if i==3 else 0.), scales[1])]
                  for i, x in enumerate((1., 2., 5., 4.))]
            y = [[v] for v in (1., 2., 4., 7.)]
            if swap: xs = [r[::-1] for r in xs]
            if reverse: xs, y = xs[::-1], y[::-1]
            cases.append(dict(id=str(len(cases)), x=xs, y=y, constant=constant,
                              reverse=reverse, swap=swap, scales=scales, delta_exponent=exponent))
        selections = [(1, 1), (1, 2), (1, 3), (4, 2)]
    elif mode == 'decomposition':
        for n, scale, constant, reverse, response in itertools.product(
                (3, 4, 5, 8), (-20, 0, 20), (False, True), (False, True), range(4)):
            xs = [[math.ldexp(float(i+1), scale), math.ldexp(float((i*3+1)%5), scale)]
                  for i in range(n)]
            y = [[float(i == response)] for i in range(n)]
            if reverse: xs, y = xs[::-1], y[::-1]
            cases.append(dict(id=str(len(cases)), x=xs, y=y, constant=constant,
                              reverse=reverse, scale=scale, response=response))
        selections = [(r, c) for r in range(1, 6) for c in range(1, 4)]
    elif mode == 'residual':
        # Zero-mean orthogonal predictors separate residual arithmetic from the
        # intercept's translation. Basis responses read the fitted projection.
        for n, scale, reverse, response in itertools.product(
                (3, 4, 5, 8, 9, 16), (-20, 0, 20), (False, True), range(5)):
            x = [math.ldexp(float(2*i-n+1), scale) for i in range(n)]
            if response < 3:
                y = [float(i == response) for i in range(n)]
            elif response == 3:
                y = [float((i*3+1)%7-3) for i in range(n)]
            else:
                y = [1.0 if i in (0, n-1) else -2.0/(n-2) for i in range(n)]
            if reverse: x, y = x[::-1], y[::-1]
            cases.append(dict(id=str(len(cases)), x=[[v] for v in x], y=[[v] for v in y],
                              constant=True, reverse=reverse, scale=scale, response=response))
        selections = [(1, 1), (1, 2), (2, 1), (2, 2), (3, 2), (4, 2), (5, 2)]
    else:
        raise ValueError('unknown regression question mode')
    return cases, selections


def service_questions(checkout, output, cache, mode='rank', offline=False):
    """Ask the common controlled design; answers never replace an evaluator."""
    try:
        from .fingerprint import load_service, CachedService, write
    except ImportError:
        from fingerprint import load_service, CachedService, write
    client, _ = load_service(checkout)
    service = CachedService(client.Client(timeout=300), cache, offline)
    output = Path(output)
    output.mkdir(exist_ok=False)
    cases, selections = service_design(mode)
    write(output/'questions.json', dict(mode=mode, cases=cases, selections=selections,
          source='Original controlled numeric inputs, Excel question-service outputs; no implementation source copied.'))
    counts = Counter()
    with (output/'observations.jsonl').open('x') as stream:
        for function in ('LINEST', 'LOGEST'):
            arguments = [[client.number_array(
                [[2.0**v[0]] for v in c['y']] if function == 'LOGEST' and mode == 'residual' else c['y']),
                client.number_array(c['x']), client.boolean(c['constant']), client.boolean(True)] for c in cases]
            for selection in selections:
                answers = service.ask_rows(function, arguments, chunk_size=500, result_index=selection)
                for case, args, answer in zip(cases, arguments, answers):
                    row = dict(function=function, case=case['id'], result_index=selection,
                               args=args, answer=answer)
                    stream.write(json.dumps(row, separators=(',', ':'))+'\n')
                    counts[function+':'+answer['type']] += 1
                stream.flush()
                print(function, selection, len(answers), flush=True)
    result = dict(complete=True, cases=len(cases), calls=sum(counts.values()), counts=counts,
                  mode=mode, receipts=service.receipts, source_sha256=digest(Path(__file__)),
                  limits='Controlled numeric array transport through CHOOSE; not a statement about range coercion. '
                         'No rank rule, decomposition or numeric kernel is inferred from match counts.')
    write(output/'summary.json', result)
    print(json.dumps(result))


if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='intercept-raw-moment':intercept_raw_moment_run(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='intercept-solve':intercept_solve_run(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='service':service_questions(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='delivered-rank':delivered_rank_constraints(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='inverse':inverse_run(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='factors':factor_run(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='statistics':statistic_run(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='rank':rank_run(*sys.argv[2:]);raise SystemExit(0)
    parser=argparse.ArgumentParser();parser.add_argument('input');parser.add_argument('destination')
    args=parser.parse_args();run(args.input,args.destination)
