"""Independent rational check of R11's entire cumulative-principal grammar.

Uses the existing nearest/even reference, not the racer or Fpu implementation.
Reconstructs each declared program's expression and checks original input bits,
every recorded node, every first contradiction, and coverage of all 90,000 IDs.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path
import sys
import time

from check import records, require, digest, tagged_value, write_result
from exact import BINARY64 as B, binary64, operation

MODELS = ('strict','stored_pc64','resident_pc64','stored_pc53','resident_pc53')
SOURCE = digest(Path(__file__))


def node(kind, model, *operands):
    return hashlib.sha256(json.dumps([kind,model,*operands],separators=(',',':')).encode()).hexdigest()


def program_catalogue():
    for i, values in enumerate(itertools.product(range(2), MODELS, range(3), MODELS,
            (False,True), range(4), MODELS, range(3), MODELS)):
        yield dict(id=i,**dict(zip(('power','growth_model','payment_form','payment_model',
            'payment_store','principal_form','principal_model','sum_order','sum_model'),values)))


def coordinates(row):
    require(row['function']=='CUMPRINC'and row['result_index'] is None,'unsupported surface')
    args=row['args'];require(len(args)==6,'wrong argument count')
    require(all(isinstance(a,str)and a.startswith('0x')and len(a)==18 for a in args),'argument type lost')
    a=list(map(binary64,args));require(a[0]>0 and a[2]>0 and a[5]in (0,1),'domain violation')
    require(all(v.denominator==1 for v in (a[1],a[3],a[4])),'noninteger loop coordinate')
    n,s,e=map(int,(a[1],a[3],a[4]));require(1<=s<=e<=n<=512,'unroll bound violation')
    return [a[0],a[2],-a[2],Fraction(0),Fraction(0)],n,s,e,bool(a[5])


@lru_cache(maxsize=8192)
def expression(serialized, n, first, last, beginning):
    p=json.loads(serialized)
    r,pv,negative_pv=(node('arg','strict',i)for i in range(3))
    one,negative_one,zero=(node('const','strict',bits)for bits in (0x3ff0000000000000,0xbff0000000000000,0))
    gm,pm,im,sm=(p[k]for k in ('growth_model','payment_model','principal_model','sum_model'))
    add=lambda a,b,m:node('add',m,a,b)
    mul=lambda a,b,m:node('mul',m,a,b)
    div=lambda a,b,m:node('div',m,a,b)
    q=add(one,r,gm)
    @lru_cache(maxsize=None)
    def power(k):
        if p['power']==0:
            value=one
            for _ in range(k):value=mul(value,q,gm)
            return value
        value=one;base=q
        while k:
            if k%2:value=mul(value,base,gm)
            k//=2
            if k:base=mul(base,base,gm)
        return value
    g=power(n)
    if p['payment_form']==0:
        payment=div(mul(mul(negative_pv,r,pm),g,pm),add(g,negative_one,pm),pm)
    else:
        denominator=add(one,div(negative_one,g,pm),pm)
        payment=(div(mul(negative_pv,r,pm),denominator,pm)if p['payment_form']==1
                 else mul(negative_pv,div(r,denominator,pm),pm))
    if beginning:payment=div(payment,q,pm)
    if p['payment_store']:payment=mul(payment,one,'strict')
    terms=[];previous=zero;balance=pv
    for i in range(1,last+1):
        if p['principal_form']==0:
            principal=payment if beginning and i==1 else div(payment,power(n-i+1),im)
        elif p['principal_form']==1:
            if i==1:
                principal=payment if beginning else add(payment,mul(pv,r,im),im)
            elif beginning and i==2:
                principal=add(mul(payment,q,im),mul(pv,r,im),im)
            else:principal=mul(previous,q,im)
        else:
            principal=payment if beginning and i==1 else add(payment,mul(balance,r,im),im)
            if i<last:
                balance=(add(balance,principal,im)if (beginning and i==1)or p['principal_form']==2
                         else add(mul(balance,q,im),payment,im))
        if i>=first:terms.append(principal)
        previous=principal
    def pairwise(items):
        if not items:return zero
        if len(items)==1:return items[0]
        middle=len(items)//2
        return add(pairwise(items[:middle]),pairwise(items[middle:]),sm)
    if p['sum_order']==2:return pairwise(terms)
    if p['sum_order']==1:terms.reverse()
    total=zero
    for term in terms:total=add(total,term,sm)
    return total


def check(trace, directory, output, expected_rows=540):
    started=time.monotonic();directory=Path(directory)
    expected_rows=int(expected_rows)
    original=list(records(trace));require(expected_rows>0 and len(original)==expected_rows,'incomplete declared discovery observations')
    cases=[coordinates(row)for row in original]
    programs=list(records(directory/'programs.jsonl'));declared=list(program_catalogue())
    require(programs==declared and len(programs)==90000,'grammar bounds/order differ')
    summaries={};counts=Counter()
    for record in records(directory/'traces.jsonl'):
        i=record['row_index'];args,*_=cases[i]
        require(record['trace_id']not in summaries,'duplicate trace ID')
        require(list(map(binary64,record['input_bits']))==args,'instantiated input bits changed')
        graph=record['graph'];values=[];shapes=[]
        require(len(graph['nodes'])==len(record['values']),'missing full-width values')
        for j,(item,text)in enumerate(zip(graph['nodes'],record['values'])):
            require(item['model']in MODELS,'undeclared arithmetic model')
            require(len(item['op'])==1,'invalid operation')
            kind,operands=next(iter(item['op'].items()));model=item['model']
            if kind=='arg':
                require(model=='strict'and operands in range(3),'undeclared argument node')
                value=args[operands];shape=node(kind,model,operands)
            elif kind=='const':
                bits=operands['bits'];require(model=='strict'and bits in (0,0x3ff0000000000000,0xbff0000000000000),'undeclared constant')
                value=binary64(bits);shape=node(kind,model,bits)
            else:
                require(kind in ('add','mul','div')and len(operands)==2 and all(0<=a<j for a in operands),'invalid DAG edge')
                a,b=operands;value=operation(values[a],values[b],kind,model)
                shape=node(kind,model,shapes[a],shapes[b])
            expected_tag='ext80-le:'if kind not in ('arg','const')and model.startswith('resident_')else 'binary64:'
            require(text.startswith(expected_tag),'wrong retained value representation')
            require(tagged_value(text)==value,f'arithmetic differs at trace{record["trace_id"]} node{j}')
            values.append(value);shapes.append(shape);counts['nodes_checked']+=1
        root=graph['output'];require(0<=root<len(values),'invalid root')
        actual=binary64(record['actual']);require(actual<0,'nonnegative/nonfinite publication outside this unit')
        require(B.round(values[root])==actual,'final binary64 publication differs')
        summaries[record['trace_id']]=(i,shapes[root],record['actual'])
        counts['trace_groups_checked']+=1
    # Several syntactically different graphs may have the same exact expression.
    # Preserve their trace IDs while sharing the verified publication lookup.
    outcomes={}
    for i,shape,actual in summaries.values():
        key=(i,shape)
        if key in outcomes:require(outcomes[key]==actual,'equal expression has unequal result')
        outcomes[key]=actual
    killed={};histogram=Counter()
    for kill in records(directory/'first-kills.jsonl'):
        pid=kill['program_id'];require(pid not in killed and pid in range(90000),'duplicate/unknown program kill')
        i=kill['row_index'];require(kill['raw']==original[i]and kill['expected']==original[i]['expected'],'observation binding differs')
        p=programs[pid];serialized=json.dumps(p,sort_keys=True,separators=(',',':'))
        for row_index in range(i+1):
            _,n,s,e,t=cases[row_index];shape=expression(serialized,n,s,e,t)
            key=(row_index,shape);require(key in outcomes,'missing earlier candidate replay')
            actual=outcomes[key]
            require((actual==original[row_index]['expected'])==(row_index<i),'not the first exact contradiction')
            if row_index==i:
                require(summaries[kill['trace_id']]==(i,shape,actual),'wrong killing trace')
                require(kill['actual']==actual,'changed killing publication')
        killed[pid]=i;histogram[i]+=1
    survivors=json.loads((directory/'survivors.json').read_text())
    require(not survivors and set(killed)==set(range(90000)),'nonzero/unverified survivors')
    declared_summary=json.loads((directory/'summary.json').read_text())
    require(declared_summary['original_rows']==expected_rows,'declared observation count differs')
    require(declared_summary['first_kill_histogram']=={str(k):v for k,v in sorted(histogram.items())},'kill counts differ')
    require(digest(Path(__file__))==SOURCE,'checker changed during execution')  # ubs:ignore — public source hashes, no secret comparison
    write_result(output,dict(complete=True,programs=90000,survivors=0,original_available_rows=expected_rows,
        **counts,first_kill_histogram=histogram,seconds=time.monotonic()-started,
        cumulative_checker_sha256=SOURCE,inputs={str(p):digest(p)for p in (Path(trace),directory/'programs.jsonl',directory/'first-kills.jsonl',directory/'traces.jsonl')},
        claim='Every declared program is independently bound to its exact original arguments, full arithmetic trace and first exact contradiction. No fitted coefficient, answer component, approximate equality or unsearched beam.',
        limitations='Only these integer-period geometric-growth/payment/principal/accumulation sketches and five arithmetic models. Other helpers, associations, stores, loops, stopping rules, profiles and input transports remain unsearched. Available observations are discovery, not fresh validation. Signed-zero/nonfinite general semantics are outside this ordinary nonzero scope.'))


if __name__=='__main__':check(*sys.argv[1:])
