"""Check latent-payment balance DAGs and monotone lattice cuts independently.

The recorded payment defines an existential publication cell. Its values are
never installed as constants in a numerical program. No racer/Fpu imports.
"""
from collections import Counter
from fractions import Fraction as Q
import itertools
import json
from pathlib import Path
import sys
import time

from check import records, require, digest, graph_trace, tagged_value, write_result
from cumulative import coordinates, node, MODELS
from exact import BINARY64 as B, EXTENDED64 as E, Span, binary64, inverse, operation

SOURCE=digest(Path(__file__))
HELPER_SOURCE=digest(Path(__file__).parent/'cumulative.py')


def catalogue():
    fields=('update','payment_store','interest_model','principal_model','balance_model','sum_order','sum_model')
    return [dict(id=i,**dict(zip(fields,values)))for i,values in enumerate(itertools.product(
        range(2),(False,True),MODELS,MODELS,MODELS,range(3),MODELS))]


def expected_shape(p,start,end,beginning,interest_sum=False):
    rate,pv,latent=(node('arg','strict',i)for i in range(3))
    one=node('const','strict',0x3ff0000000000000);zero=node('const','strict',0)
    payment=node('mul','strict',latent,one)if p['payment_store']else latent
    q=node('add',p['balance_model'],one,rate)if p['update']==1 else one
    balance=pv;terms=[]
    for period in range(1,end+1):
        interest=zero if beginning and period==1 else node('mul',p['interest_model'],balance,rate)
        principal=(payment if beginning and period==1 else
            node('add',p['principal_model'],payment,interest))
        if period>=start:terms.append(interest if interest_sum else principal)
        if period<end:
            balance=(node('add',p['balance_model'],balance,principal)if p['update']==0 or beginning and period==1
                else node('add',p['balance_model'],node('mul',p['balance_model'],balance,q),payment))
    def paired(items):
        if not items:return zero
        if len(items)==1:return items[0]
        mid=len(items)//2
        return node('add',p['sum_model'],paired(items[:mid]),paired(items[mid:]))
    if p['sum_order']==2:result=paired(terms)
    else:
        if p['sum_order']==1:terms.reverse()
        result=zero
        for term in terms:result=node('add',p['sum_model'],result,term)
    if interest_sum:
        import struct
        count=int.from_bytes(struct.pack('>d',float(end-start+1)),'big')
        payments=node('mul',p['principal_model'],node('const','strict',count),payment)
        result=node('add',p['principal_model'],payments,result)
    return result


def check(trace,structure,directory,output):
    started=time.monotonic();directory=Path(directory)
    rows=list(records(trace));require(len(rows)==540,'missing original observations')
    for r in rows:coordinates(r)
    original=json.loads(Path(structure).read_text())['groups'];require(len(original)==60,'missing original loan groups')
    loans=list(records(directory/'loans.jsonl'));require(len(loans)==60,'missing publication domains')
    domains=[];all_rows=[]
    for i,(loan,old)in enumerate(zip(loans,original)):
        require(loan['loan_index']==i and loan['payment']==old['PMT']and loan['row_indices']==old['row_indices'],'loan/PMT binding differs')
        p=old['PMT'];a=p['args'];require(p['function']=='PMT'and len(a)==5 and binary64(a[3])==0,'unsupported PMT surface')
        y=binary64(p['expected']);require(y<0,'payment sign outside scope')
        bounds=E.coordinates(B.pullback(Span(y,y)));require(bounds is not None,'empty publication cell')
        low,high=bounds;values=list(map(tagged_value,loan['latent_payment_values']))
        require(len(values)==high-low+1 and values==[E.value(k)for k in range(low,high+1)],'payment lattice is incomplete, reordered or altered')
        require(len({bytes.fromhex(v.split(':')[1])[8:]for v in loan['latent_payment_values']})==1,'outside declared one-binade cell')
        for row_index in loan['row_indices']:
            r=rows[row_index];require(r['args'][:3]==a[:3]and r['args'][5]==a[4]and r['reading']==p['reading'],'loan input/environment differs')
        all_rows.extend(loan['row_indices']);domains.append(values)
    require(sorted(all_rows)==list(range(540)),'missing/duplicate original rows')
    free_file=directory/'initial-domains.json';free_domains={};interest_sum=False
    if free_file.exists():
        doc=json.loads(free_file.read_text())
        interest_sum=doc.get('interest_sum',False)
        require(doc['trace_sha256']==digest(Path(trace))and doc['structure_sha256']==digest(Path(structure)),'initial-domain source differs')  # ubs:ignore — public source identities, no secret comparison
        for d in doc['domains']:
            li=d['loan_index'];ri=loans[li]['row_indices'][0]
            require(d['row_index']==ri,'initial domain uses wrong first singleton')
            a,_,start,end,beginning=coordinates(rows[ri]);require(start==end==1,'first observation is not first singleton')
            y=binary64(rows[ri]['expected']);required=B.pullback(Span(y,y))
            interest=Q(0)if beginning else operation(a[0],a[1],'mul',d['interest_model'])
            if interest_sum:
                if d['sum_order']!=2:interest=operation(Q(0),interest,'add',d['sum_model'])
                required=inverse(required,interest,'add',True,d['principal_model'])
                required=inverse(required,Q(1),'mul',True,d['principal_model'])
            else:
                required=inverse(required,Q(0),'add',True,d['sum_model'])
                if not beginning:required=inverse(required,interest,'add',True,d['principal_model'])
            if d['payment_store']:required=B.pullback(required)
            if not interest_sum:
                direct=B.pullback(Span(y,y))
                if not beginning:direct=inverse(direct,interest,'add',True,d['principal_model'])
                if d['payment_store']:direct=B.pullback(direct)
                require(E.coordinates(direct)==E.coordinates(required),'direct pairwise singleton has a different Ext80 domain')
            require(d['required']==(required.json()if required else None),'first-singleton inverse differs')
            bounds=E.coordinates(required)
            if bounds:
                low,high=bounds;require(d['count']==high-low+1,'initial domain count differs')
                require(tagged_value(d['lower'])==E.value(low)and tagged_value(d['upper'])==E.value(high),'initial domain omits/adds lattice endpoints')
            else:require(d['count']==0 and d['lower']is None and d['upper']is None,'empty preimage is hidden')
            key=(li,d['payment_store'],d['interest_model'],d['principal_model'],d['sum_model'])
            if interest_sum:key+=d['sum_order'],
            require(key not in free_domains,'duplicate initial domain');free_domains[key]=bounds
        factors=[range(60),(False,True),MODELS,MODELS,MODELS]
        if interest_sum:factors.append(range(3))
        keys=set(itertools.product(*factors))
        require(set(free_domains)==keys,'incomplete or widened first-shell domain catalogue')
    def domain(pid,li):
        if not free_domains:return domains[li]
        p=programs[pid];key=(li,p['payment_store'],p['interest_model'],p['principal_model'],p['sum_model'])
        if interest_sum:key+=p['sum_order'],
        bounds=free_domains[key]
        return []if bounds is None else range(bounds[0],bounds[1]+1)
    def payment(pid,li,index):
        value=domain(pid,li)[index]
        return E.value(value)if free_domains else value
    programs=list(records(directory/'profiles.jsonl'));declared=catalogue()
    require(programs==declared and len(programs)==7500,'profile grammar/order differs or incomplete run')
    graphs={};counts=Counter()
    for item in records(directory/'graphs.jsonl'):
        key=tuple(item['graph_id']);require(key not in graphs,'duplicate graph ID')
        pid,loan_index,row_index=key;require(pid in range(7500)and row_index in loans[loan_index]['row_indices'],'invalid graph identity')
        raw=rows[row_index];args,n,start,end,beginning=coordinates(raw)
        # coordinates()'s third ordinary leaf is -PV; replace it with the
        # symbolic payment at the complete cell's two extreme lattice values.
        endpoint_args=[[args[0],args[1],p,Q(0),Q(0)]for p in (payment(pid,loan_index,0),payment(pid,loan_index,-1))]
        g=item['graph'];endpoint_values=[graph_trace(g,a)for a in endpoint_args]
        shapes=[];depends=[]
        for j,instruction in enumerate(g['nodes']):
            kind,operands=next(iter(instruction['op'].items()));model=instruction['model'];require(model in MODELS,'unknown model')
            if kind=='arg':
                require(model=='strict'and operands in range(3),'undeclared argument')
                shape=node(kind,model,operands);dependent=operands==2
            elif kind=='const':
                value=binary64('0x'+format(operands['bits'],'016x'))
                require(model=='strict'and (operands['bits']in (0,0x3ff0000000000000)or interest_sum and value.denominator==1 and 1<=value<=512),'undeclared constant')
                shape=node(kind,model,operands['bits']);dependent=False
            else:
                require(kind in ('add','mul')and len(operands)==2 and all(0<=k<j for k in operands),'invalid or nonmonotone instruction')
                a,b=operands;shape=node(kind,model,shapes[a],shapes[b]);dependent=depends[a]or depends[b]
                if kind=='mul'and dependent:
                    require(depends[a]!=depends[b],'product of two variable subexpressions is outside the monotonicity proof')
                    fixed=b if depends[a]else a
                    require(endpoint_values[0][fixed]==endpoint_values[1][fixed]>0,'variable multiplied by a nonpositive/nonconstant factor')
            if not dependent:
                require(endpoint_values[0][j]==endpoint_values[1][j],'purported constant depends on payment')
            require(endpoint_values[0][j]<=endpoint_values[1][j],'endpoint order violates monotonicity')
            shapes.append(shape);depends.append(dependent)
        require(shapes[g['output']]==expected_shape(programs[pid],start,end,beginning,interest_sum),'balance expression differs from declared profile')
        graphs[key]=g;counts['monotone_graphs_checked']+=1
    traces={}
    for t in records(directory/'boundary-traces.jsonl'):
        key=tuple(t['graph_id']);require(key in graphs and t['trace_id']not in traces,'unknown/duplicate trace')
        pid,loan_index,row_index=key;index=t['latent_index'];require(index in range(len(domain(pid,loan_index))),'payment index out of range')
        a,_,_,_,_=coordinates(rows[row_index]);args=[a[0],a[1],payment(pid,loan_index,index),Q(0),Q(0)]
        require(list(map(tagged_value,t['inputs']))==args,'typed boundary inputs changed')
        g=graphs[key];expected=graph_trace(g,args);require(len(expected)==len(t['values']),'incomplete boundary trace')
        for instruction,encoded,value in zip(g['nodes'],t['values'],expected):
            kind,operands=next(iter(instruction['op'].items()))
            ext=(kind=='arg'and operands==2)or(kind not in ('arg','const')and instruction['model'].startswith('resident_'))
            require(encoded.startswith('ext80-le:'if ext else'binary64:'),'wrong typed intermediate representation')
            require(tagged_value(encoded)==value,'boundary arithmetic differs');counts['nodes_checked']+=1
        actual=binary64(t['actual']);require(actual<0 and B.round(expected[g['output']])==actual,'boundary publication differs')
        traces[t['trace_id']]=(key,index,t['actual']);counts['boundary_traces_checked']+=1
    constraints=list(records(directory/'constraints.jsonl'));position=0;survivors=[];histogram=Counter();used_traces=set();used_graphs=set()
    for p in programs:
        killed=False;feasible=[]
        for loan_index,loan in enumerate(loans):
            low,high=0,len(domain(p['id'],loan_index))
            if high==0:
                killed=True;histogram[loan['row_indices'][0]]+=1;counts['empty_initial_domains']+=1;break
            for row_index in loan['row_indices']:
                require(position<len(constraints),'missing constraint');c=constraints[position];position+=1
                key=(p['id'],loan_index,row_index);require(c['profile_id']==p['id']and c['loan_index']==loan_index and c['row_index']==row_index and tuple(c['graph_id'])==key,'constraint identity/order differs')
                require(c['before']==[low,high]and c['expected']==rows[row_index]['expected'],'constraint starts from wrong interval/observation')
                left,right=c['after'];require(low<=left<=right<=high,'invalid lattice interval cut')
                samples={}
                for b in c['boundaries']:
                    tid=b['trace_id'];require(traces[tid]==(key,b['latent_index'],b['actual']),'boundary references differ')
                    require(b['latent_index']not in samples,'duplicate boundary sample');samples[b['latent_index']]=binary64(b['actual']);used_traces.add(tid)
                target=binary64(c['expected'])
                if left>low:require(samples[left-1]<target,'left exclusion boundary failed')
                if left<high:require(samples[left]>=target,'left inclusion boundary failed')
                if right>left:require(samples[right-1]<=target,'right inclusion boundary failed')
                if right<high:require(samples[right]>target,'right exclusion boundary failed')
                # A second, direct inverse proof checks the present first-shell
                # exclusions independently of the forward bisection certificate.
                if not free_domains and rows[row_index]['args'][3:5]==['0x3ff0000000000000']*2 and binary64(rows[row_index]['args'][5])==0:
                    a,_,_,_,_=coordinates(rows[row_index]);required=B.pullback(Span(target,target))
                    if p['sum_order']!=2:required=inverse(required,Q(0),'add',True,p['sum_model'])
                    interest=operation(a[0],a[1],'mul',p['interest_model'])
                    required=inverse(required,interest,'add',True,p['principal_model'])
                    if p['payment_store']:required=B.pullback(required)
                    published=binary64(loan['payment']['expected']);publication=B.pullback(Span(published,published))
                    common=publication.meet(required)
                    if left==right:require(common is None or E.coordinates(common)is None,'direct first-shell inverse admits a payment')
                    counts['direct_first_shell_inverse_checks']+=1
                low,high=left,right;used_graphs.add(key);counts['interval_cuts_checked']+=1
                if low==high:killed=True;histogram[row_index]+=1;break
            if killed:break
            feasible.append(dict(loan_index=loan_index,latent_interval=[low,high]))
        if not killed:survivors.append(dict(profile_id=p['id'],loans=feasible))
    require(position==len(constraints)and used_graphs==set(graphs)and used_traces==set(traces),'unused/missing proof material')
    actual=list(records(directory/'survivors.jsonl'))
    require([{k:r[k]for k in ('profile_id','loans')}for r in actual]==survivors,'survivor domains differ')
    summary=json.loads((directory/'summary.json').read_text())
    require(summary['complete']and summary['completed_profiles']==7500 and summary['survivors']==len(survivors),'incomplete or wrong result summary')
    require(summary['first_kill_histogram']=={str(k):v for k,v in sorted(histogram.items())},'first-kill counts differ')
    require(digest(Path(__file__))==SOURCE,'checker changed during execution')  # ubs:ignore — public source identity, no secret comparison
    require(digest(Path(__file__).parent/'cumulative.py')==HELPER_SOURCE,'cumulative helper changed during execution')  # ubs:ignore — public source identity, no secret comparison
    write_result(output,dict(complete=True,profiles=7500,survivors=len(survivors),original_available_rows=540,
        loan_publication_domains=60,**counts,first_kill_histogram=histogram,seconds=time.monotonic()-started,
        balance_checker_sha256=SOURCE,cumulative_helper_sha256=HELPER_SOURCE,inputs={str(p):digest(p)for p in (Path(trace),Path(structure),directory/'profiles.jsonl',directory/'constraints.jsonl',directory/'graphs.jsonl',directory/'boundary-traces.jsonl',directory/'loans.jsonl')},
        initial_domain='complete first-singleton inverse; no external PMT constraint'if free_domains else'external PMT publication',
        interest_sum=interest_sum,
        first_singleton_domains_independently_checked=len(free_domains),
        initial_domains_sha256=digest(free_file)if free_domains else None,
        claim='Every initial payment lattice is complete. Independently reconstructed positive-coefficient DAGs establish monotonicity; exact full-width corner traces prove every lattice cut. No latent payment is fitted, selected or installed.',
        limitations=('Necessary consumers with one arbitrary Ext80 payment shared across the original intervals; no external PMT constraint. First-singleton inverse domains and all later cuts are independently verified. 'if free_domains else
            'Necessary consumers assume CUMPRINC internal payment publishes to external PMT. All first-singleton failures say nothing about accumulation. ')+
            'No complete payment producer, CUMPRINC kernel or global backend equivalence identified. Historical read-back audits unavailable.'))


if __name__=='__main__':check(*sys.argv[1:])
