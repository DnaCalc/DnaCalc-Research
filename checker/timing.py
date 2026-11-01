"""Independent exact checks of PMT timing and shared-helper programs.

No racer, candidate evaluator or research interval imports. The six operations
below follow the declared annuity identity, using the independent rational
rounding/preimage implementation. Original witnesses come from the one loader.
"""
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path
import json
import sys
import time

from check import (require, records, digest, validate_row, source_hashes, product_preimage,
                   tagged_value, graph_shape, expected_shape, graph_trace)
from exact import BINARY64 as B, EXTENDED64 as E, Span, binary64, inverse, operation, round_model, two

MODELS = ('strict','stored_pc64','resident_pc64')
LOADED_SOURCE = {**source_hashes(), 'timing.py':digest(Path(__file__))}


def preimage(row, models, factor_operation=None, due_form=None):
    rate,pv,_,timing,_ = map(binary64,row['leaves'])
    require(timing in (0,1), 'noncanonical timing')
    if factor_operation is not None:
        require(factor_operation in ('mul','div'),'unknown factor operation')
        target=binary64(row['raw']['expected'])
        require(target and pv,'zero target/balance outside factor scope')
        return inverse(B.pullback(Span(target,target)),pv,factor_operation,False,models[0])
    rt = rate*timing  # exact multiplication by0 or1, no arithmetic choice
    if due_form is not None:
        if due_form==0:
            tf=operation(Q(1),rt,'add',models[0])
            product=operation(rate,rt,'mul',models[1])
            correction=operation(product,tf,'div',models[2])
            factor=operation(rate,-correction,'add',models[3])
        elif due_form==1:
            tf=operation(Q(1),rt,'add',models[0])
            ratio=operation(rt,tf,'div',models[1])
            complement=operation(Q(1),-ratio,'add',models[2])
            factor=operation(rate,complement,'mul',models[3])
        else:
            require(due_form==2,'unknown due coefficient identity')
            reciprocal=operation(Q(1),rate,'div',models[0])
            denominator=operation(reciprocal,timing,'add',models[1])
            factor=operation(Q(1),denominator,'div',models[2])
        target=binary64(row['raw']['expected'])
        require(target and pv and factor,'zero target/balance/coefficient outside scope')
        result=inverse(B.pullback(Span(target,target)),factor,'mul',True,models[-1])
        return inverse(result,pv,'div',False,models[-2])
    tf = operation(Q(1),rt,'add',models[0])
    ratio = operation(rt,tf,'div',models[1])
    factor = operation(Q(1),-ratio,'add',models[2])
    target = binary64(row['raw']['expected'])
    require(target and pv and factor, 'zero constant outside this inverse scope')
    result = B.pullback(Span(target,target))
    result = inverse(result,rate,'mul',True,models[5])
    result = inverse(result,factor,'mul',True,models[4])
    return inverse(result,pv,'div',False,models[3])


def product_shells(scoped, route_path, output, discount=False):
    """Bind the existing final-product exclusion to the current route audit.

    Quantify over every finite binary64 quotient. No helper or quotient is
    selected, and a nonempty interval supplies no function match.
    """
    start=time.monotonic()
    route_rows=json.loads(route_path.read_text())
    routes={r['id']:r for r in route_rows}
    require(len(routes)==len(route_rows),'duplicate audited route ID')
    counts=Counter();proofs=[];seen=set()
    for row in records(scoped):
        validate_row(row);raw=row['raw'];route=routes[raw['id']]
        require(raw['id']not in seen,'duplicate original witness ID');seen.add(raw['id'])
        require(route['args']==raw['args']and route['expected']==raw['expected'],'route/witness identity differs')
        require(route['authored_echoes']==raw['args'],'application changed authored inputs')
        require(len({route[k]for k in ('workbook','registered','direct_body','saved_reopened')})==1,'application routes differ')
        r=abs(binary64(row['leaves'][0]));y=abs(binary64(raw['expected']))
        require(r and y>=B.value(B.half),'normal nonzero target/rate required')
        counts['rows']+=1
        counts['application_matches']+=route['workbook']==raw['expected']
        multipliers=[('',r)]
        if discount:
            timing=binary64(row['leaves'][3]);require(timing in (0,1),'noncanonical timing')
            multipliers=[]
            for factor_model in MODELS:
                factor=operation(Q(1),r*timing,'add',factor_model)
                for division_model in MODELS:
                    d=operation(r,factor,'div',division_model)
                    multipliers.append((factor_model+'/'+division_model+'/',d))
        shells={}
        for label,multiplier in multipliers:
            for model,span in (('strict',B.pullback(Span(y,y))),('stored_pc64',product_preimage(y))):
                effective=B.round(multiplier)if model=='strict'else multiplier
                interval=span.linear(1/effective);exists=B.coordinates(interval)is not None
                counts[label+model+'_unreachable']+=not exists
                neighbors=None
                if not exists:
                    k=B.floor_coordinate(interval.low)
                    low,high=B.value(k),B.value(k+1)
                    require(low<=interval.low<=interval.high<=high,'incorrect gap bracket')
                    require(not interval.contains(low)and not interval.contains(high),'a quotient reaches the target')
                    require(round_model(low*effective,model)!=y and round_model(high*effective,model)!=y,'neighbor replay contradicts gap')
                    neighbors=[[str(q.numerator),str(q.denominator)]for q in (low,high)]
                shells[label+model]=dict(finite_binary64_quotient_exists=exists,quotient_magnitude_interval=interval.json(),adjacent_magnitudes=neighbors)
                if discount:shells[label+model]['effective_multiplier']=[str(effective.numerator),str(effective.denominator)]
        proofs.append(dict(witness=raw,route=route,shells=shells))
    require(LOADED_SOURCE=={**source_hashes(),'timing.py':digest(Path(__file__))},'checker sources changed')
    receipt=dict(complete=True,counts=counts,proofs=proofs,seconds=time.monotonic()-start,
        inputs={str(p):digest(p)for p in (scoped,route_path)},checker_sha256=LOADED_SOURCE,
        claim='Every excluded normal target lies outside the final RN53(q*r) and RN53(RN64(q*r)) shells for every finite binary64 q. The original current-build observations and unchanged application routes are bound explicitly.',
        limits='This identifies a necessary structural difference at or before the final product. It does not identify an executable quotient/helper producer, prove a particular Excel precision, or score reachability as a match.')
    if discount:
        receipt.update(profiles=18,all_rows_reachable_profiles=[key[:-12]for key,count in counts.items()if key.endswith('_unreachable')and not count],
            claim='Exact product preimages quantify over every finite binary64 q for the18 declared final actuarial-discount-rate shells. The multiplier is calculated solely from rate/timing and fixed arithmetic models, never from an answer.',
            limits='A nonempty shell is only a necessary reachability condition, not an executable quotient/helper or a matched function. No fitted or selected intermediate, product change or claim about other multipliers/precision programs.')
    with output.open('x')as f:json.dump(receipt,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in receipt.items()if k!='proofs'}))


def helper_component(h, component, n, log, elementary):
    """Verify helper arithmetic; exp/log values remain explicit hypotheses."""
    t=-operation(n,log,'mul',h['t'])
    require(tagged_value(component['t'])==t,'helper exponent arithmetic differs')
    exp=tagged_value(component['exp_raw'])
    key=('exp',t)
    require(key not in elementary or elementary[key]==exp,'same elementary input has different values')
    elementary[key]=exp
    u=B.round(exp)if h['exp_store']else exp
    require(tagged_value(component['u'])==u,'helper exponential store differs')
    if u==1:
        em=t
        require(component['log_raw']is None,'unused logarithm in unit-exponential branch')
    else:
        d=operation(u,Q(-1),'add',h['difference'])
        if not h['correction']or abs(t)>=1:
            em=d
            require(component['log_raw']is None,'unused logarithm outside Kahan branch')
        else:
            log_raw=tagged_value(component['log_raw']);key=('log',u)
            require(key not in elementary or elementary[key]==log_raw,'same logarithm input has different values')
            elementary[key]=log_raw
            log_u=B.round(log_raw)if h['log_store']else log_raw
            a,b=h['arithmetic']
            if h['association']==0:em=operation(operation(d,t,'mul',a),log_u,'div',b)
            elif h['association']==1:em=operation(d,operation(t,log_u,'div',a),'mul',b)
            else:
                require(h['association']==2,'unknown helper association')
                em=operation(operation(d,log_u,'div',a),t,'mul',b)
    em=B.round(em)if h['em_store']else em
    require(tagged_value(component['em'])==em,'helper correction arithmetic differs')
    return em


def series_log1p(rate):
    """DLMF4.6.E6 in resident PC64, with an unchanged-sum stop."""
    require(0<rate<=Q(1,4),'outside declared log-provider domain')
    op=lambda a,b,kind:operation(a,b,kind,'resident_pc64')
    z=op(rate,op(Q(2),rate,'add'),'div');square=op(z,z,'mul')
    term=total=z
    for k in range(3,1024,2):
        term=op(term,square,'mul');new=op(total,op(term,Q(k),'div'),'add')
        if new==total:return op(Q(2),total,'mul')
        total=new
    raise ValueError('declared log series did not reach an unchanged sum')


def direct_component(h,component,n,log,elementary):
    """Existing negative-exponent F2XM1 helper, with a conditional raw leaf."""
    values=list(map(tagged_value,component['nodes']))
    require(len(values)==10,'incomplete direct helper trace')
    mt,mx,mf,me=h['models']
    t=operation(-n,log,'mul',mt)
    rom=Q(0xb8aa3b295c17f0bc,1<<63)
    x=operation(t,rom,'mul',mx)
    lower=x.numerator//x.denominator;fraction=x-lower
    k=Q(lower+int(fraction>Q(1,2)or(fraction==Q(1,2)and lower%2)))
    reduced=operation(x,-k,'add',mf)
    require(abs(reduced)<=Q(1,2),'F2XM1 argument outside declared reduced interval')
    require(values[:4]==[t,x,k,reduced],'direct range-reduction arithmetic differs')
    raw=values[4];key=('f2xm1',reduced)
    require(key not in elementary or elementary[key]==raw,'same F2XM1 input differs')
    elementary[key]=raw
    f=B.round(raw)if h['f2xm1_store']else raw
    # This verifier aborts on scales outside the finite-Ext80 exponent range;
    # no omitted or unsupported scale may supply a successful certificate.
    require(E.emin-E.precision<=k<=E.emax+1,'scale exponent outside finite checker scope')
    power=two(int(k));scaled=E.round(f*power);two_k=E.round(power)
    integer_part=operation(two_k,Q(-1),'add','resident_pc64')
    em=operation(scaled,integer_part,'add',me)
    require(values[5:]==[f,scaled,two_k,integer_part,em],'direct scaling/reconstruction arithmetic differs')
    require(tagged_value(component['em'])==em,'published direct helper differs')
    return em


def lossless_pc53_nodes(graph,encoded_nodes):
    """Check the local premise for a Stored/Resident-PC53 labeling change."""
    require(len(encoded_nodes)==len(graph['nodes']),'equivalence trace length differs')
    count=0;values={}
    def value(index):
        if index not in values:values[index]=tagged_value(encoded_nodes[index])
        return values[index]
    for i,node in enumerate(graph['nodes']):
        if node['model']!='stored_pc53':continue
        require(encoded_nodes[i].startswith('binary64:'),'PC53 store is not recorded as binary64')
        stored=value(i)
        require(stored!=0,'zero sign outside this bounded equivalence proof')
        kind,args=next(iter(node['op'].items()))
        require(kind in ('add','mul','div')and len(args)==2,'PC53 annotation outside binary arithmetic')
        resident=operation(*(value(a)for a in args),kind,'resident_pc53')
        require(resident==stored,'PC53 spill loses a resident value')
        count+=1
    return count


def pc53_equivalence(basis_path,directory,archive_path,output):
    """Compose a checked four-model proof with exact PC53 spill equalities."""
    import hashlib
    import tarfile
    start=time.monotonic();basis=json.loads(basis_path.read_text())
    require(basis['complete']and basis['eliminated']==basis['programs']and basis['survivors']==0,'basis is not a complete exclusion')
    require(basis['consumer_models']==list(MODELS+('stored_pc53',)),'wrong canonical model family')
    for name in ('exact.py','check.py'):
        require(basis['checker_sha256'][name]==source_hashes()[name],'basis arithmetic/checker library differs')
    with tarfile.open(archive_path,'r:gz')as archive:
        member=archive.extractfile('research/checker/timing.py')
        require(member is not None,'missing basis checker source')
        require(hashlib.sha256(member.read()).hexdigest()==basis['checker_sha256']['timing.py'],'basis checker archive differs')
    for path,wanted in basis['inputs'].items():
        require(digest(Path(path))==wanted,'basis proof input changed') # ubs:ignore — public file identities, not secrets or authentication
    for name in ('consumers.json','first-kills.jsonl','helper-profiles.json'):
        require(str(directory/name)in basis['inputs'],'basis does not bind this equivalence input')
    cs=json.loads((directory/'consumers.json').read_text())
    helper_count=len(json.loads((directory/'helper-profiles.json').read_text()))
    require(len(cs)==2*4**6 and helper_count==(810 if basis['direct']else 9900),'wrong canonical cross dimensions')
    require(basis['programs']==len(cs)*helper_count,'basis profile count differs')
    seen=bytearray(basis['programs']);canonical=expanded=groups=equalities=0
    for kill in records(directory/'first-kills.jsonl'):
        cid=kill['consumer'];require(type(cid)is int and 0<=cid<len(cs),'consumer ID outside equivalence family')
        c=cs[cid];k=c['models'].count('stored_pc53')
        require(lossless_pc53_nodes(c['graph'],kill['nodes'])==k,'PC53 site count differs')
        for hid in kill['helpers']:
            require(type(hid)is int and 0<=hid<helper_count,'helper ID outside equivalence family')
            identity=cid*helper_count+hid;require(not seen[identity],'duplicate canonical first kill')
            seen[identity]=1;canonical+=1
        expanded+=len(kill['helpers'])*(1<<k);groups+=1;equalities+=k
    require(canonical==basis['programs']and all(seen),'canonical coverage incomplete')
    require(expanded==2*5**6*helper_count,'five-model weighted coverage incomplete')
    require(LOADED_SOURCE=={**source_hashes(),'timing.py':digest(Path(__file__))},'equivalence checker sources changed')
    receipt=dict(complete=True,canonical_profiles=canonical,expanded_profiles=expanded,
        canonical_consumers=len(cs),expanded_consumers=2*5**6,helper_profiles=helper_count,
        checked_groups=groups,checked_nonzero_lossless_pc53_nodes=equalities,
        seconds=time.monotonic()-start,checker_sha256=LOADED_SOURCE,
        inputs={str(p):digest(p)for p in (basis_path,archive_path)},basis_proof_inputs=basis['inputs'],
        claim='Each completed canonical first kill also contradicts every labeling obtained by replacing any subset of its StoredPc53 sites with ResidentPc53. Exact nonzero spill equality propagates through this fixed DAG. All five-model assignments have one canonical image and weighted coverage is complete.',
        limits='A superset of the canonical four-model exclusions, never an additional count to sum. Consumer-only on the certified counterexamples and unchanged helper/log hypotheses. No global spill equivalence; finite nonzero equality guards exclude signed-zero and underflow/extended-exponent differences. No fitted value, identified Excel kernel or product gain.')
    with output.open('x')as f:json.dump(receipt,f,indent=2);f.write('\n')
    print(json.dumps(receipt))


def canonical_kahan_helpers():
    import itertools
    expected=[]
    for t,exp_store,difference,correction in itertools.product(MODELS,(True,False),MODELS,(True,False)):
        for log_store in ((True,False)if correction else(True,)):
            for association in range(3 if correction else 1):
                for a,b in itertools.product(MODELS if correction else('strict',),repeat=2):
                    for em_store in (True,False):
                        expected.append(dict(id=len(expected),t=t,exp_store=exp_store,difference=difference,
                            log_store=log_store,association=association,arithmetic=[a,b],em_store=em_store,correction=correction))
    return expected


def discount_helpers(scoped,directory,output):
    """Check the correlated actuarial log and resident-quotient cross."""
    start=time.monotonic()
    hs=json.loads((directory/'helpers.json').read_text())
    require(hs==canonical_kahan_helpers(),'discount helper grammar differs')
    profiles=[dict(id=dp*1980+h,log_provider=dp,helper=h)for dp in range(9)for h in range(1980)]
    require(json.loads((directory/'helper-profiles.json').read_text())==profiles,'discount helper profiles differ')
    cs=json.loads((directory/'consumers.json').read_text())
    require(len(cs)==81,'discount consumer count differs')
    rate_programs=[]
    ops=[{'arg':i}for i in range(5)]+[{'const':{'bits':0x3ff0000000000000}},
         {'add':[5,4]},{'mul':[2,6]},{'add':[1,7]},{'mul':[0,3]},
         {'add':[5,9]},{'div':[8,4]},{'div':[0,10]},{'mul':[11,12]}]
    for cid,c in enumerate(cs):
        code=cid;models=[]
        for _ in range(4):models.append(MODELS[code%3]);code//=3
        require(c['id']==cid and c['numerator']==0 and c['combine']==5 and c['models']==models,
                'discount consumer identity differs')
        g=c['graph'];require(len(g['nodes'])==14 and g['output']==13,'discount graph layout differs')
        require([n['op']for n in g['nodes']]==ops and graph_shape(g)==expected_shape(0,5),'discount graph shape differs')
        for i,node in enumerate(g['nodes']):
            require(node['model']==(models[i-10]if i>=10 else'strict'),'discount site annotation differs')
        rate_programs.append(MODELS.index(models[0])*3+MODELS.index(models[2]))
    require(json.loads((directory/'rate-programs.json').read_text())==rate_programs,'rate/consumer correlation differs')
    originals={r['trace_index']:r for r in records(scoped)}
    rows={};elementary={};components=0
    for entry in records(directory/'helper-values.jsonl'):
        row=originals[entry['trace_index']];validate_row(row)
        r,pv,fv,timing,_=map(binary64,row['leaves']);n=binary64(row['raw']['args'][1])
        require(0<r<=1 and pv and fv==0 and timing==1 and n>0,'outside discount producer scope')
        require(row['raw']==entry['witness']and entry['order']not in rows,'changed/duplicate discount helper row')
        logs=list(map(tagged_value,entry['log_inputs']));require(len(logs)==len(entry['rate_values'])==9,'discount rate coverage differs')
        for dp,value in enumerate(entry['rate_values']):
            tf=operation(Q(1),r,'add',MODELS[dp//3]);d=operation(r,tf,'div',MODELS[dp%3])
            argument=-B.round(d);raw_log=tagged_value(value['raw_log'])
            require(tagged_value(value['tf'])==tf and tagged_value(value['d'])==d,'discount rate arithmetic differs')
            require(value['argument'].startswith('binary64:')and tagged_value(value['argument'])==argument,
                    'discount logarithm argument changed')
            require(value['raw_log'].startswith('binary64:')and raw_log<0 and logs[dp]==-raw_log,
                    'discount logarithm result/sign differs')
            key=('cleared_log1p',argument)
            require(key not in elementary or elementary[key]==raw_log,'repeated log1p input has different values')
            elementary[key]=raw_log
        require(len(entry['components'])==len(entry['values'])==len(profiles),'discount helper row incomplete')
        for p,comp,value in zip(profiles,entry['components'],entry['values']):
            em=helper_component(hs[p['helper']],comp,n,logs[p['log_provider']],elementary)
            require(comp['em']==value and B.round(em)is not None,'discount helper publication differs')
            for model in MODELS:require(B.round(operation(Q(1),em,'add',model))is not None,'neutral numerator overflows')
            components+=1
        rows[entry['order']]=(row,entry)
    seen=set();groups=nodes=0
    for kill in records(directory/'first-kills.jsonl'):
        row,entry=rows[kill['order']];cid=kill['consumer']
        require(type(cid)is int and 0<=cid<81,'discount consumer ID outside grammar')
        require(kill['trace_index']==row['trace_index']and kill['witness']==row['raw'],'discount witness changed')
        require(kill['expected']==row['raw']['expected']and kill['got']!=kill['expected'],'discount kill is not a contradiction')
        require(kill['helpers'],'empty discount helper group')
        for hid in kill['helpers']:
            require(type(hid)is int and 0<=hid<len(profiles)and (cid,hid)not in seen,'duplicate/invalid discount profile')
            require(profiles[hid]['log_provider']==rate_programs[cid],'uncorrelated log/final discount rate')
            require(entry['values'][hid]==kill['helper_value'],'grouped discount helpers differ')
            seen.add((cid,hid))
        args=list(map(binary64,row['leaves']));args[4]=tagged_value(kill['helper_value'])
        graph=cs[cid]['graph'];trace=graph_trace(graph,args)
        require(list(map(tagged_value,kill['nodes']))==trace,'discount consumer trace differs')
        require(trace[8]==args[1]and trace[12]==tagged_value(entry['rate_values'][rate_programs[cid]]['d']),
                'neutral numerator or correlated multiplier differs')
        require(B.round(trace[13])==binary64(kill['got']),'discount consumer publication differs')
        groups+=1;nodes+=len(trace)
    require(len(seen)==160380,'discount compatible profile coverage incomplete')
    require(json.loads((directory/'survivors.json').read_text())==[],'discount survivors remain')
    require(LOADED_SOURCE=={**source_hashes(),'timing.py':digest(Path(__file__))},'discount checker source changed')
    inputs=[scoped,*[directory/n for n in ('helpers.json','helper-profiles.json','consumers.json','rate-programs.json','helper-values.jsonl','first-kills.jsonl','survivors.json','summary.json')]]
    receipt=dict(complete=True,programs=160380,eliminated=len(seen),survivors=0,
        checked_helper_program_rows=components,checked_consumer_groups=groups,checked_consumer_nodes=nodes,
        original_rows=len(rows),elementary_constraints=len(elementary),seconds=time.monotonic()-start,
        checker_sha256=LOADED_SOURCE,inputs={str(p):digest(p)for p in inputs},
        claim='All160380 compatible actuarial-discount producer/consumer programs contradict original R10 observations. The same rate-factor/division operations determine the cleared log1p argument and final multiplier. Every helper arithmetic/store and grouped consumer trace is checked exactly.',
        limits='FV0,timing1,0<rate<=1,positive period/normal targets; existing R7 helpers and three arithmetic models only. Cleared log1p/exp/log values are recorded conditional leaves with repeated-operand consistency, not identified Excel helpers. No fitted value, selected survivor or product gain.')
    with output.open('x')as sink:json.dump(receipt,sink,indent=2);sink.write('\n')
    print(json.dumps(receipt))


def timing_helpers(scoped,directory,output,multiple_logs=False,direct=False,pc53=False):
    """Check grouped first kills of the complete finite helper/consumer cross."""
    import itertools
    start=time.monotonic()
    hs=json.loads((directory/'helpers.json').read_text())
    cs=json.loads((directory/'consumers.json').read_text())
    expected=canonical_kahan_helpers()
    if direct:
        expected=[]
        for hid in range(162):
            code=hid//2;models=[]
            for _ in range(4):models.append(MODELS[code%3]);code//=3
            expected.append(dict(id=hid,models=models,f2xm1_store=bool(hid%2)))
        constants=json.loads((directory/'constants.json').read_text())
        require(constants=={'ROM_log2e':'ext80-le:bcf0175c293baab8ff3f'},'undeclared direct-helper constant')
    require(hs==expected,'helper Cartesian coverage differs')
    helper_count=162 if direct else 1980
    profiles=[dict(id=lp*helper_count+h,log_provider=lp,helper=h)
        for lp in range(5 if multiple_logs else 1)for h in range(helper_count)]
    profile_path=directory/'helper-profiles.json'
    if profile_path.exists():require(json.loads(profile_path.read_text())==profiles,'log/helper Cartesian coverage differs')
    else:require(not multiple_logs,'missing declared helper profiles')
    consumer_models=MODELS+('stored_pc53',)if pc53 else MODELS
    model_count=len(consumer_models)**6
    require(len(cs)==2*model_count,'consumer coverage differs')
    for cid,c in enumerate(cs):
        form=cid//model_count;code=cid%model_count;models=[]
        for _ in range(6):models.append(consumer_models[code%len(consumer_models)]);code//=len(consumer_models)
        require(c['id']==cid and c['numerator']==0 and c['combine']==9+form and c['models']==models,'consumer identity differs')
        graph=c['graph'];require(graph_shape(graph,True)==expected_shape(0,9+form),'consumer shape differs')
        ops=[{'arg':i}for i in range(5)]+[{'const':{'bits':0x3ff0000000000000}},
            {'add':[5,4]},{'mul':[2,6]},{'add':[1,7]},{'mul':[0,3]},
            {'add':[5,9]},{'div':[8,4]}]
        ops+=([{'div':[9,10]},{'mul':[11,12]},{'const':{'bits':0xbff0000000000000}},
               {'mul':[14,13]},{'add':[11,15]},{'mul':[16,0]}]if form==0 else
              [{'mul':[11,0]},{'div':[9,10]},{'mul':[12,13]},{'const':{'bits':0xbff0000000000000}},
               {'mul':[15,14]},{'add':[12,16]}])
        require([n['op']for n in graph['nodes']]==ops,'consumer DAG or operation order differs')
        sites=(10,11,12,13,16,17)if form==0 else(10,11,12,13,14,17)
        negation=15 if form==0 else 16
        require(len(graph['nodes'])==18 and graph['output']==17,'consumer node layout differs')
        for i,node in enumerate(graph['nodes']):
            if i in sites:model=models[sites.index(i)]
            else:model='resident_pc64'if i==negation else'strict'
            require(node['model']==model,'neutral arithmetic or exact sign change differs')
    originals={r['trace_index']:r for r in records(scoped)}
    rows={};elementary={};component_count=0
    for entry in records(directory/'helper-values.jsonl'):
        row=originals[entry['trace_index']];validate_row(row)
        require(entry['witness']==row['raw'],'original helper input changed')
        require(entry['order']not in rows,'duplicate helper row')
        rate,pv,fv,timing,_=map(binary64,row['leaves'])
        require(rate>0 and binary64(row['raw']['args'][1])>0 and fv==0 and abs(pv)>=B.value(B.half),'outside FV0 projection')
        if multiple_logs:require(rate<=Q(1,4),'rate outside hardware/series provider domain')
        encoded_logs=entry.get('log_inputs',[f'binary64:{row["log1p"]:016x}'])
        logs=list(map(tagged_value,encoded_logs))
        require(len(logs)==(5 if multiple_logs else 1)and encoded_logs[0]==f'binary64:{row["log1p"]:016x}','frozen log input differs')
        if multiple_logs:
            require(logs[2]==B.round(logs[1]),'stored hardware log differs')
            require(logs[3]==series_log1p(rate)and logs[4]==B.round(logs[3]),'resident/stored DLMF series differs')
            key=('hardware_log1p',rate)
            require(key not in elementary or elementary[key]==logs[1],'same hardware log input differs')
            elementary[key]=logs[1]
        require(len(entry['values'])==len(entry['components'])==len(profiles),'incomplete helper row')
        for p,comp,value in zip(profiles,entry['components'],entry['values']):
            verify=direct_component if direct else helper_component
            em=verify(hs[p['helper']],comp,binary64(row['raw']['args'][1]),logs[p['log_provider']],elementary)
            # Every neutral numerator/rt annotation is exact in this scope.
            # Refuse an extended helper whose strict conversion or1+em could
            # overflow;0*Inf would invalidate the FV0 reduction.
            require(B.round(em)is not None,'helper conversion overflows the neutral prefix')
            for model in MODELS:
                require(B.round(operation(Q(1),em,'add',model))is not None,'neutral1+em overflows')
            require(comp['em']==value,'published helper differs from its component audit');component_count+=1
        rows[entry['order']]=(row,entry)
    total=len(profiles)*len(cs);seen=bytearray(total);groups=0;node_count=0;eliminated=0
    for kill in records(directory/'first-kills.jsonl'):
        row,entry=rows[kill['order']];cid=kill['consumer']
        require(type(cid)is int and 0<=cid<len(cs),'consumer ID outside canonical grammar')
        graph=cs[cid]['graph']
        require(kill['trace_index']==row['trace_index']and kill['witness']==row['raw'],'first-kill witness changed')
        require(kill['expected']==row['raw']['expected']and kill['got']!=kill['expected'],'not a contradiction')
        require(kill['helpers'],'empty helper group')
        for hid in kill['helpers']:
            require(type(hid)is int and 0<=hid<len(profiles)and not seen[cid*len(profiles)+hid],'duplicate/out-of-range helper program')
            require(entry['values'][hid]==kill['helper_value'],'grouped helpers differ')
            seen[cid*len(profiles)+hid]=1;eliminated+=1
        args=list(map(binary64,row['leaves']));args[4]=tagged_value(kill['helper_value'])
        actual=list(map(tagged_value,kill['nodes']));computed=graph_trace(graph,args)
        require(actual==computed,'full-width consumer trace differs')
        require(computed[8]==args[1]and computed[9]==args[0]*args[3],'neutral numerator/timing reduction differs')
        require(B.round(computed[graph['output']])==binary64(kill['got']),'consumer publication differs')
        groups+=1;node_count+=len(computed)
    survivors=json.loads((directory/'survivors.json').read_text());survivor_ids=set()
    for s in survivors:
        cid,hid=s['consumer'],s['helper'];require(0<=cid<len(cs)and 0<=hid<len(profiles),'survivor ID outside grammar')
        identity=cid*len(profiles)+hid;require(identity not in survivor_ids and not seen[identity],'duplicate or killed survivor')
        survivor_ids.add(identity)
        for row,entry in rows.values():
            args=list(map(binary64,row['leaves']));args[4]=tagged_value(entry['values'][hid])
            result=graph_trace(cs[cid]['graph'],args)
            require(B.round(result[cs[cid]['graph']['output']])==binary64(row['raw']['expected']),'survivor contradicts an observed row')
    require(eliminated+len(survivor_ids)==total,'incomplete profile accounting')
    require(LOADED_SOURCE=={**source_hashes(),'timing.py':digest(Path(__file__))},'checker sources changed')
    elementary_kind='raw F2XM1 calls'if direct else'raw exp/log calls'
    receipt=dict(complete=True,direct=direct,consumer_models=consumer_models,programs=total,eliminated=eliminated,survivors=len(survivor_ids),
        checked_consumer_groups=groups,checked_consumer_nodes=node_count,checked_helper_program_rows=component_count,
        elementary_input_constraints=len(elementary),original_rows=len(rows),seconds=time.monotonic()-start,
        checker_sha256=LOADED_SOURCE,inputs={str(p):digest(p)for p in [scoped,*sorted(directory.glob('*.json*'))]},
        claim='All declared helper/consumer IDs are accounted for exactly once. Every grouped contradiction has exact helper arithmetic and full-width consumer nodes; only identical tagged helper values share a trace.',
        limits=f'Conditional on the declared log providers and current cleared {elementary_kind}, whose raw operands are retained and must agree on repeated inputs. Hardware FYL2XP1 is a shared measured native leaf on both backend runs; no software emulation claim. The resident DLMF series and stores are independently checked. No proof these are Excel helpers, no selected or fitted intermediate, no complete PMT implementation. FV0 positive-rate/period normal nonzero balance/target scope; multiple-log mode requires rate<=0.25.')
    with output.open('x')as f:json.dump(receipt,f,indent=2);f.write('\n')
    print(json.dumps(receipt))


def run(scoped, directory, output, factor=False, due=False):
    start=time.monotonic()
    cohorts=('timing','nonunit_timing')if factor else('all','nonunit','timing','nonunit_timing')
    profile_count=1701 if due else 6 if factor else 3**6
    family='due_rate_factor'if due else'annuity_factor'if factor else'timing_complement'
    needed=set()
    for cohort in cohorts:
        for cert in records(directory/f'{cohort}-first-kills.jsonl'):
            needed.update(w['trace_index']for w in cert['proof']['witnesses'])
    rows={}
    for row in records(scoped):
        if row['trace_index']in needed:
            validate_row(row)
            require(row['leaves'][2]in(0,1<<63),'FV outside zero-future projection')
            require(row['trace_index']not in rows,'duplicate original trace identity')
            rows[row['trace_index']]=row
    require(set(rows)==needed,'missing original observation')
    result={}
    for cohort in cohorts:
        seen=set();reasons=Counter();witnesses=0
        for cert in records(directory/f'{cohort}-first-kills.jsonl'):
            c=cert['combine'];cid=c['id'];code=cid;models=[]
            require(0<=cid<profile_count and cid not in seen,'duplicate/out-of-range profile')
            require(cert['candidate_id']==cid,'candidate identity differs')
            operation_kind=('mul','div')[cid//3]if factor else None
            due_form=None;sites=1 if factor else 6
            if due:
                due_form=0 if cid<729 else 1 if cid<1458 else 2
                code=cid-(0,729,1458)[due_form];sites=(6,6,5)[due_form]
                require(c['expression']==due_form,'changed coefficient identity')
            for _ in range(sites):models.append(MODELS[code%3]);code//=3
            require(c['family']==family and models==c['models'],'wrong helper arithmetic')
            if factor:require(c['operation']==operation_kind,'changed factor operation')
            current=None;first=True;keys=set()
            for witness in cert['proof']['witnesses']:
                row=rows[witness['trace_index']]
                require(row['raw']==witness['witness'],'original answer or input changed')
                require(f'0x{row["leaves"][0]:016x}'==cert['rate'] and row['raw']['args'][1]==cert['nper'] and row['raw']['reading']==cert['environment'],'certificate group differs')
                if 'nonunit'in cohort:require(binary64(row['raw']['args'][1])!=1,'unit-period rescue not excluded')
                key=(row['leaves'][0],row['raw']['args'][1],json.dumps(row['raw']['reading'],sort_keys=True))
                if 'timing'in cohort:
                    require(row['leaves'][3]==cert['canonical_timing'],'changed canonical timing')
                    key+=(row['leaves'][3],)
                keys.add(key)
                interval=preimage(row,models,operation_kind,due_form)
                current=interval if first else current.meet(interval)if current else None
                first=False;witnesses+=1
            require(not first and len(keys)==1,'empty or unshared witness group')
            reason='empty_real_intersection'if current is None else'empty_ext80_lattice'
            require(current is None or E.coordinates(current)is None,'a common helper remains')
            require(reason==cert['proof']['reason'],'incorrect emptiness reason')
            if current is not None:require(current.json()==cert['proof']['remaining_interval'],'incorrect endpoint certificate')
            seen.add(cid);reasons[reason]+=1
        require(seen==set(range(profile_count)),'incomplete arithmetic-family exclusion')
        result[cohort]=dict(checked_profiles=len(seen),reasons=reasons,checked_witness_constraints=witnesses,
            certificates_sha256=digest(directory/f'{cohort}-first-kills.jsonl'))
    require(LOADED_SOURCE=={**source_hashes(),'timing.py':digest(Path(__file__))},'checker sources changed')
    receipt=dict(complete=True,cohorts=result,original_killing_rows=len(rows),input_sha256=digest(scoped),
        checker_sha256=LOADED_SOURCE,seconds=time.monotonic()-start,
        family=family,
        claim=f'Each of the {profile_count} declared profiles has a verified original observation group with no shared finite Ext80 helper in every declared cohort. Real-empty intersections also exclude arbitrary real helpers under the specified arithmetic stores.',
        limits='FV0 and the specified single-use-helper expressions only, conditional on original recorded inputs; historical read-back inputs are unavailable. No chosen or fitted helper, no new calculation pass, no conclusion about repeated-q corrections or other PMT structures.')
    with output.open('x')as f:json.dump(receipt,f,indent=2);f.write('\n')
    print(json.dumps(receipt))


if __name__=='__main__':
    if sys.argv[1]=='product-shells':product_shells(*(Path(p)for p in sys.argv[2:5]))
    elif sys.argv[1]=='discount-shells':product_shells(*(Path(p)for p in sys.argv[2:5]),discount=True)
    elif sys.argv[1]=='discount-helpers':discount_helpers(*(Path(p)for p in sys.argv[2:5]))
    elif sys.argv[1]=='pc53-equivalence':pc53_equivalence(*(Path(p)for p in sys.argv[2:6]))
    elif sys.argv[1]=='timing-helpers':timing_helpers(*(Path(p)for p in sys.argv[2:5]),multiple_logs='--logs'in sys.argv[5:],direct='--direct'in sys.argv[5:],pc53='--pc53'in sys.argv[5:])
    else:run(*(Path(p)for p in sys.argv[1:4]),factor='--factor'in sys.argv[4:],due='--due-rate-factor'in sys.argv[4:])
