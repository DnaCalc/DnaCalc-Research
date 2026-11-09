"""Round 41 shape experiments; original MIT code, Tier 1 Excel observations.

Uses the one audited question-service cache. Seeds and candidate shapes are
recorded before requests. No answer table is an application implementation.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import random
import time

import numpy as np
from .fingerprint import AuditedClientMixin, CachedService, bits, value, load_service, sha
from .question_designer import Census, DecisionPath, vector_path
from . import ets

ROOT = Path(__file__).resolve().parents[1]
OUT = Path('/home/ubuntu/codex-runs/round41')
SOURCE = Path(__file__).read_bytes()


def save(name, obj):
    path = OUT / (name + '.json')
    with path.open('x') as stream:
        json.dump(obj, stream, separators=(',', ':'), allow_nan=False)
    return path


def freeze(name, obj):
    return save(name+'-frozen', dict(obj, source_sha256=hashlib.sha256(SOURCE).hexdigest(),
                                    source=SOURCE.decode()))


def service(offline=False,priority=None):
    module, _ = load_service(ets.CHECKOUT)
    class Timed(AuditedClientMixin, module.Client):
        def ask(self, function, rows, **kwargs):
            start = time.monotonic()
            record = dict(function=function, calls=len(rows), result_index=kwargs.get('result_index'))
            try:
                reply = super().ask(function, rows, **kwargs)
                record.update(id=reply['id'], excel_build=reply.get('excel_build'), env_id=reply.get('env_id'))
                return reply
            except Exception as error:
                record.update(error=type(error).__name__, detail=str(error))
                raise
            finally:
                record['seconds'] = time.monotonic()-start
                with (OUT/'service-times.jsonl').open('a') as stream:
                    stream.write(json.dumps(record)+'\n')
                if record['seconds'] > 3:
                    print('Excel wait', function, len(rows), round(record['seconds'],3), flush=True)
    return module, CachedService(None if offline else Timed(timeout=60), OUT/'service-cache', offline,
                                 compresslevel=1, allow_unread_inputs=True,priority=priority)


def histories(count, seed=410401):
    """Eight fixed families, exponents and one-bit regions; no output-guided seed."""
    rng = random.Random(seed)
    for i in range(count):
        family=i%8; n=rng.randrange(8,33); ex=rng.randrange(-240,241)
        if family==0: y=[rng.randrange(-4096,4097)/32 for _ in range(n)]
        elif family==1:
            drift=rng.randrange(-256,257)/16; curve=rng.randrange(-16,17)/128
            y=[drift*j+curve*j*j+rng.randrange(-8,9)/64 for j in range(n)]
        elif family==2:
            at=rng.randrange(n);size=rng.randrange(-256,257)/8
            y=[128.*j+(size if j>=at else 0.) for j in range(n)]
        elif family==3:
            period=rng.randrange(2,9); pattern=[rng.randrange(-128,129)/16 for _ in range(period)]
            drift=rng.randrange(-16,17)/16
            y=[pattern[j%period]+drift*j+rng.randrange(-4,5)/64 for j in range(n)]
        elif family==4:
            center=int(bits(rng.choice([1.,math.pi,1/3])),16)
            y=[value(f'{center+rng.randrange(-32,33):016x}') for _ in range(n)]
        elif family==5:
            at=rng.randrange(n);size=rng.randrange(-256,257)/8
            y=[rng.randrange(-8,9)/64+(size if j==at else 0.) for j in range(n)]
        elif family==6:
            ratio=rng.choice([.5,.75,1.125,1.5,2.]);sign=rng.choice([-1.,1.])
            y=[sign*ratio**j+rng.randrange(-4,5)/64 for j in range(n)]
        else:
            current=float(rng.randrange(-64,65));y=[]
            for j in range(n):
                current+=rng.randrange(-128,129)/16;y.append(current)
        yield dict(y=[math.ldexp(v,ex) for v in y],season=0,family=family,index=i)


def ets_rows(c, cases, variants=(1,2,3)):
    rows=[]
    for z in cases:
        y=c.number_array([[v] for v in z['y']]);t=c.number_array([[j+1.] for j in range(len(z['y']))])
        rows.extend([[y,t,c.number(k),c.number(z.get('season',0))] for k in variants])
    return rows


def census(count=990000, batch=2000, resume=False, cache_only=False):
    design=dict(new_histories=count,prior_histories=10000,seed=410401,batch_histories=batch,
        families=['noise','curvature','step','periodic','adjacent-lattice','impulse','exponential','random-walk'],
        prior='/home/ubuntu/codex-runs/round40/step4-census-answers.json',
        prior_sha256=sha('/home/ubuntu/codex-runs/round40/step4-census-answers.json'))
    if cache_only:
        freeze('ets-census-cached-prefix',dict(original_sha256=sha(OUT/'ets-census-frozen.json'),
            reason='Reconcile a contiguous cached prefix after fresh requests timed out. No new Excel questions. The requested million remains the acquisition target.'))
    elif resume:
        original=ets.read(OUT/'ets-census-frozen.json')
        if any(original[k]!=v for k,v in design.items()):raise ValueError('changed census questions')
        for attempt in range(1,1001):
            continuation='ets-census-continuation'+(''if attempt==1 else'-'+str(attempt))
            if not (OUT/(continuation+'-frozen.json')).exists():break
        else:raise ValueError('too many census continuations')
        freeze(continuation,dict(original_sha256=sha(OUT/'ets-census-frozen.json'),transport_max_calls=24,clients=16,
            reason='Continue identical frozen questions using successful immutable replies. Multi-frame attempts timed out; all positions of one failed request answered in single-frame pieces. Failed acquisitions remain outside numeric counts.'))
    else:freeze('ets-census',design)
    c,s=service(offline=cache_only);summary=Census();by_family=defaultdict(Counter);first={};arrays=[];blocks=[];pending=[]
    if resume and not cache_only:
        from concurrent.futures import ThreadPoolExecutor
        import threading
        workers=ThreadPoolExecutor(max_workers=16);local=threading.local()
        def small_request(rows):
            if not hasattr(local,'service'):local.service=service()[1]
            oracle=local.service;reply=oracle.ask('FORECAST.ETS.STAT',rows)
            return reply['rows'],oracle.receipts[-1]
    prior=ets.read('/home/ubuntu/codex-runs/round40/step4-census-answers.json')
    for i,st in enumerate(ets.grouped(prior)):
        answer=[st[str(k)].get('bits',st[str(k)]) for k in (1,2,3)]
        summary.observe(['round40',i],answer)
        first.setdefault(json.dumps(answer),prior['cases'][i])
    consumed=0
    import inspect
    generator_hash=hashlib.sha256(inspect.getsource(histories).encode()).hexdigest()
    checkpoint=OUT/'census-checkpoint.json'
    if resume and checkpoint.exists():
        state=ets.read(checkpoint);path=OUT/'census-checkpoint.npz'
        if state['generator_sha256']!=generator_hash or state['design']!=design:raise ValueError('changed checkpoint generator/design')
        if sha(path)!=state['array_sha256']:raise ValueError('incomplete derived checkpoint; original replies remain immutable')
        for r in state['receipts']:
            if sha(s.directory/r['file'])!=r['sha256']:raise ValueError('changed checkpoint reply')
        r=state['census'];summary.counts=Counter(r['counts']);summary.first=r['first_questions'];summary.observations=r['observations'];summary.discoveries=r['discovery_positions']
        first=state['first_cases'];by_family=defaultdict(Counter,{int(k):Counter(v)for k,v in state['families'].items()})
        blocks=state['blocks'];s.receipts=state['receipts'];array=np.load(path)['answers'];arrays=[array]
        consumed=summary.observations-10000
        if array.shape!=(consumed,3):raise ValueError('checkpoint shape differs')
        print('verified checkpoint histories',summary.observations,flush=True)
    def checkpoint_state():
        import uuid
        temporary=OUT/'tmp'/('checkpoint-'+uuid.uuid4().hex+'.npz')
        np.savez_compressed(temporary,answers=np.concatenate(arrays))
        state=dict(design=design,generator_sha256=generator_hash,array_sha256=sha(temporary),
            census=summary.report(),first_cases=first,families={k:dict(v)for k,v in by_family.items()},blocks=blocks,receipts=s.receipts,
            limits='Mutable derived recovery cache; never original observations or a publication artifact.')
        temporary.replace(OUT/'census-checkpoint.npz')
        temporary=OUT/'tmp'/('checkpoint-'+uuid.uuid4().hex+'.json')
        with temporary.open('x')as stream:json.dump(state,stream,separators=(',',':'))
        temporary.replace(checkpoint)
    def consume(cases):
        start=time.monotonic();rows=ets_rows(c,cases)
        if cache_only and not s.has_reply('FORECAST.ETS.STAT',rows):
            ans=[];at=0
            while at<len(rows):
                piece=next((rows[at:at+size]for size in (99,24)
                    if s.has_reply('FORECAST.ETS.STAT',rows[at:at+size])),None)
                if piece is None:break
                ans.extend(s.ask('FORECAST.ETS.STAT',piece)['rows']);at+=len(piece)
        elif resume and not s.has_reply('FORECAST.ETS.STAT',rows):
            # The original large request and a smaller retry timed out while
            # a same-input 99-call request succeeded. Keep numeric questions
            # identical, using independent audited clients in single frames.
            ans=[]
            for answered,receipt in workers.map(small_request,[rows[i:i+24]for i in range(0,len(rows),24)]):
                ans.extend(answered);s.receipts.append(receipt)
        else:ans=s.ask_rows('FORECAST.ETS.STAT',rows,chunk_size=15000)
        raw=[]
        for z,j in zip(cases,range(0,len(ans),3)):
            a=[v.get('bits',v) for v in ans[j:j+3]];summary.observe(['round41',z['index']],a)
            key=json.dumps(a);by_family[z['family']][key]+=1;first.setdefault(key,z)
            raw.append([int(v,16) if isinstance(v,str) else 0 for v in a])
        arrays.append(np.array(raw,dtype=np.uint64))
        blocks.append(dict(histories=summary.observations,distinct=len(summary.counts),seconds=time.monotonic()-start))
        if len(blocks)%10==0:
            if not cache_only:checkpoint_state()
            print('ETS histories',summary.observations,'triples',len(summary.counts),'batch seconds',round(blocks[-1]['seconds'],2),flush=True)
        return len(ans)==len(rows)
    stopped=False
    for z in histories(count):
        if z['index']<consumed:continue
        pending.append(z)
        if len(pending)==batch:
            complete=consume(pending);pending=[]
            if not complete:stopped=True;break
    if pending:stopped=not consume(pending)
    if resume and not cache_only:workers.shutdown()
    prefix='ets-census-cached-prefix'if cache_only else'ets-census'
    np.savez_compressed(OUT/(prefix+'-outputs.npz'),answers=np.concatenate(arrays))
    save(prefix+'-result',dict(census=summary.report(),first_cases=first,families={k:dict(v)for k,v in by_family.items()},blocks=blocks,receipts=s.receipts,
        requested_new_histories=count,acquisition_interrupted=stopped,
        limits='Reported STAT triples include postfit beta reset. Finite discovery does not establish a globally finite search or identify any optimum.'))


def ets_boundaries(prior=False):
    prefix='ets-prior-boundaries'if prior else'ets-boundaries'
    if prior:
        original=ets.read('/home/ubuntu/codex-runs/round40/step4-census-answers.json');first={}
        for case,stat in zip(original['cases'],ets.grouped(original)):
            key=json.dumps([stat[str(k)].get('bits',stat[str(k)])for k in [1,2,3]])
            first.setdefault(key,case)
        doc=dict(first_cases=first)
    else:doc=ets.read(OUT/('ets-census-result.json'if (OUT/'ets-census-result.json').exists()else'ets-census-cached-prefix-result.json'))
    by_length=defaultdict(list)
    for answer,z in doc['first_cases'].items():by_length[len(z['y'])].append((answer,z))
    paths=[]
    for n,group in by_length.items():
        for i,(a,z) in enumerate(group):
            b,w=group[(i+1)%len(group)] if len(group)>1 else ('unmeasured',dict(y=[float((j*7)%11-5)for j in range(n)]))
            def norm(y):
                scale=math.frexp(max(map(abs,y)))[1] if any(y) else 0
                return [math.ldexp(v,-scale)for v in y]
            if prior:
                left,right=norm(z['y']),norm(w['y'])
            else:
                # Keep the measured target exactly. A scale change is itself
                # an unproved transformation of the near-constant initializer.
                left=z['y'];scale=math.frexp(max(map(abs,left)))[1] if any(left) else 0
                right=[math.ldexp(v,scale)for v in norm(w['y'])]
            alternatives=[]
            if not prior:
                for y in [[float(j*j)for j in range(n)],
                          [float(j==n//2)for j in range(n)],
                          [float(j)for j in range(n)]]:
                    alternatives.append([math.ldexp(v,scale)for v in norm(y)])
            paths.append(dict(left=left,right=right,alternatives=alternatives,original_left=a,original_right=b))
    freeze(prefix,dict(paths=paths,last=1<<52,coarse_points=65,refinement_points=64))
    c,s=service();results=[]
    for pi,path in enumerate(paths):
        indices=[(1<<52)*i//64 for i in range(65)]
        decision=None;rounds=[];focused=None;anchor=0;unmapped=[];right=path['right']
        for iteration in range(12):
            if not indices:break
            cases=[dict(y=vector_path(path['left'],right,i),season=0)for i in indices]
            ans=s.ask_rows('FORECAST.ETS.STAT',ets_rows(c,cases),chunk_size=24 if not prior else 15000)
            obs=[(i,[a.get('bits',a)for a in ans[j:j+3]])for i,j in zip(indices,range(0,len(ans),3))]
            if decision is None:decision=DecisionPath(obs)
            else:decision.observe(obs)
            rounds.append(dict(points=len(indices),changes=len(decision.changes())))
            target=decision.observations[0]
            if focused is None:
                bracket=next(((lo,hi)for lo,hi in decision.changes()
                    if decision.observations[lo]==target),None)
                if bracket is None:
                    if anchor<len(path['alternatives']):
                        unmapped.append(dict(anchor=anchor,**decision.report()))
                        right=path['alternatives'][anchor];anchor+=1;decision=None;focused=None
                        indices=[(1<<52)*i//64 for i in range(65)]
                        continue
                    indices=[];break
            else:
                lo,hi=focused
                local=sorted(i for i in decision.observations if lo<=i<=hi)
                at=next(j for j,i in enumerate(local)if decision.observations[i]!=target)
                bracket=(local[at-1],local[at])
            focused=bracket;lo,hi=bracket
            count=min(64,hi-lo-1)
            indices=sorted({lo+(hi-lo)*j//(count+1)for j in range(1,count+1)})
        results.append(dict(path=pi,rounds=rounds,focused_edge=focused,anchor=anchor,unmapped_anchor_trials=unmapped,
            focused_adjacent=focused is not None and focused[1]-focused[0]==1,
            original_endpoint_preserved=json.loads(path['original_left'])==decision.observations[0],
            **decision.report()))
        print('ETS path',pi,'observations',len(decision.observations),'changes',len(decision.changes()),'adjacent',results[-1]['adjacent'],flush=True)
    save(prefix+'-result',dict(paths=results,receipts=s.receipts,
        limits='Only observed changes on these paths, including unresolved jagged regions. Adjacent path indices may repeat binary64 vectors. No global decision region is inferred.'))


def normal(resume=False):
    from .normal_quantile import CENTRAL_N,CENTRAL_D,tail
    rng=random.Random(410201);centers=[]
    lower=[2.**-k for k in [1,2,3,4,8,16,32,64,128,256,512,1022,1070]]
    lower += [.075,.076,.1,.2,.3,.4,.499999999999]
    lower += [math.exp(rng.uniform(-740,math.log(.5)))for _ in range(64-len(lower))]
    upper=[.5,.500000000001,.55,.6,.7,.8,.9,.924,.925,.93,.9375]
    upper += [rng.uniform(.5,.9375)for _ in range(64-len(upper))]
    for region,ps in [('lower',lower),('upper',upper)]:
        for p in ps:
            start=max(1,int(bits(p),16)-4096)
            if region=='lower':start=min(start,int(bits(.5),16)-8191)
            else:start=max(int(bits(.5),16),min(start,int(bits(.9375),16)-8191))
            centers.append(dict(region=region,start=f'{start:016x}',count=8192))
    if resume:
        if ets.read(OUT/'normal-dense-frozen.json')['windows'] != centers:raise ValueError('changed normal questions')
        freeze('normal-dense-continuation',dict(original=sha(OUT/'normal-dense-frozen.json'),
            reason='Retain explicit input-changed unread rows as acquisition failures; never numeric matches. Same frozen windows and candidates.'))
    else:freeze('normal-dense',dict(windows=centers,candidate='Published AS241 binary64 Horner, q*N/D and q*(N/D); no fit',seed=410201))
    c,s=service();results=[];raw=[]
    for i,w in enumerate(centers):
        start=int(w['start'],16);points=[value(f'{j:016x}')for j in range(start,start+w['count'])]
        def predict(p,after=False):
            if abs(p-.5)>.425:return bits(tail(p))
            q=p-.5;r=.180625-q*q;n=CENTRAL_N[0];d=CENTRAL_D[0]
            for a,b in zip(CENTRAL_N[1:],CENTRAL_D[1:]):n=n*r+a;d=d*r+b
            return bits(q*(n/d)if after else(q*n)/d)
        preds=[[predict(p,False)for p in points],[predict(p,True)for p in points]]
        answers=s.ask_rows('NORM.S.INV',[[c.number(p)]for p in points],chunk_size=16000)
        actual=[a.get('bits')for a in answers];raw.append([int(a,16)if a else 0 for a in actual])
        def ordered(b):
            x=int(b,16);return (~x & ((1<<64)-1))if x>>63 else x | (1<<63)
        jumps=Counter(ordered(b)-ordered(a)for a,b in zip(actual,actual[1:])if a and b)
        contradictions=[]
        for pred in preds:
            contradictions.append(next((dict(index=j,input=bits(points[j]),prediction=p,answer=actual[j])for j,p in enumerate(pred)if actual[j] is not None and p!=actual[j]),None))
        delta_jumps=[]
        for pred in preds:
            delta_jumps.append(sum((ordered(actual[j+1])-ordered(actual[j]))!=(ordered(pred[j+1])-ordered(pred[j]))for j in range(len(points)-1)if actual[j]and actual[j+1]))
        results.append(dict(**w,jumps=dict(jumps),candidate_first_contradictions=contradictions,
            different_step_positions=delta_jumps,non_numeric=sum(a is None for a in actual)))
        if i%16==15:print('normal windows',i+1,'inputs',(i+1)*8192,flush=True)
    np.savez_compressed(OUT/'normal-dense-outputs.npz',answers=np.array(raw,dtype=np.uint64))
    save('normal-dense-result',dict(windows=results,selections=len(centers)*8192,receipts=s.receipts,
        limits='Mapping and step positions precede any fit. Candidate failures refute these frozen forms; no extra polynomial pieces or coefficient fit added.'))


def regression():
    import itertools
    from .regression import service_design
    cases=[];designs={}
    for mode in ['rank','decomposition','residual']:
        designed,selections=service_design(mode)
        start=len(cases)
        cases.extend(dict(z,design=mode,global_index=start+i)for i,z in enumerate(designed))
        designs[mode]=dict(start=start,end=len(cases),selections=selections)
    start=len(cases)
    for n,constant,zero,reverse,swap,ex,response in itertools.product(
            [4,8,16,32],[False,True],['none','before','after'],[False,True],
            [False,True],[None,-22,-18,-12],range(3)):
        delta=0. if ex is None else 2.**ex
        x=[[float(j+1),2.*(j+1)+(delta if j==n-1 else 0.)]for j in range(n)]
        if swap:x=[r[::-1]for r in x]
        if zero=='before':x=[[0.]+r for r in x]
        if zero=='after':x=[r+[0.]for r in x]
        y=[[float(j==response)]for j in range(n)]
        if reverse:x,y=x[::-1],y[::-1]
        cases.append(dict(x=x,y=y,constant=constant,zero=zero,reverse=reverse,swap=swap,
                          delta_exponent=ex,response=response,n=n,design='zero-rank',global_index=len(cases)))
    designs['zero-rank']=dict(start=start,end=len(cases),selections=[(1,1),(1,2),(1,3),(4,2),(5,2)])
    freeze('regression',dict(cases=cases,designs=designs,
        method='Common rank/reflector/residual designs plus zero-column position, reversal and near-collinear basis responses. LOGEST uses exact powers of two for all response values.'))
    c,s=service();results={};counts=Counter()
    for function in ['LINEST','LOGEST']:
        for mode,design in designs.items():
            selected=cases[design['start']:design['end']]
            args=[[c.number_array([[2.**v[0]]for v in z['y']]if function=='LOGEST'else z['y']),
                c.number_array(z['x']),c.boolean(z['constant']),c.boolean(True)]for z in selected]
            results[function+'/'+mode]={}
            for selection in design['selections']:
                ans=s.ask_rows(function,args,chunk_size=4000,result_index=selection)
                results[function+'/'+mode][str(selection)]=ans
                counts.update(function+'/'+a['type']for a in ans)
            print('regression',function,mode,len(selected),'selections',len(design['selections']),flush=True)
    save('regression-result',dict(results=results,counts=counts,receipts=s.receipts,
        limits='Rank and residual mapping, with every numeric/error/unread row retained. LOGEST log transforms and full QR are not inferred from agreement counts.'))


def replay():
    from .ets_proof import application_replay
    application_replay(ROOT/'target/debug/symbolic-racer',OUT/'registered-replay',
        directories=[OUT/'service-cache'],functions=['LINEST','LOGEST','BINOM.DIST','BINOM.INV','CRITBINOM','YIELD','ODDFYIELD'])


def residual_isolation():
    """Coordinate predictor leaves every residual operand exactly unchanged."""
    from .x87_programs import Builder, Native
    rng=np.random.default_rng(410302);raw=rng.integers(0x3fe0000000000000,0x4000000000000000,size=2000000,dtype=np.uint64)
    xs=raw.view(np.float64);strict=xs*xs
    wide=(xs.astype(np.longdouble)*xs.astype(np.longdouble)).astype(np.float64)
    selected=np.flatnonzero(strict.view(np.uint64)!=wide.view(np.uint64))
    points=[float(xs[i])for i in selected[:384]]+[float(x)for x in xs[:128]]
    native=Native('/home/ubuntu/codex-runs/round40/x87-round40.so');g=Builder(1)
    square=g.op('mul',0,0,store=True);graph=g.graph(square,square)
    pred=[f'{int(v):016x}'for v in native.grid(graph,np.array([[int(bits(x),16)]for x in points],dtype=np.uint64))[:,0]]
    cases=[dict(x=[[1.],[0.],[0.]],y=[[0.],[x],[0.]],constant=False)for x in points]
    freeze('residual-isolation',dict(cases=cases,strict=[bits(x*x)for x in points],stored_pc64=pred,
        native_sha256=sha(native.path),graph=graph,seed=410302,
        method='Exact coordinate column: reflector changes only the first response coordinate (zero). Tail residual is the original x,0, so selected SSresid reads its square without a QR approximation.'))
    c,s=service();ans=s.ask_rows('LINEST',[[c.number_array(z['y']),c.number_array(z['x']),c.boolean(False),c.boolean(True)]for z in cases],result_index=(5,2))
    save('residual-isolation-result',dict(answers=ans,receipts=s.receipts,
        stored_pc64_contradictions=[i for i,(p,a)in enumerate(zip(pred,ans))if p!=a.get('bits')],
        strict_contradictions=[i for i,(x,a)in enumerate(zip(points,ans))if bits(x*x)!=a.get('bits')]))


def residual_sums():
    from .x87_programs import Builder,Native
    from .question_designer import choose_split
    native=Native('/home/ubuntu/codex-runs/round40/x87-round40.so')
    rng=np.random.default_rng(410303)
    pool=np.ldexp(rng.uniform(.5,1.,size=(100000,4)),rng.integers(-20,21,size=(100000,4)))
    raw=pool.view(np.uint64);predictions={};graphs={}
    for model,cw,stored in [('strict',0x123f,True),('stored_pc64',0x133f,True),('resident_pc64',0x133f,False)]:
        for product_store in [False,True]:
            if model=='strict' and not product_store:continue
            for order in ['forward','reverse','pair']:
                g=Builder(4);squares=[g.op('mul',i,i,store=product_store)for i in range(4)]
                zero=g.constant(bits(0.))
                def total(ids):
                    if order=='reverse':ids=ids[::-1]
                    if order=='pair'and len(ids)>1:
                        middle=len(ids)//2
                        return g.op('add',total(ids[:middle]),total(ids[middle:]),store=stored)
                    result=zero
                    for k in ids:result=g.op('add',result,k,store=stored)
                    return result
                resid=total(squares[1:]);full=total(squares)
                reg=g.op('sub',full,resid,store=True);graph=g.graph(resid,reg)
                label=f'{model}/products{product_store}/{order}'
                result=native.grid(graph,raw,cw,(1<<g.sites)-1)
                predictions[label]=result;graphs[label]=dict(graph=graph,cw=cw,mask=(1<<g.sites)-1)
    # Select a bounded union of pre-answer disagreements, plus independent controls.
    ids=set(range(64));baseline=predictions['stored_pc64/productsTrue/forward']
    for name,pred in predictions.items():
        ids.update(map(int,np.flatnonzero(np.any(pred!=baseline,axis=1))[:64]))
    ids=sorted(ids);points=pool[ids].tolist()
    frozen_pred={name:[[f'{int(v):016x}'for v in row]for row in pred[ids]]for name,pred in predictions.items()}
    freeze('residual-sums',dict(points=points,predictions=frozen_pred,graphs=graphs,seed=410303,
        native_sha256=sha(native.path),method='Coordinate column e0 isolates original response tail. Fixed forward, reverse and pair sums, fixed precision/store alternatives; pre-answer distinguishing union plus controls.'))
    c,s=service();args=[[c.number_array([[v]for v in y]),c.number_array([[1.],[0.],[0.],[0.]]),c.boolean(False),c.boolean(True)]for y in points]
    answers=[s.ask_rows('LINEST',args,result_index=selection)for selection in [(5,2),(5,1)]]
    actual=[[answers[j][i].get('bits')for j in range(2)]for i in range(len(points))]
    save('residual-sums-result',dict(answers=answers,receipts=s.receipts,
        programs={name:dict(first_contradiction=next((i for i,(p,a)in enumerate(zip(pred,actual))if p!=a),None),
            failures=sum(p!=a for p,a in zip(pred,actual)))for name,pred in frozen_pred.items()}))


def residual_fresh():
    from fractions import Fraction as Q
    from checker.exact import operation
    from .x87_programs import Builder,Native
    rng=random.Random(410304);cases=[];pred=[]
    def total(y):
        result=Q(0)
        for v in y:
            square=operation(Q(v),Q(v),'mul','stored_pc64')
            result=operation(result,square,'add','stored_pc64')
        return result
    for i in range(1024):
        n=rng.randrange(3,33);p=rng.choice([1,2]);y=[math.ldexp(rng.uniform(-1.,1.),rng.randrange(-200,201))for _ in range(n)]
        x=[[float(j==k)for k in range(p)]for j in range(n)]
        residual=total(y[p:]);full=total(y);reg=operation(full,-residual,'add','stored_pc64')
        pred.append([bits(float(residual)),bits(float(reg))]);cases.append(dict(y=y,x=x,p=p))
    # A second implementation checks every chosen native model on the original
    # four-coordinate distinguishing design before any new answer arrives.
    original=ets.read(OUT/'residual-sums-frozen.json')
    for y,wanted in zip(original['points'],original['predictions']['stored_pc64/productsTrue/forward']):
        residual=total(y[1:]);reg=operation(total(y),-residual,'add','stored_pc64')
        if [bits(float(residual)),bits(float(reg))]!=wanted:raise ValueError('independent rational prediction disagrees')
    freeze('residual-fresh',dict(cases=cases,predictions=pred,seed=410304,
        independent_replay_cases=len(original['points']),independent_source_sha256=sha(ROOT/'checker/exact.py'),
        method='Fresh lengths 3..32, one/two exact coordinate predictors, independent Fraction RN53(RN64) at every square and forward sum. All response exponents -200..200.'))
    c,s=service();args=[[c.number_array([[v]for v in z['y']]),c.number_array(z['x']),c.boolean(False),c.boolean(True)]for z in cases]
    ans=[s.ask_rows('LINEST',args,result_index=selection)for selection in [(5,2),(5,1)]]
    actual=[[ans[j][i].get('bits')for j in range(2)]for i in range(len(cases))]
    save('residual-fresh-result',dict(answers=ans,receipts=s.receipts,
        failures=[dict(case=i,prediction=p,answer=a)for i,(p,a)in enumerate(zip(pred,actual))if p!=a],
        selections=2*len(cases)))


def residual_additions():
    from fractions import Fraction as Q
    from checker.exact import operation
    from .x87_programs import Builder,Native
    native=Native('/home/ubuntu/codex-runs/round40/x87-round40.so')
    g=Builder(4);squares=[g.op('mul',i,i,store=True)for i in range(4)]
    zero=g.constant(bits(0.));g.phase=1
    def total(ids):
        result=zero
        for k in ids:result=g.op('add',result,k,store=True)
        return result
    residual=total(squares[1:]);reg=g.op('sub',total(squares),residual,store=True)
    graph=g.graph(residual,reg);mask=(1<<g.sites)-1
    rng=np.random.default_rng(410305)
    pool=np.ldexp(rng.uniform(.5,1.,size=(500000,4)),rng.integers(-20,21,size=(500000,4)))
    good=native.grid(graph,pool.view(np.uint64),0x133f,mask,inner_cw=0x133f)
    strict_add=native.grid(graph,pool.view(np.uint64),0x133f,mask,inner_cw=0x123f)
    ids=sorted(set(map(int,np.flatnonzero(good[:,0]!=strict_add[:,0])[:128]))|
               set(map(int,np.flatnonzero(good[:,1]!=strict_add[:,1])[:128]))|set(range(64)))
    points=pool[ids].tolist();pred=[[f'{int(v):016x}'for v in row]for row in good[ids]]
    other=[[f'{int(v):016x}'for v in row]for row in strict_add[ids]]
    # Independent exact arithmetic validates the native per-phase precision.
    for y,wanted,alternative in zip(points,pred,other):
        products=[operation(Q(v),Q(v),'mul','stored_pc64')for v in y]
        results=[]
        for mode in ['stored_pc64','strict']:
            def reduce(items):
                result=Q(0)
                for v in items:result=operation(result,v,'add',mode)
                return result
            residual=reduce(products[1:]);reg=operation(reduce(products),-residual,'add',mode)
            results.append([bits(float(residual)),bits(float(reg))])
        if results!=[wanted,alternative]:raise ValueError('mixed-precision rational replay differs')
    freeze('residual-additions',dict(points=points,predictions=pred,strict_additions=other,
        graph=graph,mask=mask,product_cw=0x133f,addition_cws=[0x133f,0x123f],seed=410305,
        native_sha256=sha(native.path),exact_source_sha256=sha(ROOT/'checker/exact.py'),
        method='Pre-answer mixed-precision separators with identical stored PC64 products, plus independent controls. Independent Fraction rounding checks both programs.'))
    c,s=service();args=[[c.number_array([[v]for v in y]),c.number_array([[1.],[0.],[0.],[0.]]),c.boolean(False),c.boolean(True)]for y in points]
    ans=[s.ask_rows('LINEST',args,result_index=selection)for selection in [(5,2),(5,1)]]
    actual=[[ans[j][i].get('bits')for j in range(2)]for i in range(len(points))]
    save('residual-additions-result',dict(answers=ans,receipts=s.receipts,selections=len(points)*2,
        failures=[i for i,(p,a)in enumerate(zip(pred,actual))if p!=a],
        strict_addition_contradictions=[i for i,(p,a)in enumerate(zip(other,actual))if p!=a]))


def binomial():
    import itertools
    from types import SimpleNamespace
    rng=random.Random(410401)
    groups=[dict(n=n,p=p)for n,p in itertools.product([6,10,20,32,64,128],[.125,.3,.5,.7,.875])]
    groups += [dict(n=rng.randrange(3,129),p=rng.uniform(.02,.98))for _ in range(16)]
    freeze('binomial',dict(groups=groups,seed=410401,
        method='Full PMF/CDF prefixes, lower/upper endpoint powers, and inverse edges near public CDF. Conditional mass/power leaves remain explicit.'))
    c,s=service();labels=[(i,k)for i,z in enumerate(groups)for k in range(z['n']+1)]
    args=[[c.number(k),c.number(groups[i]['n']),c.number(groups[i]['p'])]for i,k in labels]
    pmf=s.ask_rows('BINOM.DIST',[a+[c.boolean(False)]for a in args],chunk_size=12000)
    cdf=s.ask_rows('BINOM.DIST',[a+[c.boolean(True)]for a in args],chunk_size=12000)
    powerargs=[[c.number(base),c.number(z['n'])]for z in groups for base in [1-z['p'],z['p']]]
    powers=s.ask_rows('POWER',powerargs)
    mass=[{}for _ in groups];tails=[{}for _ in groups]
    for (i,k),a,b in zip(labels,pmf,cdf):mass[i][k]=a;tails[i][k]=b
    selected=[]
    for i,z in enumerate(groups):
        n,p=z['n'],z['p'];sd=math.sqrt(n*p*(1-p))
        ks=sorted({0,max(0,int(n*p-sd)),min(n-1,int(n*p)),min(n-1,int(n*p+sd))})
        for k in ks:
            a=tails[i][k]
            if a['type']=='number' and 0<value(a['bits'])<1:selected.append(dict(group=i,k=k,cdf=a['bits']))
    # Conditional forward/backward recurrence programs are fixed before edges.
    predictions={};prefix_checks=[]
    def published_edge(x,inclusive):
        f=float(x)
        if not 0<f<1:return None
        raw=int(bits(f),16)
        if np.longdouble(f)<x or (inclusive and np.longdouble(f)==x):raw+=1
        return f'{raw:016x}'
    for mode in ['strict','stored_pc64','resident_pc64']:
        def op(kind,a,b):
            if mode=='strict':
                a,b=float(a),float(b)
                return a+b if kind=='add' else a*b if kind=='mul' else a/b
            a,b=np.longdouble(a),np.longdouble(b)
            v=a+b if kind=='add' else a*b if kind=='mul' else a/b
            return float(v)if mode=='stored_pc64'else v
        for seed in ['power','public_mass']:
            for shape in range(6):
                for direction in ['lower','upper']:
                    label=f'{mode}/{seed}/{shape}/{direction}';thresholds={}
                    for i,z in enumerate(groups):
                        n,p=z['n'],z['p'];q=1-p
                        term=value(powers[2*i+(direction=='upper')]['bits'])if seed=='power'else value(mass[i][n if direction=='upper'else 0]['bits'])
                        cumulative=0.
                        for k in (range(n+1)if direction=='lower'else range(n,-1,-1)):
                            cumulative=op('add',cumulative,term)
                            if direction=='lower':thresholds[i,k]=cumulative
                            else:thresholds[i,k-1]=op('add',1.,-cumulative)
                            a,b=(n-k,k+1)if direction=='lower'else(k,n-k+1)
                            pp,qq=(p,q)if direction=='lower'else(q,p)
                            if not a:continue
                            if shape==0:term=op('div',op('mul',op('div',op('mul',term,a),b),pp),qq)
                            elif shape==1:term=op('div',op('mul',op('div',op('mul',term,pp),qq),a),b)
                            elif shape==2:term=op('div',op('mul',term,op('mul',a,pp)),op('mul',b,qq))
                            elif shape==3:term=op('mul',term,op('div',op('mul',a,pp),op('mul',b,qq)))
                            elif shape==4:term=op('mul',op('mul',op('div',term,b),a),op('div',pp,qq))
                            else:term=op('mul',op('mul',term,op('div',a,b)),op('div',pp,qq))
                    for inclusive in [False,True]:
                        predictions[label+f'/equal{inclusive}']=[published_edge(thresholds[z['group'],z['k']],inclusive)for z in selected]
        # Public-mass accumulation is an explicitly conditional alternative.
        for direction in ['lower','upper']:
            for reverse in [False,True]:
                values_=[]
                for z in selected:
                    i,k=z['group'],z['k'];n=groups[i]['n']
                    ids=list(range(k+1))if direction=='lower'else list(range(k+1,n+1))
                    if reverse:ids.reverse()
                    total=0.
                    for j in ids:total=op('add',total,value(mass[i][j]['bits']))
                    if direction=='upper':total=op('add',1.,-total)
                    values_.append(total)
                for inclusive in [False,True]:predictions[f'{mode}/public_sum/{direction}/reverse{reverse}/equal{inclusive}']=[published_edge(v,inclusive)for v in values_]
    freeze('binomial-edges',dict(selected=selected,predictions=predictions,
        component_receipts=list(s.receipts),
        limits='A public POWER or PMF is a conditional seed, not a private primitive. Resident comparison edges use the binary64 ceiling of the retained sum; equal comparison alternatives remain separate.'))
    edges=[];proxy=SimpleNamespace(uses_socket=True,ask=s.ask)
    for j,z in enumerate(selected):
        i,k=z['group'],z['k'];g=groups[i];raw=int(z['cdf'],16)
        predicate=lambda _,a,k=k: value(a['bits'])>k if a['type']=='number'else None
        for radius in [4,64,4096,1<<20]:
            lo=max(1,raw-radius);hi=min(int(bits(1.),16)-1,raw+radius)
            boundary_args=[c.number(g['n']),c.number(g['p'])]
            ends=s.ask('BINOM.INV',[boundary_args+[c.number(f'{v:016x}')]for v in [lo,hi]])['rows']
            if predicate(None,ends[0])!=predicate(None,ends[1]):
                edge=c.Client.bisect_edge(proxy,'BINOM.INV',f'{lo:016x}',f'{hi:016x}',predicate,
                                         arguments=boundary_args,index=2,points=64)
                edge.update(group=i,k=k,cdf=z['cdf'],cdf_offset=int(edge['high'],16)-raw);break
        else:edge=dict(group=i,k=k,cdf=z['cdf'],unresolved=True)
        edges.append(edge)
        if j%32==31:print('binomial inverse edges',j+1,'/',len(selected),flush=True)
    aliasargs=[[c.number(groups[e['group']]['n']),c.number(groups[e['group']]['p']),c.number(e[side])]
               for e in edges if 'high'in e for side in ['low','high']]
    aliases=s.ask_rows('CRITBINOM',aliasargs,chunk_size=12000)
    expected=[e[side+'_answer']for e in edges if 'high'in e for side in ['low','high']]
    programs={name:dict(first_contradiction=next((i for i,(v,e)in enumerate(zip(pred,edges))if e.get('high')!=v),None),
        contradictions=sum(e.get('high')!=v for e,v in zip(edges,pred)))for name,pred in predictions.items()}
    save('binomial-result',dict(groups=groups,labels=labels,pmf=pmf,cdf=cdf,powers=powers,edges=edges,
        programs=programs,alias_answers=aliases,alias_disagreements=[i for i,(a,b)in enumerate(zip(aliases,expected))if a!=b],
        offsets=dict(Counter(e['cdf_offset']for e in edges if 'cdf_offset'in e)),receipts=s.receipts,
        limits='Every edge is an observed adjacent-alpha transition for one fixed n,p,k. Monotonic quantile comparison is the bisection premise. Conditional programs and unresolved brackets never become matches.'))


def binomial_fresh():
    rng=random.Random(410402);groups=[]
    for i in range(128):
        n=rng.randrange(2,257)
        if i%4==0:p=math.ldexp(1.,-rng.randrange(1,19))
        elif i%4==1:p=1-math.ldexp(1.,-rng.randrange(1,19))
        elif i%4==2:p=value(f'{int(bits(rng.choice([.125,.3,.5,.7,.875])),16)+rng.randrange(-64,65):016x}')
        else:p=rng.uniform(.01,.99)
        groups.append(dict(n=n,p=p))
    freeze('binomial-fresh-design',dict(groups=groups,seed=410402,
        shape='Initial count-zero mass uses strict alpha < CDF(0). Thereafter select first k>=1 with CDF(k)>=alpha. Public CDF is an independent conditional leaf.'))
    c,s=service();labels=[(i,k)for i,z in enumerate(groups)for k in range(z['n']+1)]
    cdfs=s.ask_rows('BINOM.DIST',[[c.number(k),c.number(groups[i]['n']),c.number(groups[i]['p']),c.boolean(True)]for i,k in labels],chunk_size=15000)
    rows=[[]for _ in groups]
    for (i,k),a in zip(labels,cdfs):rows[i].append(a)
    cases=[];pred=[]
    for i,(z,cs)in enumerate(zip(groups,rows)):
        candidates=[]
        for k,a in enumerate(cs[:-1]):
            if a['type']!='number' or not 0<value(a['bits'])<1:continue
            raw=int(a['bits'],16)
            if k and raw-int(cs[k-1].get('bits','0'),16)<16:continue
            if int(cs[k+1].get('bits','0'),16)-raw<16:continue
            candidates.append(k)
        if len(candidates)>5:candidates=[candidates[j*(len(candidates)-1)//4]for j in range(5)]
        for k in candidates:
            raw=int(cs[k]['bits'],16)
            for step in [-2,-1,0,1,2]:
                alpha=value(f'{raw+step:016x}')
                if alpha<value(cs[0]['bits']):answer=0
                else:answer=next(j for j,a in enumerate(cs[1:],1)if value(a['bits'])>=alpha)
                cases.append(dict(group=i,k=k,alpha=bits(alpha),offset=step));pred.append(bits(answer))
    freeze('binomial-fresh',dict(cases=cases,predictions=pred,component_receipts=s.receipts.copy(),
        shape='One independent CDF, strict initial endpoint then inclusive positive counts; no probability-dependent rule or fitted coefficient.'))
    args=[[c.number(groups[z['group']]['n']),c.number(groups[z['group']]['p']),c.number(z['alpha'])]for z in cases]
    answers={f:s.ask_rows(f,args,chunk_size=15000)for f in ['BINOM.INV','CRITBINOM']}
    save('binomial-fresh-result',dict(groups=groups,labels=labels,cdf=cdfs,answers=answers,receipts=s.receipts,
        failures={f:[dict(case=i,prediction=p,answer=a)for i,(p,a)in enumerate(zip(pred,ans))if p!=a.get('bits')]for f,ans in answers.items()},
        selections=len(cases)*2,limits='Complete inverse prediction conditional on the separately measured public CDF. CDF arithmetic, hidden internal reuse and search order remain unidentified.'))


def yield_plateaus():
    from datetime import date
    from types import SimpleNamespace
    serial=lambda d:float((d-date(1899,12,30)).days)
    cases=[]
    for n in [2,4,8]:
        for frequency in [1,2]:
            for coupon in [0.,.0625]:
                year,month=divmod(2024*12+n*12//frequency,12)
                first_year,first_month=divmod(2024*12+12//frequency,12)
                normal=[serial(date(2024,1,15)),serial(date(year,month+1,15)),coupon,100.,100.,float(frequency),0.]
                odd=[normal[0],normal[1],serial(date(2023,1,15)),serial(date(first_year,first_month+1,15)),coupon,100.,100.,float(frequency),0.]
                cases.append(dict(n=n,frequency=frequency,coupon=coupon,normal=normal,odd=odd))
    offsets=sorted({0.,*[sign*2.**e for sign in [-1,1]for e in range(-40,6)]})
    freeze('yield-plateaus',dict(cases=cases,price_offsets=offsets,
        method='Ordinary bonds settle on coupon-calendar dates; ODDFYIELD has a long first coupon and different positive-coupon cashflows. Counts 2/4/8. Map each function separately by price exponent ladder. A return plateau alone identifies no iterate count.'))
    c,s=service();proxy=SimpleNamespace(uses_socket=True,ask=s.ask);results={}
    for function in ['YIELD','ODDFYIELD']:
        price_index=3 if function=='YIELD'else 5;key='normal'if function=='YIELD'else'odd'
        rows=[]
        for z in cases:
            for offset in offsets:
                a=z[key].copy();a[price_index]=100.+offset;rows.append([c.number(v)for v in a])
        ans=s.ask_rows(function,rows,chunk_size=12000);per=[]
        for i,z in enumerate(cases):
            a=ans[i*len(offsets):(i+1)*len(offsets)];middle=offsets.index(0.);at_par=a[middle];edges=[]
            for direction in [-1,1]:
                ids=list(range(middle-1,-1,-1))if direction<0 else list(range(middle+1,len(offsets)))
                unlike=next((j for j in ids if a[j]!=at_par),None)
                if unlike is None:continue
                low,high=sorted([100.,100.+offsets[unlike]])
                arguments=[c.number(v)for j,v in enumerate(z[key])if j!=price_index]
                predicate=lambda _,answer,at_par=at_par: answer==at_par
                edge=c.Client.bisect_edge(proxy,function,bits(low),bits(high),predicate,
                    arguments=arguments,index=price_index,points=64)
                edge['direction']=direction;edges.append(edge)
            per.append(dict(case=i,at_par=at_par,plateau_edges=edges))
        results[function]=dict(answers=ans,cases=per)
    save('yield-plateaus-result',dict(results=results,receipts=s.receipts,
        limits='Adjacent price edges are measured, without assuming an iteration producer or assigning iterate counts. The previous refuted RATE/secant/Brent/fixed-count programs are not expanded to absorb these observations. No stop tolerance is installed.'))


def binomial_beta():
    cases=[]
    for name in ['binomial-result','binomial-fresh-result']:
        d=ets.read(OUT/(name+'.json'))
        for (i,k),a in zip(d['labels'],d['cdf']):
            z=d['groups'][i]
            if 0<k<z['n']:
                cases.append(dict(source=name,group=i,k=k,n=z['n'],p=z['p'],cdf=a))
    freeze('binomial-beta',dict(cases=cases,
        relation='BINOM.DIST(k,n,p,TRUE) = BETA.DIST(binary64(1-p),n-k,k+1,TRUE); positive interior counts only.'))
    c,s=service();answers=s.ask_rows('BETA.DIST',[[c.number(1-z['p']),c.number(z['n']-z['k']),c.number(z['k']+1),c.boolean(True)]for z in cases],chunk_size=15000)
    save('binomial-beta-result',dict(answers=answers,receipts=s.receipts,selections=len(cases),
        contradictions=[dict(case=i,binomial=z['cdf'],beta=a)for i,(z,a)in enumerate(zip(cases,answers))if z['cdf']!=a],
        limits='Observable public-function coupling, not proof of private code reuse. Conditional inverse results do not repair the current beta producer.'))


def ets_search():
    """Rank two fixed structural families, retaining their first contradictions."""
    from .ets_optimizer import objective,GRID
    from .ets_round40 import fit_with_trend
    components=ets.read(OUT/'ets-search-components.json')
    reported=[tuple(value(v)for v in row[:2])for row in components['triples']]
    # This is an explicitly learned, frozen hypothesis, never an application
    # table. Reported beta=.001 may be a postfit reset, so test both complete
    # reported pairs and nonreset pairs plus the prior frozen zero-beta entries.
    tables={'original23':list(GRID),'reported30':reported,
            'nonreset25':sorted(set([v for v in reported if v[1]>.001]+[v for v in GRID if v[1]==.001]))}
    cases=list(histories(128,410406));predictions={};raw_predictions={}
    def finish(y,pair):
        a,b=pair
        _,_,_,_,reset=fit_with_trend(y,a,b,ets.initial_trend(y))
        return [bits(a),bits(.001 if reset else b),bits(2.**-52)]
    for table,pool in tables.items():
        for order in ['first','last']:
            label=f'table/{table}/{order}';raw=[]
            for z in cases:
                y=z['y'];ex=math.frexp(max(map(abs,y)))[1]if any(y)else 0
                scaled=[math.ldexp(v,-ex)for v in y]
                scores=[objective(scaled,a,b)for a,b in pool];best=min(scores)
                ids=[i for i,score in enumerate(scores)if score==best]
                raw.append(pool[ids[0]if order=='first'else ids[-1]])
            raw_predictions[label]=raw;predictions[label]=[finish(z['y'],p)for z,p in zip(cases,raw)]
    def golden(y,a,b,coordinate,steps):
        low,high=(max(.002,b+.001),.998)if coordinate==0 else(.001,a-.001)
        ratio=(math.sqrt(5.)-1.)/2
        def score(x):return objective(y,x,b)if coordinate==0 else objective(y,a,x)
        left=high-ratio*(high-low);right=low+ratio*(high-low)
        fl,fr=score(left),score(right)
        for _ in range(steps):
            if fl<=fr:
                high,right,fr=right,left,fl;left=high-ratio*(high-low);fl=score(left)
            else:
                low,left,fl=left,right,fr;right=low+ratio*(high-low);fr=score(right)
        return left if fl<=fr else right
    for order in [(0,1),(1,0)]:
        for steps in [8,24]:
            for sweeps in [1,3]:
                label=f'golden/{order}/{steps}/{sweeps}';raw=[]
                for z in cases:
                    y=z['y'];ex=math.frexp(max(map(abs,y)))[1]if any(y)else 0
                    y=[math.ldexp(v,-ex)for v in y];a,b=.5,.25
                    for _ in range(sweeps):
                        for coordinate in order:
                            point=golden(y,a,b,coordinate,steps)
                            if coordinate==0:a=point
                            else:b=point
                    raw.append((float(format(a,'.6g')),float(format(b,'.6g'))))
                raw_predictions[label]=raw;predictions[label]=[finish(z['y'],p)for z,p in zip(cases,raw)]
    freeze('ets-search',dict(cases=cases,tables=tables,predictions=predictions,raw_parameters=raw_predictions,
        component_inputs=components,component_sha256=sha(OUT/'ets-search-components.json'),seed=410406,
        ranking=['Fixed table, stable first/last exact-score tie rule','Coordinate golden section, fixed seed/box/order/iteration counts'],
        limits='Candidate tables are explicit discovery hypotheses, not recovered private table loads. Existing binary64 initial trend, SSE and postfit-reset programs are conditional. No point or iteration count will be appended after fresh answers.'))
    c,s=service();ans=s.ask_rows('FORECAST.ETS.STAT',ets_rows(c,cases),chunk_size=15000)
    actual=[[a.get('bits',a)for a in ans[j:j+3]]for j in range(0,len(ans),3)]
    save('ets-search-result',dict(answers=ans,receipts=s.receipts,
        programs={name:dict(first_contradiction=next((dict(case=i,prediction=p,answer=a)for i,(p,a)in enumerate(zip(pred,actual))if p!=a),None),
            contradictions=sum(p!=a for p,a in zip(pred,actual)))for name,pred in predictions.items()},
        stop='Two new fixed structural families tested on fresh histories; failure refutes the declared complete programs, not all conceivable tables/objectives/optimizers. No larger family is fitted to residuals.'))


def ets_boundary_analysis(prior=True):
    from .ets_optimizer import objective
    from .ets_round40 import fit_with_trend
    prefix='ets-prior-boundaries'if prior else'ets-boundaries'
    design=ets.read(OUT/(prefix+'-frozen.json'));mapped=ets.read(OUT/(prefix+'-result.json'))
    pool=ets.read(OUT/'ets-search-frozen.json')['tables']['nonreset25']
    rows=[];constraints=[];fractional=[]
    fractions={.166667:1/6,.333333:1/3,.666667:2/3,.833333:5/6,.0833333:1/12}
    for path in mapped['paths']:
        edge=path['focused_edge']
        if edge is None:continue
        z=design['paths'][path['path']]
        right=z['right']if not path.get('anchor')else z['alternatives'][path['anchor']-1]
        observations=dict(path['observations'])
        for index in edge:
            y=vector_path(z['left'],right,index);actual=observations[index];r=[]
            for use_fractions in [False,True]:
                scores=[objective(y,*([fractions.get(v,v)for v in p]if use_fractions else p))for p in pool]
                best=min(scores);minimum=[i for i,v in enumerate(scores)if v==best];compatible=[]
                for i,(a,b)in enumerate(pool):
                    _,_,_,_,reset=fit_with_trend(y,a,b,ets.initial_trend(y))
                    if [bits(a),bits(.001 if reset else b),bits(2.**-52)]==actual:compatible.append(i)
                allowed=sorted(set(minimum)&set(compatible))
                r.append(dict(fractional_search_loads=use_fractions,minimum=minimum,compatible_outputs=compatible,
                    allowed_minimum=allowed,score_gap=None if not compatible else min(scores[i]for i in compatible)-best))
                if not use_fractions and len(allowed)==1:
                    constraints.extend((allowed[0],j)for j in minimum if j!=allowed[0])
            rows.append(dict(path=path['path'],index=index,answer=actual,programs=r))
    # Necessary stable-first constraints from unambiguous ties. Multiple good
    # outputs leave disjunctions unresolved; never manufacture an ordering.
    after={i:set()for i in range(len(pool))}
    for a,b in constraints:after[a].add(b)
    pending=set(after);order=[]
    while pending:
        ready=next((i for i in sorted(pending)if not any(i in after[j]for j in pending)),None)
        if ready is None:break
        pending.remove(ready);order.append(ready)
    save(prefix+'-analysis',dict(rows=rows,pool=pool,stable_first_constraints=constraints,
        necessary_order=order,cycle=bool(pending),
        no_order_can_explain_decimal_score=sum(not r['programs'][0]['allowed_minimum']for r in rows),
        no_order_can_explain_fractional_score=sum(not r['programs'][1]['allowed_minimum']for r in rows),
        limits='Retained boundary replay, not new answers. Exact score equality only. Missing compatible minima refute the complete stated objective/load/reset tree independently of table tie order. Unresolved initialization remains conditional.'))


def ets_stage_analysis(prior=True):
    """Keep optimizer constraints distinct from postfit reset contradictions."""
    prefix='ets-prior-boundaries'if prior else'ets-boundaries'
    path=OUT/(prefix+'-analysis.json');document=ets.read(path);pool=document['pool'];rows=[]
    for r in document['rows']:
        alpha,beta,_=r['answer'];compatible=[i for i,(a,b)in enumerate(pool)
            if bits(a)==alpha and (bits(b)==beta or beta==bits(.001))]
        rows.append(dict(path=r['path'],index=r['index'],answer=r['answer'],
            search_compatible=compatible,
            programs=[dict(fractional_search_loads=p['fractional_search_loads'],
                allowed_search_minimum=sorted(set(compatible)&set(p['minimum'])),
                complete_program_allowed_minimum=p['allowed_minimum'])for p in r['programs']]))
    freeze(prefix+'-stage-analysis-plan',dict(input_sha256=sha(path),
        rule='Nonreset output must name the selected table pair. Beta=.001 permits an unobserved raw beta with the same alpha; it is only a necessary search condition.'))
    save(prefix+'-stage-analysis',dict(rows=rows,input_sha256=sha(path),
        decimal_search_contradictions=sum(not r['programs'][0]['allowed_search_minimum']for r in rows),
        fractional_search_contradictions=sum(not r['programs'][1]['allowed_search_minimum']for r in rows),
        reset_or_other_consumer_contradictions=sum(bool(r['programs'][0]['allowed_search_minimum'])and not r['programs'][0]['complete_program_allowed_minimum']for r in rows),
        limits='Necessary optimizer constraints only. Hidden beta at reset, initialization, exact table-load values and tie order are not recovered by these constraints. No product table is expanded.'))


def normal_shape():
    from fractions import Fraction as Q
    from checker.exact import Span,pull_rounding
    from .normal_quantile import central_graph
    from .x87_programs import Native
    windows=ets.read(OUT/'normal-dense-frozen.json')['windows']
    answers=np.load(OUT/'normal-dense-outputs.npz')['answers'];valid=np.load(OUT/'normal-dense-validity.npz')['numeric']
    collisions=[];collision_count=0
    for i,w in enumerate(windows):
        start=int(w['start'],16);previous=None
        for j in range(w['count']):
            p=value(f'{start+j:016x}');q=bits(p-.5);a=f'{int(answers[i,j]):016x}'
            if previous and valid[i,j]and previous['numeric']and q==previous['q']and a!=previous['answer']:
                collision_count+=1
                if len(collisions)<16:collisions.append(dict(window=i,index=j,p=bits(p),q=q,answer=a,previous=previous))
            previous=dict(p=bits(p),q=q,answer=a,numeric=bool(valid[i,j]))
    linear=[]
    for model in ['strict','stored_pc64']:
        interval=None;failure=None;count=0
        for j,word in enumerate(answers[0]):
            p=value(f'{int(windows[0]["start"],16)+j:016x}');q=abs(Q(p)-Q(.5));z=abs(Q(value(f'{int(word):016x}')))
            if not q:continue
            span=pull_rounding(Span(z,z),model).linear(1/q)
            new=span if interval is None else interval.meet(span);count+=1
            if new is None:
                failure=dict(window=0,index=j,p=bits(p),answer=f'{int(word):016x}',
                    preceding_interval=interval.json(),required_interval=span.json());break
            interval=new
        linear.append(dict(model=model,selections=count,first_empty_intersection=failure,
            limits='Any single real constant C in publication((p-.5)*C), with exact near-half subtraction. Empty intersection refutes the complete multiplication tree, not a numerator-product/denominator tree.'))
    native=Native('/home/ubuntu/codex-runs/round40/x87-round40.so');graph=central_graph()
    raw=np.array([[int(windows[0]['start'],16)+j]for j in range(12)],dtype=np.uint64)
    expected=answers[0,:12];profiles=[]
    for cw in [0x123f,0x133f]:
        for mask in range(1<<graph['sites']):
            pred=native.grid(graph,raw,cw,mask)
            for output in [0,1]:
                if np.array_equal(pred[:,output],expected):profiles.append(dict(cw=cw,mask=mask,output=output))
    rng=random.Random(410203);points=[value(f'{int(bits(.5),16)+rng.choice([-1,1])*rng.randrange(1<<16,1<<18):016x}')for _ in range(256)]
    rows=np.array([[int(bits(p),16)]for p in points],dtype=np.uint64)
    predictions={str(i):[f'{int(v):016x}'for v in native.grid(graph,rows,p['cw'],p['mask'])[:,p['output']]]for i,p in enumerate(profiles)}
    freeze('normal-shape',dict(stored_q_collisions=collision_count,collision_examples=collisions,
        arbitrary_constant_product=linear,profiles=profiles,graph=graph,points=points,predictions=predictions,
        native_sha256=sha(native.path),exact_source_sha256=sha(ROOT/'checker/exact.py'),seed=410203,
        method='Mapping collisions refute functions solely of stored q. Exact publication preimages test one arbitrary linear coefficient. Fixed published AS241 instruction/store programs then receive independent near-half questions; no coefficient or degree is fitted.'))
    c,s=service();ans=s.ask_rows('NORM.S.INV',[[c.number(p)]for p in points],chunk_size=99)
    actual=[a.get('bits')for a in ans]
    save('normal-shape-result',dict(answers=ans,receipts=s.receipts,
        programs={i:dict(first_contradiction=next((j for j,(p,a)in enumerate(zip(pred,actual))if p!=a),None),
            contradictions=sum(p!=a for p,a in zip(pred,actual)))for i,pred in predictions.items()},
        limits='A surviving near-half program is only local evidence with published conditional coefficients. R40 general-central and lower counterexamples remain unchanged; no new normal kernel is installed.'))


def normal_expanded():
    from .normal_quantile import central_graph,expanded_central_coefficients
    from .x87_programs import Native
    native=Native('/home/ubuntu/codex-runs/round40/x87-round40.so')
    windows=ets.read(OUT/'normal-dense-frozen.json')['windows'];old=np.load(OUT/'normal-dense-outputs.npz')['answers']
    raw=np.array([[int(windows[0]['start'],16)+j]for j in range(12)],dtype=np.uint64)
    programs=[];tables={}
    for decimal_tables in [True,False]:
        for decimal_shift in [True,False]:
            label=str((decimal_tables,decimal_shift));coefficients=expanded_central_coefficients(decimal_tables,decimal_shift)
            tables[label]=coefficients;graph=central_graph(*coefficients,squared_coordinate=True)
            for cw in [0x123f,0x133f]:
                for mask in range(1<<graph['sites']):
                    pred=native.grid(graph,raw,cw,mask)
                    for output in [0,1]:
                        if np.array_equal(pred[:,output],old[0,:12]):programs.append(dict(table=label,cw=cw,mask=mask,output=output))
    rng=random.Random(410204);points=[rng.uniform(.0751,.9249)for _ in range(256)];predictions={}
    for i,p in enumerate(programs):
        graph=central_graph(*tables[p['table']],squared_coordinate=True)
        arr=np.array([[int(bits(v),16)]for v in points],dtype=np.uint64)
        predictions[str(i)]=[f'{int(v):016x}'for v in native.grid(graph,arr,p['cw'],p['mask'])[:,p['output']]]
    freeze('normal-expanded',dict(tables=tables,profiles=programs,points=points,predictions=predictions,
        method='Fixed degree7/7 AS241 polynomials transformed algebraically to the q-squared basis. Four explicit source-load interpretations; no fitted coefficient, degree or new piece. Near-half discovery precedes any new central answer.'))
    c,s=service();ans=s.ask_rows('NORM.S.INV',[[c.number(p)]for p in points],chunk_size=99)if predictions else[]
    actual=[a.get('bits')for a in ans]
    save('normal-expanded-result',dict(answers=ans,receipts=s.receipts,
        programs={i:dict(contradictions=sum(p!=a for p,a in zip(pred,actual)),first_contradiction=next((j for j,(p,a)in enumerate(zip(pred,actual))if p!=a),None))for i,pred in predictions.items()},
        stop='An empty near-half profile set refutes this fixed expanded-table family before fresh validation; no questions or fit are claimed in that case.'))


def final_audit():
    """Reconcile the published census against original replies and exact inputs."""
    import struct
    from . import ets_round40
    prefix='ets-census'if (OUT/'ets-census-result.json').exists()else'ets-census-cached-prefix'
    document=ets.read(OUT/(prefix+'-result.json'));design=ets.read(OUT/'ets-census-frozen.json')
    prior=ets.read(design['prior'])
    if sha(design['prior'])!=design['prior_sha256']:raise ValueError('changed prior census')
    expected=np.empty((design['new_histories'],32),dtype=np.uint8)
    def digest(y):return hashlib.sha256(bytes([len(y)])+struct.pack('>'+str(len(y))+'d',*y)).digest()
    wanted={q[1]for q in document['census']['first_questions'].values()if q[0]=='round41'};representatives={}
    for z in histories(design['new_histories'],design['seed']):
        expected[z['index']]=np.frombuffer(digest(z['y']),dtype=np.uint8)
        if z['index']in wanted:representatives[z['index']]=z
    published=np.load(OUT/(prefix+'-outputs.npz'))['answers']
    if published.shape[1:]!=(3,)or len(published)>design['new_histories']:raise ValueError('published census shape differs')
    if len(published)!=design['new_histories']and not document.get('acquisition_interrupted'):raise ValueError('unacknowledged acquisition shortfall')
    positions=defaultdict(list);offset=0
    for r in document['receipts']:
        if r['function']!='FORECAST.ETS.STAT' or r['calls']%3:raise ValueError('unexpected census receipt')
        positions[r['file']].append((offset,r));offset+=r['calls']//3
    if offset!=len(published):raise ValueError('census receipt population differs')
    visited=np.zeros(len(published),dtype=bool);new_types=Counter();seen_files=set()
    timelines={n:dict(type='array',values=[[dict(type='number',bits=bits(j+1.))]for j in range(n)])for n in range(8,33)}
    def visit(path,raw,record):
        if path.name not in positions:return
        request,reply=raw['request'],raw['reply'];rows=request['rows'];ans=reply['rows']
        if request['function']!='FORECAST.ETS.STAT' or len(rows)%3:raise ValueError('changed census request')
        seen_files.add(path.name)
        for start,receipt in positions[path.name]:
            if record['sha256']!=receipt['sha256'] or len(rows)!=receipt['calls']:raise ValueError('changed original census reply')
            for j in range(0,len(rows),3):
                index=start+j//3
                if visited[index]:raise ValueError('duplicate logical census position')
                y=rows[j][0]
                if y.get('type')!='array' or not all(len(r)==1 and r[0].get('type')=='number'for r in y['values']):raise ValueError('changed census shape')
                n=len(y['values']);raw_y=b''.join(bytes.fromhex(r[0]['bits'])for r in y['values'])
                actual=hashlib.sha256(bytes([n])+raw_y).digest()
                if not np.array_equal(np.frombuffer(actual,dtype=np.uint8),expected[index]):raise ValueError('census numeric history changed')
                for k in range(3):
                    if rows[j+k]!=[y,timelines[n],dict(type='number',bits=bits(k+1.)),dict(type='number',bits=bits(0.))]:raise ValueError('changed timeline/statistic/season')
                    a=ans[j+k];new_types[a['type']]+=1
                    if a['type']!='number':raise ValueError('a nonnumeric census component cannot be called a parameter triple')
                    if int(a['bits'],16)!=int(published[index,k]):raise ValueError('NPZ differs from original census answer')
                visited[index]=True
    # Reuse the existing one transport audit and slow-request ledger route.
    ets_round40.OUT=OUT;ets_round40.audit_evidence(allow_unread_inputs=True,visitor=visit)
    if not visited.all() or set(positions)!=seen_files:raise ValueError('incomplete immutable census reconciliation')
    counts=Counter();first={};discoveries=[];families=defaultdict(Counter)
    rows=[tuple(st[str(k)]['bits']for k in [1,2,3])for st in ets.grouped(prior)]
    def observe(question,answer,position):
        key=json.dumps(answer,separators=(',',':'))
        if key not in counts:first[key]=question;discoveries.append(position)
        counts[key]+=1
    for i,answer in enumerate(rows):observe(['round40',i],answer,i+1)
    for i,row in enumerate(published):
        answer=[f'{int(v):016x}'for v in row];observe(['round41',i],answer,10000+i+1)
        families[str(i%8)][json.dumps(answer)]+=1
    report=document['census']
    if report['observations']!=len(rows)+len(published)or report['distinct']!=len(counts):raise ValueError('census population/distinct fields differ')
    if dict(counts)!=report['counts'] or first!=report['first_questions']or discoveries!=report['discovery_positions']:raise ValueError('census count/discovery map differs')
    if {k:dict(v)for k,v in families.items()}!=document['families']:raise ValueError('census family map differs')
    cases={tuple(json.loads(key)):z for key,z in document['first_cases'].items()}
    if set(cases)!={tuple(json.loads(key))for key in first}:raise ValueError('first-case keys differ')
    for key,q in first.items():
        actual=cases[tuple(json.loads(key))];original=prior['cases'][q[1]]if q[0]=='round40'else representatives[q[1]]
        if set(actual)!=set(original)or digest(actual['y'])!=digest(original['y'])or any(actual[k]!=original[k]for k in original if k!='y'):raise ValueError('first-case representative changed')  # ubs:ignore -- public numeric history digests and fields, no secrets
    save('final-census-reconciliation',dict(histories=report['observations'],new_histories=len(published),prior_histories=len(rows),
        parameter_triples=len(counts),new_component_types=new_types,original_reply_files=len(seen_files),
        exact_history_input_digests=True,reply_bits_equal_npz=True,count_discovery_family_maps_equal=True,first_cases_equal=True,
        requested_histories=design['new_histories']+10000,acquisition_interrupted=document.get('acquisition_interrupted',False),
        arrays_sha256=sha(OUT/(prefix+'-outputs.npz')),source_sha256=sha(Path(__file__)),
        limits='Rebuilt from immutable original reply contents and deterministic exact binary64 inputs, independently of mutable checkpoint counts. Prior 10,000 histories are bound to the original R40 document hash.'))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('step',choices=['census','resume-census','census-cached-prefix','ets-boundaries','ets-prior-boundaries','ets-search','ets-prior-analysis','ets-boundary-analysis','ets-prior-stage','ets-stage','normal','resume-normal','normal-shape','normal-expanded','regression','replay','residual-isolation','residual-sums','residual-fresh','residual-additions','binomial','binomial-fresh','yield-plateaus','binomial-beta','final-audit']);parser.add_argument('--count',type=int,default=990000);parser.add_argument('--batch',type=int,default=2000)
    args=parser.parse_args();OUT.mkdir(exist_ok=True)
    if args.step=='census':census(args.count,args.batch)
    elif args.step=='resume-census':census(args.count,args.batch,True)
    elif args.step=='census-cached-prefix':census(args.count,args.batch,True,True)
    elif args.step=='ets-boundaries':ets_boundaries()
    elif args.step=='ets-prior-boundaries':ets_boundaries(True)
    elif args.step=='ets-search':ets_search()
    elif args.step=='ets-prior-analysis':ets_boundary_analysis(True)
    elif args.step=='ets-boundary-analysis':ets_boundary_analysis(False)
    elif args.step=='ets-prior-stage':ets_stage_analysis(True)
    elif args.step=='ets-stage':ets_stage_analysis(False)
    elif args.step=='normal-shape':normal_shape()
    elif args.step=='normal-expanded':normal_expanded()
    elif args.step=='final-audit':final_audit()
    elif args.step=='regression':regression()
    elif args.step=='replay':replay()
    elif args.step=='residual-isolation':residual_isolation()
    elif args.step=='residual-sums':residual_sums()
    elif args.step=='residual-fresh':residual_fresh()
    elif args.step=='residual-additions':residual_additions()
    elif args.step=='binomial':binomial()
    elif args.step=='binomial-fresh':binomial_fresh()
    elif args.step=='yield-plateaus':yield_plateaus()
    elif args.step=='binomial-beta':binomial_beta()
    else:normal(args.step=='resume-normal')

if __name__=='__main__':main()
