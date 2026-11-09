"""Frozen ETS expression predictions, fresh service answers and exact comparison.

Original MIT research. A conditional proof names its measured parameter inputs.
Every frozen document carries its producer source, not only its digest.
"""
import argparse
from collections import Counter
import json
import math
import random
from pathlib import Path
from .ets import (OUT, read, write, sha, bits, value, service, ask_cases, grouped,
                  numeric_fit, numeric_forecast, numeric_statistics, initial_trend,
                  slope, SLOPES, sequential_total)

def freeze(name,document):
    p=OUT/(name+'-frozen.json')
    document['sources']={str(Path(v).name):Path(v).read_text()for v in [__file__,str(Path(__file__).with_name('ets.py'))]}
    document['producer_sha256']=sha(__file__)
    if p.exists():
        if read(p)['predictions']!=document['predictions']:raise ValueError('changed frozen predictions')
    else:write(p,document)

def finish(name,cases,predictions,function,variants):
    c,s=service();d=ask_cases(c,s,cases,function,variants)
    write(OUT/(name+'-answers.json'),d)
    actual=[r.get('bits')for r in d['answers']]
    failures=[dict(index=i,predicted=p,actual=a)for i,(p,a)in enumerate(zip(predictions,actual))if p!=a]
    result=dict(cases=len(cases),selections=len(actual),failures=failures,producer_sha256=sha(__file__))
    write(OUT/(name+'-proof.json'),result)
    print(name,result['cases'],result['selections'],'failures',len(failures),failures[:12],flush=True)
    return d

def recurrence():
    name='step2-final-recurrence';rng=random.Random(392771);cases=[]
    for i in range(64):
        n=rng.randrange(7,33);ex=rng.randrange(-200,201);trend=math.ldexp(128.,ex)
        half=[rng.randrange(-8,9)for _ in range((n+1)//2)]
        noise=half+half[-2::-1]if n%2 else half+half[::-1];total=sum(noise)
        c=rng.choice([1/256,1/64,1/16]);r=rng.choice([0.,1/64,1/16,.25])
        y=[math.ldexp(128*j+c*(3*(2*j-n+1)**2-(n*n-1))+r*(n*noise[j]-total),ex)for j in range(n)]
        if any(slope(y,k)!=trend for k in SLOPES):raise ValueError('constructed trend not exact')
        cases.append(dict(y=y,season=0,trend=trend))
    c,s=service();stat=ask_cases(c,s,cases,'FORECAST.ETS.STAT',[1,2,3]);write(OUT/(name+'-parameters.json'),stat)
    leaves=grouped(stat);pred=[]
    for z,p in zip(cases,leaves):
        for h in [0.,1.,2.,5.]:pred.append(bits(numeric_forecast(z['y'],value(p['1']['bits']),value(p['2']['bits']),h,z['trend'])))
    freeze(name,dict(cases=cases,predictions=pred,conditional_inputs='STAT 1-2; exact constructed OLS trend; binary64 innovations tree with two epsilon seasons and final centering.'))
    finish(name,cases,pred,'FORECAST.ETS',[0.,1.,2.,5.])

def three():
    name='step2-three-general';rng=random.Random(393181);cases=[];pred=[]
    for i in range(100):
        ex=rng.randrange(-200,201)
        if i<70:y=[math.ldexp(rng.randrange(-10000,10001)/32,ex)for _ in range(3)]
        else:
            base=math.ldexp(rng.choice([1.,1.5,math.pi]),ex);raw=int(bits(base),16)
            y=[value(f'{raw+rng.randrange(-12,13):016x}')for _ in range(3)]
        cases.append(dict(y=y,season=0))
        t=(y[2]-y[0])/2
        for h in [0.,1.,2.,5.]:pred.append(bits(numeric_forecast(y,.1,.099,h,t)))
    freeze(name,dict(cases=cases,predictions=pred,method='Unconditional three-observation candidate: initial trend (last-first)/2, default .1/.099; includes single-bit neighbours and exponent changes.'))
    finish(name,cases,pred,'FORECAST.ETS',[0.,1.,2.,5.])
    c,s=service();write(OUT/(name+'-stat.json'),ask_cases(c,s,cases,'FORECAST.ETS.STAT',range(1,9)))

def complete(final=False):
    from .ets_optimizer import parameters,corrected
    from .ets_confidence import variance
    name='step45-final-validation'if final else'step45-complete';rng=random.Random(399401 if final else 399193);cases=[];stats=[];forecasts=[];confidence=[];raw=[];z=value(read(OUT/'step6-normal-quantiles.json')['reply']['rows'][0]['bits'])
    for i in range(128):
        n=rng.randrange(8,33);ex=rng.randrange(-350,351)
        if i<100:
            trend=rng.randrange(-32,33)/8
            y=[math.ldexp(rng.randrange(-512,513)/32+trend*j,ex)for j in range(n)]if i<75 else[math.ldexp(j*j/32+trend*j+rng.randrange(-8,9)/64,ex)for j in range(n)]
        else:
            base=math.ldexp(rng.choice([1.,1.5,math.pi]),ex);ib=int(bits(base),16)
            y=[value(f'{ib+rng.randrange(-16,17):016x}')for _ in range(n)]
        cases.append(dict(y=y,season=0));a,b=parameters(y);l,t,se,pb,e,p,corr=corrected(y,a,b);raw.append(dict(alpha=a,beta=b,correction=corr))
        stats.extend(bits(v)for v in [a,pb,2.**-52,*numeric_statistics(y,a,b,initial_trend(y)),1.])
        forecasts.extend(bits((l+h*t)+se[(n-1+int(h))%2])for h in [0.,1.,2.,5.])
        confidence.extend(bits(z*math.sqrt(sequential_total(v*v for v in e)/n*variance(a,pb,h,39)))for h in [1.,2.,5.,10.])
    extra={k:Path(__file__).with_name(k).read_text()for k in ['ets_optimizer.py','ets_confidence.py']}
    for kind,pred in [('stat',stats),('forecast',forecasts),('confidence',confidence)]:freeze(name+'-'+kind,dict(cases=cases,predictions=pred,raw_parameters=raw,extra_sources=extra,ordinary_cases=100,adjacent_cases=28))
    finish(name+'-stat',cases,stats,'FORECAST.ETS.STAT',range(1,9))
    finish(name+'-forecast',cases,forecasts,'FORECAST.ETS',[0.,1.,2.,5.])
    finish(name+'-confidence',cases,confidence,'FORECAST.ETS.CONFINT',[(h,.95)for h in [1.,2.,5.,10.]])

def catalogue(seasonal=False):
    from .corpus_cases import ets_selected_calls
    rng=random.Random(397991 if seasonal else 399540);cases=[]
    for i in range(40 if seasonal else 80):
        if seasonal:
            p=rng.choice([2,4,8]);n=p*rng.randrange(max(2,(8+p-1)//p),32//p+1)
            pattern=[rng.randrange(-64,65)/8 for _ in range(p)]
            drift=rng.randrange(-16,17)/16;offset=rng.randrange(-128,129)/4
            y=[offset+pattern[j%p]+drift*j for j in range(n)]
            kind='exact_seasonal_affine'
        else:
            p=0;n=rng.randrange(8,33)if i<64 else 3+(i-64)%5
            trend=rng.randrange(-16,17)/8
            y=[rng.randrange(-128,129)/8+trend*j for j in range(n)]if i%2 else[j*j/16+trend*j+rng.randrange(-8,9)/16 for j in range(n)]
            kind='random'if i%2 else'quadratic'
        cases.append(dict(y=y,timeline=list(range(1,n+1)),seasonality=p,data_completion=1,aggregation=1,pattern=kind))
    prefix='ets-seasonal-catalogue'if seasonal else'ets-catalogue'
    doc=dict(sheet='dcj9k_r39_seasonal'if seasonal else'dcj9k_r39_ets',expected_cases=len(cases),statistics=list(range(1,9)),evidence_round='round39',cases=cases,producer_source=Path(__file__).read_text())
    if seasonal:doc['evidence_section']='ets-seasonal-catalogue'
    design=OUT/(prefix+'-design.json')
    if design.exists():
        if read(design)['cases']!=cases:raise ValueError('catalogue design changed')
    else:write(design,doc)
    c,s=service();groups={};labels={}
    for i,z in enumerate(cases):
        for name,args in ets_selected_calls(z,doc['statistics'],c.number,c.number_array):groups.setdefault(name,[]).append(args);labels.setdefault(name,[]).append(i)
    for name,rows in groups.items():
        answers=s.ask_rows(name,rows,chunk_size=100);write(OUT/(prefix+'-'+name+'.json'),dict(labels=labels[name],answers=answers,receipts=s.receipts));print('catalogue',name,len(answers),flush=True)

def application_replay(executable,output,include_round38=False,directories=None,functions=None):
    """Continue R38's registered-route replay over original audited questions."""
    import gzip
    import subprocess
    from collections import defaultdict
    from .fingerprint import accept
    executable=Path(executable).resolve();output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False)
    provenance=json.loads(subprocess.check_output([str(executable),'provenance'],timeout=30))
    directories=[OUT/'service-cache']if directories is None else list(map(Path,directories))
    if include_round38:directories+=sorted(p for p in OUT.parent.joinpath('round38').glob('*cache')if p.is_dir())
    groups=defaultdict(list);bindings={};seen=set()
    def decode(a):
        if a['type']=='number':return '0x'+a['bits']
        if a['type']=='array':return [[decode(v)for v in row]for row in a['values']]
        if a['type']=='boolean':return a['value']
        raise ValueError('unsupported audited ETS operand')
    def expected(a):
        if a['type']=='number':return '0x'+a['bits']
        if a['type']=='error':return 'error:'+{'#NUM!':'Num','#N/A':'NA','#DIV/0!':'Div0','#VALUE!':'Value','#REF!':'Ref','#NAME?':'Name','#SPILL!':'Spill'}[a['code']]
        raise ValueError('unsupported ETS answer')
    for directory in directories:
        for path in sorted(directory.glob('*.json.gz')):
            if '.refusal-'in path.name:continue
            d=read(path);request,reply=d['request'],d['reply'];name=request['function']
            selected = name in functions if functions is not None else name.startswith('FORECAST.ETS')
            if not selected or not reply.get('ok')or reply['id']in seen:continue
            accept(reply,name,request['rows'],audit=d.get('audit'),result_index=request.get('result_index'));seen.add(reply['id']);source=str(path.relative_to(OUT.parent));bindings[source]=sha(path)
            for i,(args,answer)in enumerate(zip(request['rows'],reply['rows'])):
                row=dict(id=reply['id']+':'+str(i),args=[decode(a)for a in args],expected=expected(answer),source=source)
                if 'result_index' in request:row['result_index']=request['result_index']
                groups[name].append(row)
    counts={};contradictions=[]
    for name,rows in groups.items():
        q=output/(name+'.inputs.json');p=output/(name+'.predictions.jsonl');write(q,dict(function=name,witnesses=rows))
        subprocess.run([str(executable),'invoke-inputs',str(q),str(p)],check=True,timeout=600)
        actual=[json.loads(line)for line in p.read_text().splitlines()]
        if len(actual)!=len(rows):raise ValueError('missing registered-route selections')
        tally=Counter();by_round=defaultdict(Counter);shown=0
        for row,a in zip(rows,actual):
            if row['id']!=a['id']or row['args']!=a['args']:raise ValueError('changed registered-route arguments')
            prediction=a['calc_func_prediction'];kind='matched'if prediction==row['expected']else'refused'if prediction.startswith('fault:')and'NotHeld'in prediction else'failed'if prediction.startswith('fault:')else'mismatched'
            tally[kind]+=1;by_round[row['source'].split('/')[0]][kind]+=1
            if kind=='mismatched'and shown<8:contradictions.append(dict(function=name,**row,prediction=prediction));shown+=1
        counts[name]=dict(total=len(rows),counts=tally,by_round=by_round,pass_rate=tally['matched']/len(rows))
        for path in[q,p]:
            path.with_suffix(path.suffix+'.gz').write_bytes(gzip.compress(path.read_bytes(),mtime=0));path.unlink()
        print(name,counts[name],flush=True)
    doc=dict(executable_sha256=sha(executable),compiled_provenance=provenance,raw_files=bindings,functions=counts,first_contradictions=contradictions,limits='Registered CalcFunc invocation of the original audited typed operands. Worksheet errors may match; every refusal, other fault and numeric discrepancy fails. Workbook corpus results use the existing canonical corpus command separately.')
    write(output/'summary.json',doc)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['recurrence','three','complete','catalogue','catalogue-seasonal','validate','replay']);p.add_argument('--executable',default='target/debug/symbolic-racer');p.add_argument('--output');p.add_argument('--include-round38',action='store_true');a=p.parse_args()
    if a.command=='recurrence':recurrence()
    elif a.command=='three':three()
    elif a.command=='catalogue':catalogue()
    elif a.command=='catalogue-seasonal':catalogue(True)
    elif a.command=='replay':application_replay(a.executable,a.output or OUT/'application-replay',a.include_round38)
    else:complete(a.command=='validate')
