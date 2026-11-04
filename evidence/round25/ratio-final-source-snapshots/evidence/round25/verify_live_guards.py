"""Independent rational replay of retained live scalar observations."""
from collections import Counter
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import struct
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'checker'))
from variance import arithmetic,evaluate
from normal import encoding,sqrt_round

raw,output=map(Path,sys.argv[1:3]);files=set();readings=[]


def read(path):
    files.add(path.resolve());return json.loads(path.read_text())


def accept(function,questions,reply):
    assert reply['ok']and len(questions)==len(reply['rows'])==len(reply['inputs'])
    for q,a,echo in zip(questions,reply['rows'],reply['inputs']):
        args=q.get('service_args')
        if args is None:args=[dict(type='number',bits=x.removeprefix('0x'))for x in q['args']]
        assert echo==[{k:v for k,v in x.items()if k!='mode'}for x in args]
        readings.append(dict(function=function,args=args,answer=a,request_id=reply['id'],excel_build=reply['excel_build']))


for packet in read(raw/'live-variance/questions.json'):
    accept(packet['function'],packet['rows'],read(raw/'live-variance'/(packet['function']+'.json')))
for directory in ['live-variance-edges','live-variance-relative','live-variance-cancellation','live-variance-cancellation-v2','live-variance-sqrt','live-bessel','live-discrete']:
    for path in sorted((raw/directory).glob('*.json')):
        if path.name in ['binding.json','states.json','summary.json']:continue
        v=read(path)
        if 'reply'not in v:continue
        function=path.stem.split('-')[-1]
        accept(function,v['questions'],v['reply'])
for directory in ['live-guard-operations','live-relative-operations','live-unique-relative-operations']:
    q=read(raw/directory/'questions.json');accept('VAR.S',q['questions'],read(raw/directory/'reply.json'))
for path in sorted((raw/'live-variance-context').glob('*.json')):
    v=read(path);accept(path.stem,[dict(args=x)for x in v['args']],v['reply'])
accept('GEOMEAN',read(raw/'live-geomean/questions.json'),read(raw/'live-geomean/reply.json'))

counts=Counter();identity=set();variance_counts=Counter();root_checks=[]
for reading in readings:
    function=reading['function'];counts[function]+=1
    identity.add(json.dumps([function,reading['args']],sort_keys=True))
    if function not in ['VAR','VAR.S','VAR.P','STDEV.S','STDEV.P']:continue
    xs=[Q(struct.unpack('>d',bytes.fromhex(x['bits']))[0])for x in reading['args']]
    assert all(x['type']=='number'for x in reading['args'])
    op=lambda k,a,b:arithmetic(k,'stored_pc64',a,b)
    s=q=Q(0)
    for x in xs:s=op('add',s,x);q=op('add',q,op('mul',x,x))
    n=Q(len(xs));square=op('mul',s,s);mean=op('div',s,n)
    pop=function.endswith('.P');vf='VAR.P'if pop else'VAR.S'
    profile=dict(algorithm=8 if pop else 7,linear_order=0,square_order=0,models=['stored_pc64']*5)
    variance=evaluate(profile,xs,vf)
    cancellation=op('div',op('div',square,n),q)>Q(.99)if q else square>0
    magnitude=op('mul',variance,Q(1000))<=abs(mean)
    if cancellation or magnitude:
        profile['algorithm']=0;variance=evaluate(profile,xs,vf)
    calculated=sqrt_round(variance,'stored_pc64')if function.startswith('STDEV')else variance
    got=encoding(calculated);expected='0x'+reading['answer']['bits']
    assert got==expected,(function,reading,got)
    variance_counts['readings']+=1;variance_counts['matched']+=1
    variance_counts['negative_mean']+=mean<0
    variance_counts['cancellation']+=cancellation;variance_counts['magnitude']+=magnitude
for path in sorted((raw/'live-variance-sqrt').glob('STDEV.*.json')):
    v=read(path)
    for q,a in zip(v['questions'],v['reply']['rows']):
        variance=Q(struct.unpack('>d',bytes.fromhex(q['variance'].removeprefix('0x')))[0])
        staged=encoding(sqrt_round(variance,'stored_pc64'));strict=encoding(sqrt_round(variance,'strict'))
        assert staged==q['staged_root']=='0x'+a['bits'] and strict==q['strict_root']!=staged
        root_checks.append(dict(id=q['id'],variance=q['variance'],strict=strict,staged=staged))
sources=[Path(__file__),ROOT/'checker/variance.py',ROOT/'checker/exact.py',ROOT/'checker/normal.py']
files.update(p.resolve()for p in sources)
result=dict(complete=True,route='Separate Value2 cell references; scalarized questions are not workbook/range captures.',
    live_readings=len(readings),unique_function_argument_modes=len(identity),functions=counts,
    variance_replay=variance_counts,independent_sqrt_checks=root_checks,
    inputs_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in sorted(files)},
    limits='Observed scalar scope only. Stored(v*1000)<=abs(mean) and stored(stored(S*S/n)/Q)>.99 fit. The16 later unique constructions distinguish the ratio from both multiplied and correction-divided forms. The mean/1000 magnitude rearrangement remains indistinguishable. Counts are diagnostics, not proof over the parameter space.')
with output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({k:v for k,v in result.items()if k not in ['inputs_sha256','independent_sqrt_checks']}))
