"""Locate cancellation switches independently of the magnitude guard."""
import hashlib
import json
from pathlib import Path
import sys
from fractions import Fraction as Q
from live_variance_edges import bits,value
from live_variance import Client,number,evaluate,ROOT
from normal import encoding

output=Path(sys.argv[1]);output.mkdir(exist_ok=False)
states=[dict(function=f,exponent=e,low=bits(1.12),high=bits(1.13),active=True,history=[])
        for f in ('VAR.S','VAR.P')for e in (0,20)]


def construct(state,t):
    x=value(t);xs=[(2-x)*2.0**state['exponent'],2.0**state['exponent'],x*2.0**state['exponent']]
    result={}
    for name,algorithm in [('raw',7 if state['function']=='VAR.S'else 8),('centered',0)]:
        result[name]=encoding(evaluate(dict(algorithm=algorithm,linear_order=0,square_order=0,models=['stored_pc64']*5),list(map(Q,xs)),state['function'])).removeprefix('0x')
    return dict(t=f'{t:016x}',args=[f'{bits(x):016x}'for x in xs],predictions=result)


for iteration in range(58):
    for function in ('VAR.S','VAR.P'):
        batch=[]
        for state in states:
            if state['function']!=function or not state['active']:continue
            if iteration==0:
                for t in (state['low'],state['high']):batch.append((state,construct(state,t)))
                continue
            lo,hi=state['low'],state['high']
            if hi-lo<=1:state['active']=False;continue
            middle=(lo+hi)//2;selected=None
            for offset in [0]+[k for d in range(1,65)for k in(-d,d)]:
                t=middle+offset
                if not lo<t<hi:continue
                row=construct(state,t)
                if len(set(row['predictions'].values()))==2:selected=row;break
            if selected is None:state['active']=False;state['limit']='indistinguishable interior';continue
            batch.append((state,selected))
        if not batch:continue
        reply=Client(timeout=900).ask(function,[[number(x)for x in r['args']]for _,r in batch])
        (output/f'{iteration:02d}-{function}.json').write_text(json.dumps(dict(questions=[r for _,r in batch],reply=reply),indent=2)+'\n')
        assert reply['inputs']==[[number(x)for x in r['args']]for _,r in batch]
        for (state,row),answer in zip(batch,reply['rows']):
            matches=[k for k,v in row['predictions'].items()if v==answer.get('bits')]
            state['history'].append(dict(row,answer=answer))
            if len(matches)!=1:state['active']=False;state['limit']='neither or indistinguishable';continue
            key='low'if matches[0]=='centered'else'high';t=int(row['t'],16)
            if iteration==0:assert state[key]==t,'measured bracket endpoints do not straddle'
            else:state[key]=t
    (output/'states.json').write_text(json.dumps(states,indent=2)+'\n')
    print(iteration,[s['high']-s['low']for s in states],flush=True)
    if not any(s['active']for s in states):break
sources=[Path(__file__),Path(__file__).with_name('live_variance_edges.py'),Path(__file__).with_name('live_variance.py'),ROOT/'checker/variance.py',ROOT/'checker/exact.py',ROOT/'checker/normal.py']
(output/'binding.json').write_text(json.dumps({str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest()for p in sources},indent=2)+'\n')
