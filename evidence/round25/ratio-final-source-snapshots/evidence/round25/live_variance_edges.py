"""Bisect observed raw/centered decisions, retaining indistinguishable gaps."""
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import struct
import sys
from fractions import Fraction as Q

from live_variance import Client, number, evaluate, ROOT
from normal import encoding


def bits(x):
    return struct.unpack('>Q', struct.pack('>d', x))[0]


def value(x):
    return struct.unpack('>d', struct.pack('>Q', x))[0]


def construct(shape, scale, function):
    xs = [x*value(scale) for x in shape]
    predictions = {}
    for name, algorithm in [('raw', 7 if function == 'VAR.S' else 8), ('centered', 0)]:
        z = evaluate(dict(algorithm=algorithm,linear_order=0,square_order=0,models=['stored_pc64']*5),
                     [Q(x) for x in xs], function)
        predictions[name] = encoding(z).removeprefix('0x')
    return dict(scale=f'{scale:016x}',args=[f'{bits(x):016x}'for x in xs],predictions=predictions)


def main():
    coarse, output = map(Path, sys.argv[1:3])
    output.mkdir(exist_ok=False)
    sources=[Path(__file__),Path(__file__).with_name('live_variance.py'),ROOT/'checker/variance.py',
             ROOT/'checker/exact.py',ROOT/'checker/normal.py']
    binding={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest()for p in sources}
    questions=json.loads((coarse/'questions.json').read_text())
    states=[]
    for packet in questions:
        function=packet['function'];reply=json.loads((coarse/(function+'.json')).read_text())
        assert reply['ok'] and reply['inputs']==[[number(x)for x in r['args']]for r in packet['rows']]
        for n in (2,3,4,8):
            rows=[(r,a)for r,a in zip(packet['rows'],reply['rows'])if r['n']==n]
            for i,(r,a)in enumerate(rows):
                if '0x'+a.get('bits','')==r['predictions']['raw_pc64']:
                    low,low_answer=rows[i-1]
                    assert '0x'+low_answer['bits']==low['predictions']['centered_pc64']
                    shape=[math.ldexp(value(int(x,16)),-r['exponent'])for x in r['args']]
                    states.append(dict(function=function,n=n,shape=shape,
                        low=bits(2.0**low['exponent']),high=bits(2.0**r['exponent']),
                        coarse_low=low,coarse_high=r,active=True,history=[]))
                    break
    for iteration in range(58):
        for function in ('VAR.S','VAR.P'):
            batch=[]
            for state in states:
                if state['function']!=function or not state['active']:continue
                lo,hi=state['low'],state['high']
                if hi-lo<=1:state['active']=False;continue
                middle=(lo+hi)//2;selected=None
                for offset in [0]+[k for d in range(1,65)for k in (-d,d)]:
                    coordinate=middle+offset
                    if not lo<coordinate<hi:continue
                    row=construct(state['shape'],coordinate,function)
                    if len(set(row['predictions'].values()))==2:
                        selected=(state,row);break
                if selected is None:
                    state['active']=False;state['limit']='No differing predictions inside final bracket';continue
                batch.append(selected)
            if not batch:continue
            reply=Client(timeout=900).ask(function,[[number(x)for x in r['args']]for _,r in batch])
            (output/f'{iteration:02d}-{function}.json').write_text(json.dumps(dict(questions=[r for _,r in batch],reply=reply),indent=2)+'\n')
            assert reply['inputs']==[[number(x)for x in r['args']]for _,r in batch]
            for (state,row),answer in zip(batch,reply['rows']):
                got=answer.get('bits');choices=[k for k,x in row['predictions'].items()if x==got]
                state['history'].append(dict(row,answer=answer))
                if len(choices)!=1:
                    state['active']=False;state['limit']='Excel differs from both declared programs';continue
                key='low'if choices[0]=='centered'else'high';state[key]=int(row['scale'],16)
        (output/'states.json').write_text(json.dumps(states,indent=2)+'\n')
        print('iteration',iteration,'active',sum(s['active']for s in states),
              'gaps',[s['high']-s['low']for s in states],flush=True)
        if not any(s['active']for s in states):break
    assert binding=={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest()for p in sources}
    (output/'binding.json').write_text(json.dumps(binding,indent=2)+'\n')


if __name__=='__main__':main()
