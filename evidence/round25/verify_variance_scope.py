"""Compare complete canonical identities, not positional bits across changed scopes."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import sys

executable,checkout,before,output=map(lambda p:Path(p).resolve(),sys.argv[1:5])
output.mkdir(exist_ok=False)
functions=['VAR','VAR.S','VAR.P','VARA','VARP','VARPA','STDEV','STDEV.S','STDEV.P','STDEVA','STDEVP','STDEVPA']
sources=[Path(__file__),executable,before]
binding={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in sources}
after=[]
for function in functions:
    trace=output/(function+'.jsonl')
    with(output/(function+'.log')).open('x')as log:
        subprocess.run([str(executable),'witness','judge','--bank','full','--func',function,'--trace',str(trace)],cwd=checkout,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=180)  # ubs:ignore — operator selects the local canonical judge; explicit argv, owned output and bounded runtime
    after.extend(json.loads(line)for line in trace.read_text().splitlines())
old=[json.loads(line)for line in before.read_text().splitlines()]
def key(row):return json.dumps({k:v for k,v in row.items()if k!='actual'},sort_keys=True,separators=(',',':'))
a={key(row):row for row in old};b={key(row):row for row in after}
assert len(a)==len(old)and len(b)==len(after),'duplicate complete canonical identity'
assert a.keys()==b.keys(),'input/answer/environment/source identity changed'
counts=Counter()
for identity,x in a.items():
    y=b[identity];was=x['actual']==x['expected'];now=y['actual']==y['expected']
    for name,flag in [('rows',True),('before_matched',was),('after_matched',now),('gains',now and not was),('losses',was and not now),('changed_outputs',x['actual']!=y['actual'])]:counts[name]+=flag
assert binding=={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in sources}
(output/'after.jsonl').write_text(''.join(json.dumps(b[k],sort_keys=True)+'\n'for k in sorted(b)))
result=dict(complete=True,counts=counts,identity_answer_environment_source_sha256=hashlib.sha256(('\n'.join(sorted(a))+'\n').encode()).hexdigest(),inputs_sha256=binding,
    limits='All twelve original function cohorts, including nonnumeric/error controls, through the existing canonical judge with source invocation context. Exact identity equality is required. No changing global denominator is used for gains or losses.')
(output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
