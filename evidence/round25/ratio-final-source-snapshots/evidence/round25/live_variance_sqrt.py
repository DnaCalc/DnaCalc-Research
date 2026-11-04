"""Observe the fixed square-root store discriminants through scalar cells."""
import hashlib
import json
from pathlib import Path
import sys
from live_variance import Client,number

source,output=map(Path,sys.argv[1:3]);output.mkdir(exist_ok=False)
questions=json.loads(source.read_text())['questions']
for function in ('STDEV.S','STDEV.P'):
    selected=[q for q in questions if q['function']==function]
    reply=Client(timeout=900).ask(function,[[number(x)for x in q['args']]for q in selected])
    (output/(function+'.json')).write_text(json.dumps(dict(questions=selected,reply=reply),indent=2)+'\n')
    assert reply['inputs']==[[number(x)for x in q['args']]for q in selected]
    print(function,len(selected),'observed',reply['ok'],flush=True)
(output/'binding.json').write_text(json.dumps({str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest()for p in[Path(__file__),source,Path(__file__).with_name('live_variance.py')]},indent=2)+'\n')
