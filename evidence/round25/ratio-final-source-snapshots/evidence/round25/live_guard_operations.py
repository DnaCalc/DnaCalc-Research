"""Observe complete-result guard discriminants; retain original unknowns."""
import hashlib
import json
from pathlib import Path
import sys
from live_variance import Client,number

source,output=map(Path,sys.argv[1:3]);output.mkdir(exist_ok=False)
questions=json.loads(source.read_text())['questions']
reply=Client(timeout=900).ask('VAR.S',[[number(x)for x in q['args']]for q in questions])
(output/'reply.json').write_text(json.dumps(reply,indent=2)+'\n')
assert reply['inputs']==[[number(x)for x in q['args']]for q in questions]
(output/'questions.json').write_bytes(source.read_bytes())
(output/'binding.json').write_text(json.dumps({str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest()for p in[Path(__file__),source,Path(__file__).with_name('live_variance.py')]},indent=2)+'\n')
print('VAR.S',len(questions),'observed',reply['ok'],flush=True)
