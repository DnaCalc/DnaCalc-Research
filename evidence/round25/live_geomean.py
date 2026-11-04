"""Read exact scalar ranges from the existing packet; ask the canonical worker."""
import hashlib
import json
from pathlib import Path
import re
import sys
from live_variance import Client, number

packet, output = map(Path, sys.argv[1:3])
output.mkdir(exist_ok=False)
source = packet/'workbook-questions.json'
questions = json.loads(source.read_text())['questions']
echoes = {q['row']:q['intended_bits'].removeprefix('0x') for q in questions if q['role']=='operand_echo'}
selected = []
for q in questions:
    match = re.fullmatch(r'GEOMEAN\(D(\d+):D(\d+)\)',q['formula'])
    if match:
        first,last=map(int,match.groups())
        selected.append(dict(q,args=[echoes[r] for r in range(first,last+1)]))
(output/'questions.json').write_text(json.dumps(selected,indent=2)+'\n')
reply=Client(timeout=900).ask('GEOMEAN',[[number(x)for x in q['args']]for q in selected])
(output/'reply.json').write_text(json.dumps(reply,indent=2)+'\n')
assert reply['inputs']==[[number(x)for x in q['args']]for q in selected]
binding={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest()for p in [Path(__file__),source,Path(__file__).with_name('live_variance.py')]}
(output/'binding.json').write_text(json.dumps(binding,indent=2)+'\n')
print('GEOMEAN',len(selected),'observed',reply['ok'])
