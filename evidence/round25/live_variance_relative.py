"""Ask the packet's adjacent cancellation-edge questions, using exact echoes."""
import hashlib
import json
from pathlib import Path
import re
import sys
from live_variance import Client, number

packet, output=map(Path,sys.argv[1:3]);output.mkdir(exist_ok=False)
source=packet/'workbook-questions.json';questions=json.loads(source.read_text())['questions']
echoes={q['row']:q['intended_bits'].removeprefix('0x')for q in questions if q['role']=='operand_echo'}
for function in ('VAR.S','VAR.P','STDEV.S','STDEV.P'):
    selected=[]
    for q in questions:
        match=re.fullmatch(r'_xlfn\.'+re.escape(function)+r'\(D(\d+):D(\d+)\)',q['formula'])
        if match and 'candidate_predictions'in q:
            first,last=map(int,match.groups());selected.append(dict(q,args=[echoes[r]for r in range(first,last+1)]))
    reply=Client(timeout=900).ask(function,[[number(x)for x in q['args']]for q in selected])
    (output/(function+'.json')).write_text(json.dumps(dict(questions=selected,reply=reply),indent=2)+'\n')
    assert reply['inputs']==[[number(x)for x in q['args']]for q in selected]
    print(function,len(selected),'observed',reply['ok'],flush=True)
binding={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest()for p in [Path(__file__),source,Path(__file__).with_name('live_variance.py')]}
(output/'binding.json').write_text(json.dumps(binding,indent=2)+'\n')
