"""New scalarized Bessel/discrete questions, not captures of range workbooks."""
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import re
import sys
from live_variance import Client,number
from client import boolean

packet,output=map(Path,sys.argv[1:3]);output.mkdir(exist_ok=False)
source=packet/'workbook-questions.json';questions=json.loads(source.read_text())['questions']
echoes={q['row']:q['intended_bits'].removeprefix('0x')for q in questions if q['role']=='operand_echo'}
groups=defaultdict(list)
for q in questions:
    match=re.fullmatch(r'(?:_xlfn\.)?(BESSELJ|BESSELY|BINOM\.DIST|POISSON\.DIST)\(([^()]*)\)',q['formula'])
    if not match:continue
    function,arguments=match.groups();args=[]
    for token in arguments.split(','):
        if re.fullmatch(r'D\d+',token):args.append(number(echoes[int(token[1:])]))
        elif token in ('TRUE','FALSE'):args.append(boolean(token=='TRUE',literal=True))
        elif re.fullmatch(r'-?\d+',token):args.append(number(float(token)))
        else:raise ValueError('unparsed scalar token '+token)
    groups[function].append(dict(q,service_args=args))
for function,selected in groups.items():
    reply=Client(timeout=900).ask(function,[q['service_args']for q in selected])
    (output/(function+'.json')).write_text(json.dumps(dict(questions=selected,reply=reply),indent=2)+'\n')
    # Worker inputs return typed values without literal-mode annotations.
    expected=[[{k:v for k,v in a.items()if k!='mode'}for a in q['service_args']]for q in selected]
    assert reply['inputs']==expected
    print(function,len(selected),'observed',reply['ok'],flush=True)
binding={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest()for p in [Path(__file__),source,Path(__file__).with_name('live_variance.py')]}
(output/'binding.json').write_text(json.dumps(binding,indent=2)+'\n')
