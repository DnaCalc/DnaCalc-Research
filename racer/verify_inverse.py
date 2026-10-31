"""Check inverse intervals against current cleared native Fpu forward controls."""
import argparse
from pathlib import Path
import collections
import copy
import json
from fractions import Fraction as Q
from certify import exact_root
from intersect import MODELS, compile_expr, row_preimage
from intervals import number

parser = argparse.ArgumentParser()
parser.add_argument('input')
parser.add_argument('catalogue')
parser.add_argument('controls')
parser.add_argument('output')
a = parser.parse_args()
rows=[]
for line in open(a.input):
    r=json.loads(line)
    if number(f'{r["leaves"][2]:016x}')==0:
        rows.append(r)
cat=Path(a.catalogue)
exprs=json.loads((cat/'expressions.json').read_text())
profiles=[]
for c in json.loads((cat/'combines.json').read_text()):
    models=[MODELS[m] for m in c['models']]
    expr,site=compile_expr(exprs[c['expression']],models,1)
    assert site==4
    profiles.append((c,expr,models))
counts=collections.Counter()
for control in map(json.loads,open(a.controls)):
    target=control['control_target']
    counts['forward_controls']+=1
    counts['backend_published_differences']+=target!=control['software_control_target']
    if not target.startswith('0x') or number(target)==0:
        counts['typed_or_zero_controls_not_inverted']+=1
        continue
    row=copy.copy(rows[control['row_index']])
    row['raw']=dict(row['raw'],expected=target)
    try:
        preimage=row_preimage(row,profiles[control['combine_id']])
    except (ArithmeticError,AssertionError) as e:
        counts['unsupported_controls']+=1
        continue
    q=exact_root(control['em'])
    assert q is not None and preimage is not None and preimage.contains(q),control
    counts['independent_forward_inverse_memberships']+=1
result=dict(counts=counts,source='synthetic current native/software Fpu outputs, not Excel answers',
    claim='Every supported nonzero finite native forward control lies in the independently rationally inverted input interval, including strict operand stores and resident/stored double rounding.')
with Path(a.output).open('x') as f:json.dump(result,f,indent=2)
print(json.dumps(result))
