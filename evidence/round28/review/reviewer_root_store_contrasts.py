#!/usr/bin/env python3
"""Independent four-model single-publication contrast audit; no enumeration."""
from collections import Counter
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys

root, library, output = map(Path, sys.argv[1:])
sys.path.insert(0, str(root / 'checker'))
from x87 import Exact, Instructions, encode, decode, load64, published_bits

evidence = root / 'evidence/round28'
study_path = evidence / 'root-store-contrasts.json'
prior_path = evidence / 'root-public-stages.json'
doc, prior = [json.loads(p.read_text()) for p in (study_path, prior_path)]
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
frozen = {str(p): digest(p) for p in (study_path, prior_path, library, Path(__file__),
                  root / 'checker/x87.py', root / 'checker/exact.py')}
assert doc['complete'] and prior['complete'] and doc['original_sha256'] == digest(prior_path)
masks = dict(stored=7, resident_log=6, resident_reciprocal=5, resident_argument=3)
assert list(doc['masks'].items()) == list(masks.items())


def grammar(n):
    nodes = [dict(op='input', a=0), dict(op='const', bits='3ff0000000000000'),
        dict(op='fldln2', a=0, b=0), dict(op='fyl2x', a=2, b=0, site=0),
        dict(op='const', bits=f'{published_bits(encode(Q(n))):016x}'),
        dict(op='div', a=1, b=4, site=1), dict(op='mul', a=3, b=5, site=2),
        dict(op='fldl2e', a=0, b=0), dict(op='mul', a=6, b=7),
        dict(op='frndint', a=8, b=0), dict(op='sub', a=8, b=9),
        dict(op='f2xm1', a=10, b=0), dict(op='add', a=1, b=11), dict(op='fscale', a=12, b=9)]
    return dict(nodes=nodes, outputs=[13, 13], sites=3)


expected_arities = sorted({r['count'] for r in prior['records']})
assert doc['graphs'] == {str(n): grammar(n) for n in expected_arities}
assert len(doc['records']) == len(prior['records']) == 6016
exact = Exact(Instructions(library))
zero, one = encode(Q(0)), encode(Q(1))


def operation(kind, a=zero, b=zero, stored=False):
    return exact.operation(kind, a, b, 0x133f, stored)


ln2, log2e = operation('fldln2'), operation('fldl2e')
counts = {name: Counter() for name in masks}
reciprocal_counts = {name: Counter() for name in masks}
branch_differences = {name: [] for name in masks}
for index, (row, original) in enumerate(zip(doc['records'], prior['records'])):
    n, x = row['count'], row['input']
    assert (n, x, row['excel']) == (original['count'], original['input'], original['excel'])
    assert list(row['predictions']) == list(masks)
    for name, mask in masks.items():
        log = operation('fyl2x', ln2, load64(int(x, 16)), stored=bool(mask & 1))
        reciprocal = operation('div', one, encode(Q(n)), stored=bool(mask & 2))
        argument = operation('mul', log, reciprocal, stored=bool(mask & 4))
        t = operation('mul', argument, log2e)
        whole = operation('frndint', t)
        small = operation('f2xm1', operation('sub', t, whole))
        result = operation('fscale', operation('add', one, small), whole)
        predicted = f'{published_bits(result):016x}'
        assert row['predictions'][name] == predicted
        counts[name]['matched' if row['excel']['bits'] == predicted else 'mismatched'] += 1
        alternative = predicted
        if decode(t) < 0:
            positive = operation('neg', t)
            integer = operation('frndint', positive)
            fraction = operation('sub', positive, integer)
            positive_root = operation('fscale', operation('add', one, operation('f2xm1', fraction)), integer)
            alternative = f'{published_bits(operation("div", one, positive_root)):016x}'
        reciprocal_counts[name]['matched' if row['excel']['bits'] == alternative else 'mismatched'] += 1
        if alternative != predicted:
            branch_differences[name].append(dict(index=index, count=n, input=x, plain=predicted,
                                                  reciprocal=alternative, excel=row['excel']))
assert counts == doc['counts']
assert frozen == {path: digest(Path(path)) for path in frozen}
receipt = dict(complete=True, verified_graphs=len(expected_arities), verified_outputs=6016 * 4,
    masks=masks, counts=counts, reciprocal_negative_branch_counts=reciprocal_counts,
    negative_branch_differences=branch_differences, source=str(Path(__file__).resolve()),
    hashes=frozen, limits='Exactly four fixed single-stage contrasts, all other graph nodes unchanged. The declared graph applies signed FRNDINT directly; reciprocal-negative EXP is an additional diagnostic of branch scope, not another fitted stage. No private algorithm identification follows.')
with output.open('x') as stream:
    json.dump(receipt, stream, indent=2)
    stream.write('\n')
print(json.dumps(dict(counts=counts, reciprocal_negative_branch_counts=reciprocal_counts,
                     negative_branch_differences={k: len(v) for k, v in branch_differences.items()})))
