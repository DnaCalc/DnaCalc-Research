#!/usr/bin/env python3
"""Replay exactly the inferred stored-LN/stored-reciprocal/stored-product model.

No model enumeration or product implementation is imported. Native leaves are
separate single instructions; other operations use independent exact arithmetic.
Until paired controls are read, these are conditional native-model predictions.
"""
from collections import Counter
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys

root, library, output = map(Path, sys.argv[1:])
sys.path.insert(0, str(root / 'checker'))
from exact import BINARY64
from geomean import verified_study_answers
from x87 import Exact, Instructions, encode, decode, load64, published_bits

evidence = root / 'evidence/round28'
path = evidence / 'root-normalization.json'
study = json.loads(path.read_text())
answers, bindings, _ = verified_study_answers(study, evidence / 'root-normalization-service')
expected = {}
for record in study['records']:
    item = next(c for c in record['conventions'] if c['convention'] == 'adjusted_frexp_truncation')
    expected[record['count'], item['reduced']] = item['measured_normalized_root']
assert len(expected) == 6016
exact = Exact(Instructions(library))
zero, one = encode(Q(0)), encode(Q(1))


def operation(kind, a=zero, b=zero, stored=False):
    return exact.operation(kind, a, b, 0x133f, stored)


ln2, log2e = operation('fldln2'), operation('fldl2e')
counts, records = Counter(), []
for (n, bits), answer in sorted(expected.items()):
    args = [dict(type='number', bits=bits)] + [dict(type='number', bits='3ff0000000000000')] * (n - 1)
    assert answers['GEOMEAN', json.dumps(args, sort_keys=True)] == answer
    log = operation('fyl2x', ln2, load64(int(bits, 16)), stored=True)
    reciprocal = encode(BINARY64.round(Q(1, n)))
    argument = operation('mul', log, reciprocal, stored=True)
    t = operation('mul', argument, log2e)
    negative = decode(t) < 0
    if negative:
        t = operation('neg', t)
    integer = operation('frndint', t)
    fraction = operation('sub', t, integer)
    small = operation('f2xm1', fraction)
    result = operation('fscale', operation('add', one, small), integer)
    if negative:
        result = operation('div', one, result)
    predicted = dict(type='number', bits=f'{published_bits(result):016x}')
    counts['matched' if predicted == answer else 'mismatched'] += 1
    records.append(dict(count=n, input=bits, log=f'{published_bits(log):016x}',
                        reciprocal=f'{published_bits(reciprocal):016x}',
                        argument=f'{published_bits(argument):016x}',
                        predicted=predicted, excel=answer))
assert counts == Counter(matched=6016)
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
receipt = dict(complete=True, model='EXP(stored_PC64(LN(M)*RN53(1/n))), with LN(M) published tobinary64 first',
    counts=counts, records=records, source=str(Path(__file__).resolve()),
    hashes={str(p): digest(p) for p in (path, Path(__file__), library, root / 'checker/x87.py', root / 'checker/exact.py')} | bindings,
    limits='Fixed inferred stage model on6016 reduced controls. Public LN/EXP instruction chains and their native CPU leaves remain conditional until independent paired public-leaf Excel readings verify them. No private polynomial/table or universal code identification follows.')
with output.open('x') as stream:
    json.dump(receipt, stream, indent=2)
    stream.write('\n')
print(json.dumps(dict(counts=counts, receipt=str(output))))
