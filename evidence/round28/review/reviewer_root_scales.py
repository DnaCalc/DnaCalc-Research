#!/usr/bin/env python3
"""Validate the entire frozen scale question grid and relative records."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

root, output = map(Path, sys.argv[1:])
sys.path.insert(0, str(root / 'checker'))
from geomean import verified_study_answers

path = root / 'evidence/round28/root-scales.json'
doc = json.loads(path.read_text())
_, bindings, ordered = verified_study_answers(doc, path.parent / 'root-scales-service')
assert doc['consecutive_inputs_per_center'] == 128
centers = (0x3ff4cccccccccccd, 0x401c000000000000, 0x408f400000000000)
xs = [f'{b:016x}' for center in centers for b in range(center - 64, center + 64)]
arities = (2, 3, 4, 5, 7, 9, 10, 11, 16, 32)
scales = (0, -20, -1, 1, 20)
one = dict(type='number', bits='3ff0000000000000')


def shift(bits, exponent):
    x = int(bits, 16)
    e = ((x >> 52) & 2047) + exponent
    assert x >> 63 == 0 and 0 < e < 2047
    return f'{(x & ((1 << 52) - 1)) | (e << 52):016x}'


expected = [('GEOMEAN', [dict(type='number', bits=shift(x, n * s))] + [one] * (n - 1))
            for n in arities for s in scales for x in xs]
assert len(ordered) == len(expected) == 19200
assert [(name, args) for name, args, _ in ordered] == expected
answers = {(len(args), args[0]['bits']): answer for _, args, answer in ordered}
records = []
counts = {str(n): Counter() for n in arities}
by_scale = {str(s): Counter() for s in scales[1:]}
for n in arities:
    for s in scales[1:]:
        for x in xs:
            scaled_input = shift(x, n * s)
            base, scaled = answers[n, x], answers[n, scaled_input]
            assert base['type'] == scaled['type'] == 'number'
            predicted = shift(base['bits'], s)
            name = 'equal' if scaled['bits'] == predicted else 'different'
            counts[str(n)][name] += 1
            by_scale[str(s)][name] += 1
            records.append(dict(count=n, output_scale=s, input=x,
                scaled_input=scaled_input, base=base, scaled=scaled,
                exact_rescaling_of_observed_base=predicted))
assert records == doc['records'] and counts == doc['counts']
assert by_scale['1'] == by_scale['20'] == Counter(equal=3840)
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
result = dict(complete=True, original_calls=len(ordered), relative_records=len(records),
    counts_by_scale=by_scale, source=str(Path(__file__).resolve()),
    hashes={str(p): digest(p) for p in (path, Path(__file__), root / 'checker/geomean.py')} | bindings,
    limits='Observed normal unit-factor grid only. Exact positive transfer is not a universal algorithm identification; negative disagreements motivate, rather than fit, the separate normalization control study.')
with output.open('x') as stream:
    json.dump(result, stream, indent=2)
    stream.write('\n')
print(json.dumps(dict(original_calls=len(ordered), relative_records=len(records), counts_by_scale=by_scale)))
