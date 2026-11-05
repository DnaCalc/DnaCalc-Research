#!/usr/bin/env python3
"""Independent raw-reading and rational midpoint audit of all root edges."""
from collections import Counter
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import struct
import sys

root, output = map(Path, sys.argv[1:])
sys.path.insert(0, str(root / 'checker'))
from geomean import verified_study_answers
from exact import binary64, BINARY64

path = root / 'evidence/round28/root-output-edges-corrected.json'
doc = json.loads(path.read_text())
answers, bindings, ordered = verified_study_answers(doc, path.parent / 'root-output-edge-service')
one = dict(type='number', bits='3ff0000000000000')
used = set()
targets = [f'{0x3ff0000000000000 + k * (1 << 49):016x}' for k in range(9)]
arities = (2, 3, 4, 5, 7, 9, 10, 11, 16, 32)
assert [(r['count'], r['target']) for r in doc['records']] == [(n, t) for n in arities for t in targets]


def observed(n, bits):
    args = [dict(type='number', bits=bits)] + [one] * (n - 1)
    key = 'GEOMEAN', json.dumps(args, sort_keys=True)
    used.add(key)
    result = answers[key]
    assert result['type'] == 'number'
    return result


counts = Counter()
for record in doc['records']:
    n, target = record['count'], record['target']
    target_bits = int(target, 16)
    midpoint = (binary64(target) + binary64(f'{target_bits - 1:016x}')) / 2
    ideal = midpoint ** n
    assert Q(record['midpoint']) == midpoint and Q(record['ideal_input']) == ideal
    center = int.from_bytes(struct.pack('>d', float(BINARY64.round(ideal))), 'big')
    low, high = center - 1024, center + 1024
    trace = record['trace']
    assert [r['input'] for r in trace[:2]] == [f'{low:016x}', f'{high:016x}']
    assert not trace[0]['side'] and trace[1]['side']
    for index, row in enumerate(trace):
        assert row['answer'] == observed(n, row['input'])
        assert row['side'] is (int(row['answer']['bits'], 16) >= target_bits)
        if index >= 2:
            point = (low + high) // 2
            assert row['input'] == f'{point:016x}'
            if row['side']:
                high = point
            else:
                low = point
    assert high - low == 1 and record['low'] == f'{low:016x}' and record['high'] == f'{high:016x}'
    window = record['window']
    xs = [f'{b:016x}' for b in range(high - 64, high + 65)]
    ys = [observed(n, x) for x in xs]
    assert window['inputs'] == xs and window['answers'] == ys
    transitions, decreases = [], []
    for index in range(1, len(xs)):
        a, b = ys[index - 1:index + 1]
        row = dict(low=xs[index - 1], high=xs[index], low_answer=a, high_answer=b)
        if (int(a['bits'], 16) >= target_bits) != (int(b['bits'], 16) >= target_bits):
            transitions.append(row)
        if int(a['bits'], 16) > int(b['bits'], 16):
            decreases.append(row)
    assert window['transitions'] == transitions and window['decreases'] == decreases
    counts['unique_local_crossing' if len(transitions) == 1 else 'multiple_local_crossings'] += 1
    counts['monotone_local_window' if not decreases else 'nonmonotone_local_window'] += 1
    exponent = ((high >> 52) & 2047) - 1023 - 52
    unit = Q(2 ** exponent) if exponent >= 0 else Q(1, 2 ** -exponent)
    assert Q(record['reference_input_ulp']) == unit
    offsets = [(binary64(f'{b:016x}') - ideal) / unit for b in (low, high)]
    assert list(map(Q, record['input_offset_from_exact_midpoint_power_in_input_ulps'])) == offsets
assert used == set(answers) and counts == doc['counts']
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
receipt = dict(complete=True, original_calls=len(ordered), unique_original_calls=len(answers),
    edge_count=len(doc['records']), counts=counts, source=str(Path(__file__).resolve()),
    hashes={str(p): digest(p) for p in (path, Path(__file__), root / 'checker/geomean.py', root / 'checker/exact.py')} | bindings,
    limits='Every original answer, bisection update and consecutive window verified. Rational midpoint powers locate questions only. Local observations do not prove global monotonicity or identify private numerical primitives.')
with output.open('x') as stream:
    json.dump(receipt, stream, indent=2)
    stream.write('\n')
print(json.dumps(dict(original_calls=len(ordered), edge_count=len(doc['records']), counts=counts)))
