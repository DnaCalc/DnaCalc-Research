#!/usr/bin/env python3
"""Raw PMF identities, independent EXP inverse boundaries and B64 projection."""
from collections import Counter
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys

root, library, output = map(Path, sys.argv[1:])
sys.path.insert(0, str(root / 'checker'))
from geomean import verified_study_answers
from exact import BINARY64 as B, EXTENDED64 as E, binary64
from x87 import Exact, Instructions, encode, decode, load64, published_bits

evidence = root / 'evidence/round28'
study_path, original_path = [evidence / p for p in ('poisson-exp-arguments.json', 'poisson-seeds.json')]
doc, original = [json.loads(p.read_text()) for p in (study_path, original_path)]
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
sources = (study_path, original_path, library, Path(__file__), root / 'checker/geomean.py',
           root / 'checker/exact.py', root / 'checker/x87.py', root / 'racer/fingerprint.py')
frozen = {str(p): digest(p) for p in sources}
assert doc['complete'] and original['complete']
assert doc['original_sha256'] == digest(original_path) and doc['original_receipts'] == original['receipts']
raw, original_bindings, original_order = verified_study_answers(original, evidence / 'poisson-service')
public, public_bindings, public_order = verified_study_answers(doc, evidence / 'poisson-exp-service')
nodes = [dict(op='input', a=0), dict(op='const', bits='3ff0000000000000'),
    dict(op='fldl2e', a=0, b=0), dict(op='mul', a=0, b=2),
    dict(op='frndint', a=3, b=0), dict(op='sub', a=3, b=4),
    dict(op='f2xm1', a=5, b=0), dict(op='add', a=1, b=6),
    dict(op='fscale', a=7, b=4), dict(op='div', a=1, b=8)]
assert doc['graph'] == dict(nodes=nodes, outputs=[9, 9], sites=0)
exact = Exact(Instructions(library))
zero, one = encode(Q(0)), encode(Q(1))


def op(kind, a=zero, b=zero):
    return exact.operation(kind, a, b, 0x133f, False)


log2e = op('fldl2e')


def exponential(magnitude):
    t = op('mul', magnitude, log2e)
    integer = op('frndint', t)
    small = op('f2xm1', op('sub', t, integer))
    positive = op('fscale', op('add', one, small), integer)
    return f'{published_bits(op("div", one, positive)):016x}'


def number(bits):
    return dict(type='number', bits=bits)


def observation(name, args):
    return raw[name, json.dumps(args, sort_keys=True)]


assert len(doc['records']) == len(original['records']) == 2048
assert len({(r['input'], r['stage']) for r in doc['records']}) == 2048
xs = {r['input'] for r in original['records']}
assert len(xs) == 1024
assert {(r['input'], r['stage']) for r in doc['records']} == {(x, s) for x in xs for s in ('zero_to_one', 'one_to_two')}
counts, projection_sizes, used = Counter(), Counter(), set()
boundaries, matched_neighbors = 0, 0
for index, (row, source) in enumerate(zip(doc['records'], original['records'])):
    x, stage = row['input'], row['stage']
    assert (row['index'], x, stage, row['excel']) == (index, source['input'], source['stage'], source['published_output'])
    k = 1 if stage == 'zero_to_one' else 2
    args = [number(f'{published_bits(encode(Q(k))):016x}'), number(x), dict(type='boolean', value=False, mode='cell')]
    assert row['excel'] == observation('POISSON.DIST', args) == observation('POISSON', args)
    assert row['excel']['type'] == 'number' and binary64(row['excel']['bits']) > 0
    cell = row['cell']
    lo, hi = [bytes.fromhex(cell[z]) for z in ('low', 'high')]
    low_value, high_value = decode(lo), decode(hi)
    assert 0 < low_value <= high_value
    low_coordinate, high_coordinate = E.floor_coordinate(low_value), E.floor_coordinate(high_value)
    assert E.value(low_coordinate) == low_value and E.value(high_coordinate) == high_value
    arguments = [encode(E.value(c)) for c in (low_coordinate - 1, low_coordinate, high_coordinate, high_coordinate + 1)]
    got = [exponential(a) for a in arguments]
    boundaries += 4
    target = row['excel']['bits']
    assert got[1] == got[2] == target and int(got[0], 16) > int(target, 16) > int(got[3], 16)
    assert [cell['predecessor_output'], cell['successor_output']] == [got[0], got[3]]
    a, b = B.ceiling_coordinate(low_value), B.floor_coordinate(high_value)
    count = max(0, b - a + 1)
    projection_sizes[count] += 1
    first = [f'{published_bits(encode(B.value(c))):016x}' for c in range(a, min(b + 1, a + 8))]
    assert (row['binary64_argument_count'], row['first_binary64_magnitudes']) == (count, first)
    status = 'binary64 argument possible' if count else 'no binary64 argument'
    assert row['status'] == status
    counts[status] += 1
    if not count:
        predecessor, successor = B.floor_coordinate(low_value), B.ceiling_coordinate(high_value)
        assert successor - predecessor == 1
        neighbors = [f'{published_bits(encode(B.value(c))):016x}' for c in (predecessor, successor)]
        assert row['magnitude_neighbors'] == neighbors
        measured = []
        for bits in neighbors:
            key = 'EXP', json.dumps([number(f'{int(bits, 16) | (1 << 63):016x}')], sort_keys=True)
            used.add(key)
            answer = public[key]
            assert answer == number(exponential(load64(int(bits, 16))))
            matched_neighbors += 1
            measured.append(answer)
        assert row['public_EXP_neighbors'] == measured
        assert row['public_bracket_verified'] is True
        assert int(measured[1]['bits'], 16) < int(target, 16) < int(measured[0]['bits'], 16)
        counts['public bracket verified'] += 1
assert used == set(public) and counts == doc['counts']
assert frozen == {p: digest(Path(p)) for p in frozen}
bindings = original_bindings | public_bindings
assert all(digest(Path(p)) == h for p, h in bindings.items())
receipt = dict(complete=True, original_source_calls=len(original_order), original_PMF_targets=2048,
    inverse_boundary_checks=boundaries, counts=counts, projection_sizes=projection_sizes,
    public_calls=len(public_order), public_unique_calls=len(public), repeated_neighbor_checks=matched_neighbors,
    source=str(Path(__file__).resolve()), hashes=frozen | bindings,
    limits='Original PMF1/2 and alias identities, all closed PC64 inverse endpoints and immediate lattice neighbors, exact binary64 projection, complete public EXP gap controls. Conditional on the declared native EXP primitive and monotone-preimage premise. This excludes a stored argument within that model for682 targets; it identifies no resident argument, private exponential, or full Poisson algorithm.')
with output.open('x') as stream:
    json.dump(receipt, stream, indent=2)
    stream.write('\n')
print(json.dumps({k: v for k, v in receipt.items() if k not in ('hashes', 'limits', 'source')}))
