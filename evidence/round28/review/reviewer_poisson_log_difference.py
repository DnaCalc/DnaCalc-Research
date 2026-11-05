#!/usr/bin/env python3
"""Manual RN64 subtraction-cell inversion and original LN comparison audit."""
from collections import Counter
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys

root, library, run, output = map(Path, sys.argv[1:])
sys.path.insert(0, str(root / 'checker'))
from geomean import verified_study_answers
from exact import BINARY64 as B, EXTENDED64 as E, binary64
from x87 import Exact, Instructions, encode, decode, load64, published_bits

evidence = root / 'evidence/round28'
study_path, prior_path = [evidence / name for name in ('poisson-log-difference.json', 'poisson-exp-arguments.json')]
doc, prior = [json.loads(p.read_text()) for p in (study_path, prior_path)]
proof_path = run / 'reviewer-poisson-exp-arguments.json'
proof = json.loads(proof_path.read_text())
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
frozen = {str(p): digest(p) for p in (study_path, prior_path, proof_path, Path(__file__), library,
                   root / 'checker/geomean.py', root / 'checker/exact.py', root / 'checker/x87.py')}
assert proof['complete'] and proof['hashes'][str(prior_path)] == digest(prior_path)
assert doc['complete'] and prior['complete'] and doc['original_sha256'] == digest(prior_path)
assert doc['original_receipts'] == prior['original_receipts']
public, bindings, ordered = verified_study_answers(doc, evidence / 'poisson-log-service')
sources = [r for r in prior['records'] if r['stage'] == 'zero_to_one']
assert len(doc['records']) == len(sources) == len(public) == len(ordered) == 1024
exact = Exact(Instructions(library))
zero, one = encode(Q(0)), encode(Q(1))


def op(kind, a=zero, b=zero):
    return exact.operation(kind, a, b, 0x133f, False)


ln2, log2e = op('fldln2'), op('fldl2e')


def negative_exp(argument):
    assert argument < 0
    t = op('mul', encode(-argument), log2e)
    integer = op('frndint', t)
    small = op('f2xm1', op('sub', t, integer))
    positive = op('fscale', op('add', one, small), integer)
    return f'{published_bits(op("div", one, positive)):016x}'


def lattice_count(fmt, lower, upper, lower_closed, upper_closed):
    a, b = fmt.ceiling_coordinate(lower), fmt.floor_coordinate(upper)
    if not lower_closed and fmt.value(a) == lower:
        a += 1
    if not upper_closed and fmt.value(b) == upper:
        b -= 1
    return max(0, b - a + 1)


counts, used = Counter(), set()
direct_countermodels = Counter()
for row, source in zip(doc['records'], sources):
    assert (row['index'], row['input'], row['excel'], row['cell']) == (source['index'], source['input'], source['excel'], source['cell'])
    x, cell = row['input'], row['cell']
    mean = binary64(x)
    magnitude_low, magnitude_high = [decode(bytes.fromhex(cell[k])) for k in ('low', 'high')]
    low_t, high_t = -magnitude_high, -magnitude_low
    low_coordinate, high_coordinate = E.floor_coordinate(low_t), E.floor_coordinate(high_t)
    assert E.value(low_coordinate) == low_t and E.value(high_coordinate) == high_t
    lower = (E.value(low_coordinate - 1) + low_t) / 2 + mean
    upper = (high_t + E.value(high_coordinate + 1)) / 2 + mean
    lower_closed, upper_closed = low_coordinate % 2 == 0, high_coordinate % 2 == 0
    expected = dict(lower=[str(lower.numerator), str(lower.denominator)],
                    upper=[str(upper.numerator), str(upper.denominator)],
                    lower_closed=lower_closed, upper_closed=upper_closed)
    assert row['private_log_interval'] == expected
    def contains(value):
        return (value > lower or value == lower and lower_closed) and (value < upper or value == upper and upper_closed)
    key = 'LN', json.dumps([dict(type='number', bits=x)], sort_keys=True)
    used.add(key)
    answer = public[key]
    assert row['public_LN'] == answer and answer['type'] == 'number'
    resident = op('fyl2x', ln2, load64(int(x, 16)))
    assert row['resident_LN'] == resident.hex()
    public_value, resident_value = binary64(answer['bits']), decode(resident)
    public_possible, resident_possible = contains(public_value), contains(resident_value)
    assert (row['public_LN_possible'], row['resident_LN_possible']) == (public_possible, resident_possible)
    pc64 = lattice_count(E, lower, upper, lower_closed, upper_closed)
    binary64_count = lattice_count(B, lower, upper, lower_closed, upper_closed)
    assert (row['pc64_log_count'], row['binary64_log_count']) == (pc64, binary64_count)
    counts['public_LN_possible' if public_possible else 'public_LN_excluded'] += 1
    counts['resident_LN_possible' if resident_possible else 'resident_LN_excluded'] += 1
    counts['binary64_log_possible' if binary64_count else 'binary64_log_impossible'] += 1
    for name, value, possible in (('public_LN', public_value, public_possible), ('resident_LN', resident_value, resident_possible)):
        argument = E.round(value - mean)
        assert possible == (low_t <= argument <= high_t)
        predicted = negative_exp(argument)
        assert (predicted == row['excel']['bits']) == possible
        direct_countermodels[name + ('_matched' if possible else '_contradicted')] += 1
assert used == set(public) and counts == doc['counts']
assert frozen == {p: digest(Path(p)) for p in frozen}
assert all(digest(Path(p)) == h for p, h in bindings.items())
receipt = dict(complete=True, original_PMF1_targets=1024, public_LN_calls=len(ordered),
    exact_subtraction_inverse_domains=1024, counts=counts, direct_native_countermodel_outputs=direct_countermodels,
    source=str(Path(__file__).resolve()), hashes=frozen | bindings,
    limits='Uses independently checked original EXP preimages and manually constructed signed PC64 midpoint ownership, then exact translation by the observed mean. Registered public LN and native resident FYL2X are comparisons only.166 domains exclude every binary64 log for this single RN64 subtraction/declared EXP model; no private logarithm or full PMF implementation is identified.')
with output.open('x') as stream:
    json.dump(receipt, stream, indent=2); stream.write('\n')
print(json.dumps({k: v for k, v in receipt.items() if k not in ('hashes', 'limits', 'source')}))
