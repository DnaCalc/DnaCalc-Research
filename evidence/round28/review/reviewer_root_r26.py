#!/usr/bin/env python3
"""Fixed stored-root stage replay on every original R26 GEOMEAN call."""
from collections import Counter
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys

root, library, baseline, output = map(Path, sys.argv[1:])
sys.path.insert(0, str(root / 'checker'))
from exact import BINARY64, EXTENDED64
from geomean import verified_study_answers
from x87 import Exact, Instructions, encode, decode, load64, published_bits

study_path = root / 'evidence/round26/roots.json'
study = json.loads(study_path.read_text())
_, bindings, ordered = verified_study_answers(study, root / 'evidence/round26/fingerprints')
assert len(ordered) == 192512
ordered = [r for r in ordered if r[0] == 'GEOMEAN']
assert len(ordered) == 31744
actual_path = baseline / 'GEOMEAN.actuals.jsonl'
input_path = baseline / 'GEOMEAN.inputs.json'
receipt_path = baseline / 'receipt.json'
actuals = [json.loads(s) for s in actual_path.read_text().splitlines()]
inputs = json.loads(input_path.read_text())
assert inputs['function'] == 'GEOMEAN'
assert len(actuals) == len(inputs['witnesses']) == len(ordered)
exact = Exact(Instructions(library))
zero, one = encode(Q(0)), encode(Q(1))


def operation(kind, a=zero, b=zero, stored=False):
    return exact.operation(kind, a, b, 0x133f, stored)


def shift(bits, exponent):
    b = int(bits, 16)
    e = ((b >> 52) & 2047) + exponent
    assert b >> 63 == 0 and 0 < e < 2047
    return f'{(b & ((1 << 52) - 1)) | (e << 52):016x}'


ln2, log2e = operation('fldln2'), operation('fldl2e')
counts, arities, changes, losses = Counter(), Counter(), [], []
for index, ((function, args, excel), before, witness) in enumerate(zip(ordered, actuals, inputs['witnesses'])):
    assert function == 'GEOMEAN' and excel['type'] == 'number'
    assert all(arg['type'] == 'number' for arg in args)
    assert all(arg['bits'] == '3ff0000000000000' for arg in args[1:])
    assert before['not_excel_evidence'] is True
    assert before['id'] == witness['id'] == f'original-service-{index}'
    assert before['args'] == witness['args'] == ['0x' + arg['bits'] for arg in args]
    bits, n = args[0]['bits'], len(args)
    exponent = ((int(bits, 16) >> 52) & 2047) - 1023
    q = (abs(exponent) // n) * (-1 if exponent < 0 else 1)
    reduced = shift(bits, -q * n)
    if n == 1:
        prediction, log_bits, argument_bits = bits, None, None
    else:
        log = operation('fyl2x', ln2, load64(int(reduced, 16)), stored=True)
        reciprocal = encode(BINARY64.round(EXTENDED64.round(Q(1, n))))
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
        prediction = shift(f'{published_bits(result):016x}', q)
        log_bits = f'{published_bits(log):016x}'
        argument_bits = f'{published_bits(argument):016x}'
    expected, old, new = '0x' + excel['bits'], before['calc_func_prediction'], '0x' + prediction
    was, now = old == expected, new == expected
    counts['total'] += 1
    counts['before_pass'] += was
    counts['after_pass'] += now
    counts['gains'] += now and not was
    counts['losses'] += was and not now
    counts['changed_outputs'] += old != new
    arities[n] += 1
    if old != new or not now:
        row = dict(index=index, count=n, input=bits, reduced=reduced, quotient=q,
                   stored_log=log_bits, stored_argument=argument_bits,
                   expected=expected, before=old, after=new,
                   classification='gained_exact_match' if now and not was else
                                  'lost_exact_match' if was and not now else 'changed_still_failing')
        changes.append(row)
        if was and not now:
            losses.append(row)
baseline_receipt = json.loads(receipt_path.read_text())
assert baseline_receipt['counts']['GEOMEAN']['matched'] == counts['before_pass'] == 29310
assert counts['total'] - counts['before_pass'] == 2434
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
for name, expected in baseline_receipt['output_sha256'].items():
    assert digest(baseline / name) == expected
receipt = dict(complete=True, counts=counts, arities=arities,
    changed_or_still_failing=changes, losses=losses,
    source=str(Path(__file__).resolve()),
    model='Exact q normalization; original-count singleton identity; LN tobinary64; reciprocal PC64->binary64; log*reciprocal PC64->binary64; existing public nativeEXP chain; exact exponent restoration.',
    hashes={str(p): digest(p) for p in (study_path, actual_path, input_path, receipt_path, Path(__file__), library,
             root / 'checker/geomean.py', root / 'checker/exact.py', root / 'checker/x87.py')} | bindings,
    limits='Complete independent replay of original R26 GEOMEAN observations. This is a fixed native-model counterfactual, not new Excel evidence or actual product execution. Private implementation internals and other routes remain open.')
with output.open('x') as stream:
    json.dump(receipt, stream, indent=2)
    stream.write('\n')
print(json.dumps(dict(counts=counts, arities=arities, losses=len(losses), receipt=str(output))))
