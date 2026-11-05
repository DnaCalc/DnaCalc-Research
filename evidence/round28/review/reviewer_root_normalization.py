#!/usr/bin/env python3
"""Read-only review: exact exponent reconstruction and full call/loss audit.

Usage: python3 reviewer_root_normalization.py RESEARCH_ROOT RUN_DIRECTORY OUTPUT
No product/search implementation is imported. The maintained independent raw
receipt verifier validates all Excel inputs, outputs and worker provenance.
Counterfactual arithmetic below only rewrites positive normal exponent bits.
"""
import hashlib
import json
from collections import Counter
from pathlib import Path
import sys

research, run, output = map(Path, sys.argv[1:])
sys.path.insert(0, str(research / 'checker'))
from geomean import verified_study_answers

evidence = research / 'evidence/round28'
sources = [evidence / 'root-scales.json', evidence / 'root-normalization.json']
sources += [run / name for name in (
    'root-original-inputs.json', 'root-normalized-inputs.json',
    'root-original-before.jsonl', 'root-normalized-before.jsonl',
    'root-original-after.jsonl', 'root-normalized-after.jsonl',
    'root-counterfactual-delta.json', 'root-actual-delta.json')]
sources += [Path(__file__), research / 'checker/geomean.py']
digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
hashes = {str(path): digest(path) for path in sources}
read = lambda path: json.loads(path.read_text())
lines = lambda path: [json.loads(line) for line in path.read_text().splitlines()]
originals = read(evidence / 'root-scales.json')
study = read(evidence / 'root-normalization.json')
_, original_bindings, original_order = verified_study_answers(
    originals, evidence / 'root-scales-service')
_, control_bindings, control_order = verified_study_answers(
    study, evidence / 'root-normalization-service')
before = lines(run / 'root-original-before.jsonl')
after = lines(run / 'root-original-after.jsonl')
core = lines(run / 'root-normalized-before.jsonl')
core_after = lines(run / 'root-normalized-after.jsonl')
original_input = read(run / 'root-original-inputs.json')
core_input = read(run / 'root-normalized-inputs.json')
assert original_input['function'] == core_input['function'] == 'GEOMEAN'
assert len(original_order) == len(before) == len(after) == len(original_input['witnesses']) == 19200
assert len(core) == len(core_input['witnesses']) == 6016
assert len(control_order) == 7680 and core_after == core
one = '0x3ff0000000000000'
core_map = {}
for row, witness in zip(core, core_input['witnesses']):
    assert row['id'] == witness['id'] and row['args'] == witness['args']
    assert row['not_excel_evidence'] is True
    n = len(row['args'])
    assert all(x == one for x in row['args'][1:])
    key = n, row['args'][0][2:]
    assert key not in core_map
    core_map[key] = row['calc_func_prediction'][2:]


def shift(bits, exponent):
    value = int(bits, 16)
    encoded = ((value >> 52) & 2047) + exponent
    assert value >> 63 == 0 and 0 < encoded < 2047
    return f'{(value & ((1 << 52) - 1)) | (encoded << 52):016x}'


counts, changes, losses, used = Counter(), [], [], set()
for index, ((function, args, answer), old, new, witness, record) in enumerate(zip(
        original_order, before, after, original_input['witnesses'], study['records'])):
    assert function == 'GEOMEAN' and answer['type'] == 'number'
    assert old['not_excel_evidence'] is True and new['not_excel_evidence'] is True
    assert old['args'] == new['args'] == witness['args'] == ['0x' + arg['bits'] for arg in args]
    assert old['id'] == new['id'] == witness['id'] == f'r28-root-original-{index}'
    assert all(x == one for x in old['args'][1:])
    n, bits = len(args), args[0]['bits']
    exponent = ((int(bits, 16) >> 52) & 2047) - 1023
    quotient = (abs(exponent) // n) * (-1 if exponent < 0 else 1)
    reduced = shift(bits, -quotient * n)
    used.add((n, reduced))
    model = next(c for c in record['conventions']
                 if c['convention'] == 'adjusted_frexp_truncation')
    assert model['reduced'] == reduced and model['quotient'] == quotient
    expected = '0x' + answer['bits']
    got = old['calc_func_prediction']
    calculated_core = core_map[n, reduced]
    predicted = '0x' + shift(calculated_core, quotient)
    assert predicted == new['calc_func_prediction']
    assert shift(model['measured_normalized_root']['bits'], quotient) == answer['bits']
    counts['total'] += 1
    counts['before_pass'] += got == expected
    counts['after_pass'] += predicted == expected
    if got != predicted:
        counts['changed_outputs'] += 1
        counts['gains'] += got != expected and predicted == expected
        counts['losses'] += got == expected and predicted != expected
        row = dict(index=index, count=n, input=bits, expected=expected, before=got,
                   after=predicted, reduced=reduced, quotient=quotient,
                   measured_reduced_root=model['measured_normalized_root'],
                   calculated_reduced_root='0x' + calculated_core)
        changes.append(row)
        if got == expected and predicted != expected:
            assert calculated_core != model['measured_normalized_root']['bits']
            losses.append(dict(row, classification='current private core differs from independently measured root at the exact reduced product'))
assert used == set(core_map)
delta = read(run / 'root-counterfactual-delta.json')
actual = read(run / 'root-actual-delta.json')
assert counts == delta['counts'] == actual['counts']
assert changes == delta['changed'] and losses == delta['losses']
assert actual['all_changed_predictions_equal_counterfactual'] and actual['normalized_core_unchanged']
for path, expected in actual['input_sha256'].items():
    assert digest(Path(path)) == expected
assert hashes == {str(path): digest(path) for path in sources}
receipt = dict(complete=True,
    scope='Complete original/control provenance, integer exponent reconstruction, every actual output and every classified loss. No private-kernel identification.',
    counts=counts, classified_losses=losses, verified_original_calls=len(original_order),
    verified_normalized_calls=len(control_order), unique_adjusted_core_app_calls=len(core_map),
    source=str(Path(__file__).resolve()), hashes=hashes | original_bindings | control_bindings)
with output.open('x') as stream:
    json.dump(receipt, stream, indent=2)
    stream.write('\n')
print(json.dumps(dict(counts=counts, classified_losses=len(losses), receipt=str(output))))
