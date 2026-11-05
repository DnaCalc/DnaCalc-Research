#!/usr/bin/env python3
"""Check every actual new core call against original raw Excel observations."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

research, run, output = map(Path, sys.argv[1:])
sys.path.insert(0, str(research / 'checker'))
from geomean import verified_study_answers
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
frozen = {}


def read(path):
    frozen[str(path)] = digest(path)
    return json.loads(path.read_text())


def lines(path):
    frozen[str(path)] = digest(path)
    return [json.loads(s) for s in path.read_text().splitlines()]


def verify_receipt(directory):
    doc = read(directory / 'receipt.json')
    assert doc['complete']
    for key in ('inputs_sha256', 'input_sha256'):
        for path, expected in doc.get(key, {}).items():
            assert digest(Path(path)) == expected
            frozen[path] = expected
    for name, expected in doc.get('output_sha256', {}).items():
        assert digest(directory / name) == expected
        frozen[str(directory / name)] = expected
    return doc


original = research / 'evidence/round26/roots.json'
raw, bindings, ordered = verified_study_answers(read(original), original.parent / 'fingerprints')
frozen.update(bindings)
assert len(ordered) == 192512
ordered = [r for r in ordered if r[0] == 'GEOMEAN']
assert len(ordered) == 31744
before = lines(run / 'r26-root-before-core/GEOMEAN.actuals.jsonl')
native = lines(run / 'r26-root-after-core/GEOMEAN.actuals.jsonl')
soft = lines(run / 'r26-root-after-core-soft/GEOMEAN.actuals.jsonl')
assert native == soft and len(native) == len(before) == len(ordered)
r26 = Counter()
for index, ((name, args, answer), old, new) in enumerate(zip(ordered, before, native)):
    assert old['id'] == new['id'] == f'original-service-{index}'
    assert old['args'] == new['args'] == ['0x' + a['bits'] for a in args]
    assert old['not_excel_evidence'] is True and new['not_excel_evidence'] is True
    expected = '0x' + answer['bits']
    was, now = old['calc_func_prediction'] == expected, new['calc_func_prediction'] == expected
    r26['total'] += 1; r26['before_pass'] += was; r26['after_pass'] += now
    r26['gains'] += now and not was; r26['losses'] += was and not now
assert r26 == Counter(total=31744, before_pass=29310, after_pass=31744, gains=2434, losses=0)
for directory in ('r26-root-before-core', 'r26-root-after-core', 'r26-root-after-core-soft'):
    verify_receipt(run / directory)
scale_study = research / 'evidence/round28/root-scales.json'
_, bindings, scales = verified_study_answers(read(scale_study), scale_study.parent / 'root-scales-service')
frozen.update(bindings)
scale_native = lines(run / 'root-original-after-core.jsonl')
scale_soft = lines(run / 'root-original-after-core-soft.jsonl')
initial = lines(run / 'root-original-before.jsonl')
normalized = lines(run / 'root-original-after.jsonl')
assert scale_native == scale_soft and len(scales) == len(scale_native) == 19200
scale_counts = Counter()
for index, ((name, args, answer), old, mid, new) in enumerate(zip(scales, initial, normalized, scale_native)):
    assert name == 'GEOMEAN'
    assert all(r['id'] == f'r28-root-original-{index}' for r in (old, mid, new))
    assert all(r['args'] == ['0x' + a['bits'] for a in args] and r['not_excel_evidence'] is True for r in (old, mid, new))
    expected = '0x' + answer['bits']
    scale_counts['total'] += 1
    scale_counts['initial_pass'] += old['calc_func_prediction'] == expected
    scale_counts['normalized_only_pass'] += mid['calc_func_prediction'] == expected
    scale_counts['after_pass'] += new['calc_func_prediction'] == expected
assert scale_counts == Counter(total=19200, initial_pass=17435, normalized_only_pass=17438, after_pass=19200)
control_study = research / 'evidence/round28/root-normalization.json'
controls, bindings, _ = verified_study_answers(read(control_study), control_study.parent / 'root-normalization-service')
frozen.update(bindings)
reduced_native = lines(run / 'root-normalized-after-core.jsonl')
reduced_soft = lines(run / 'root-normalized-after-core-soft.jsonl')
assert reduced_native == reduced_soft and len(reduced_native) == 6016
for row in reduced_native:
    args = [dict(type='number', bits=a[2:]) for a in row['args']]
    answer = controls['GEOMEAN', json.dumps(args, sort_keys=True)]
    assert row['not_excel_evidence'] is True and row['calc_func_prediction'] == '0x' + answer['bits']
assert frozen == {p: digest(Path(p)) for p in frozen}
receipt = dict(complete=True, r26=r26, scales=scale_counts,
    reduced_controls=6016, all_native_software_outputs_equal=True,
    source=str(Path(__file__).resolve()), hashes=frozen | {str(Path(__file__)): digest(Path(__file__))},
    limits='Complete original numerical obligations for the declared three cohorts. Cohorts overlap and their gains must not be added as independent bank/corpus passes.')
with output.open('x') as stream:
    json.dump(receipt, stream, indent=2); stream.write('\n')
print(json.dumps({k: v for k, v in receipt.items() if k not in ('hashes', 'limits', 'source')}))
