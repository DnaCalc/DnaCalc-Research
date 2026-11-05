#!/usr/bin/env python3
"""Verify only8 obsolete calculation pins disappear; retain complete statuses."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sys

main, run, output = map(Path, sys.argv[1:])
file = main / 'crates/calc-func/tests/witnesses/GEOMEAN.json'
before_file = run / 'reviewer-geomean-pins-before.json'
old, new = [json.loads(p.read_text()) for p in (before_file, file)]
assert {k: v for k, v in old.items() if k != 'witnesses'} == {k: v for k, v in new.items() if k != 'witnesses'}
assert len(old['witnesses']) == len(new['witnesses']) == 120
changed = []
for a, b in zip(old['witnesses'], new['witnesses']):
    if a != b:
        assert 'calc_func' in a and 'calc_func' not in b
        assert {k: v for k, v in a.items() if k != 'calc_func'} == b
        changed.append(b['id'])
assert len(changed) == 8
changed_path = run / 'bank-after-core-final/changed-outcomes.jsonl.gz'
bank = {r['id']: r for r in map(json.loads, gzip.decompress(changed_path.read_bytes()).splitlines())}
for row in new['witnesses']:
    if row['id'] in changed:
        assert row['args'] == bank[row['id']]['args']
        assert row['expected_bits'] == bank[row['id']]['actual'] == bank[row['id']]['expected']
remaining = [r for r in new['witnesses'] if 'calc_func' in r]
assert len(remaining) == 1 and remaining[0]['calc_func'] == 'fault:NotImplemented'
assert remaining[0]['id'] == 'w111broad-20260929-GEOMEAN-00008'
reconciliation = json.loads((run / 'pin-refresh-reconciliation.json').read_text())
assert reconciliation['only_removed_calculation_pins'] == changed
assert reconciliation['new_sha256'] == hashlib.sha256(file.read_bytes()).hexdigest()
read_table = lambda p: list(csv.DictReader(p.read_text().splitlines(), delimiter='\t'))
baseline_table = run / 'full-status-before-pin-refresh.tsv'
current_table = main / 'crates/calc-func/witness-status.tsv'
a, b = read_table(baseline_table), read_table(current_table)
assert len(a) == len(b) == 574
changed_status = []
for old_row, new_row in zip(a, b):
    assert (old_row['function'], old_row['bank']) == (new_row['function'], new_row['bank'])
    if old_row != new_row:
        assert (new_row['function'], new_row['bank']) == ('GEOMEAN', 'pinned')
        assert new_row['matched'] == '119' and new_row['rows'] == '120' and new_row['worst_ulp'] == '-'
        assert new_row['calc_func_commit'] == next(r['calc_func_commit'] for r in b if r['function'] == 'GEOMEAN' and r['bank'] == 'full')
        changed_status.append((old_row, new_row))
assert len(changed_status) == 1
bank_path = run / 'bank-after-core-final/summary.json'
summary = json.loads(bank_path.read_text())
full = [r for r in b if r['bank'] == 'full']
assert len(full) == len(summary['functions']) == 295
unread = {}
for row in full:
    actual = summary['functions'][row['function']]
    assert int(row['matched']) == actual['matched']
    extra = int(row['rows']) - actual['rows']
    assert extra >= 0
    if extra:
        unread[row['function']] = extra
assert unread == {'MINVERSE': 112, 'POISSON.DIST': 4}
assert sum(unread.values()) == summary['unread']
assert sum(int(r['rows']) for r in full) == summary['total']
assert sum(int(r['matched']) for r in full) == summary['matched']
files = (file, before_file, baseline_table, current_table, changed_path, bank_path,
         run / 'pin-refresh-reconciliation.json', run / 'geomean-pin-refresh.log', Path(__file__))
receipt = dict(complete=True, unchanged_witness_rows=120, only_removed_calculation_pins=changed,
    all8_now_match_original_Excel=True, remaining_refusal_unchanged=True,
    all295_full_status_rows_exactly_restored=True, refreshed_pinned_status=changed_status,
    full_status_total=summary['total'], full_status_matched=summary['matched'], unread_by_function=unread,
    source=str(Path(__file__).resolve()), hashes={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
with output.open('x') as stream:
    json.dump(receipt, stream, indent=2); stream.write('\n')
print(json.dumps({k: v for k, v in receipt.items() if k not in ('hashes', 'source', 'only_removed_calculation_pins', 'refreshed_pinned_status')}))
