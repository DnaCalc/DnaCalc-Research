#!/usr/bin/env python3
"""Complete bank bit/digest audit and explicit remaining GEOMEAN refusal."""
import gzip
import hashlib
import json
from pathlib import Path
import sys

research, run, output = map(Path, sys.argv[1:])
sys.path.insert(0, str(research / 'checker'))
from exact import BINARY64 as B, EXTENDED64 as E, binary64
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
frozen = {}


def read(path):
    frozen[str(path)] = digest(path)
    return json.loads(path.read_text())


left, right = run / 'bank-before-core-final', run / 'bank-after-core-final'
before, after = read(left / 'summary.json'), read(right / 'summary.json')
assert before['complete'] and after['complete']
for key in ('decoded_rows', 'total', 'unread', 'identity_answer_environment_sha256'):
    assert before[key] == after[key]
assert (after['total'], before['matched'], after['matched'], after['failed'], after['unread']) == (3910007, 3769062, 3769119, 140888, 116)
assert after['changed_functions'] == ['GEOMEAN'] and after['gains'] == after['changed_outcomes'] == 57 and after['losses'] == 0
for name, old in before['functions'].items():
    new = after['functions'][name]
    assert old['rows'] == new['rows']
    assert new['matched'] - old['matched'] == (57 if name == 'GEOMEAN' else 0)
    assert new['gains'] == (57 if name == 'GEOMEAN' else 0) and new['losses'] == 0
assert [f for f in before['actual_sha256_by_function'] if before['actual_sha256_by_function'][f] != after['actual_sha256_by_function'][f]] == ['GEOMEAN']
old_bits, new_bits = (left / 'matched.bits').read_bytes(), (right / 'matched.bits').read_bytes()
for path, bits, summary in ((left / 'matched.bits', old_bits, before), (right / 'matched.bits', new_bits, after)):
    assert hashlib.sha256(bits).hexdigest() == summary['match_bits_sha256']
    frozen[str(path)] = digest(path)
    assert int.from_bytes(bits, 'little') >> summary['decoded_rows'] == 0
old, new = int.from_bytes(old_bits, 'little'), int.from_bytes(new_bits, 'little')
assert old.bit_count() == before['matched'] and new.bit_count() == after['matched']
assert (~old & new).bit_count() == 57 and (old & ~new).bit_count() == 0
changed_path = right / 'changed-outcomes.jsonl.gz'
changed = [json.loads(line) for line in gzip.decompress(changed_path.read_bytes()).splitlines()]
assert len(changed) == 57
assert all(r['function'] == 'GEOMEAN' and r['before_matched'] is False and r['actual'] == r['expected'] for r in changed)
frozen[str(changed_path)] = digest(changed_path)
for summary in (before, after):
    for path, expected in summary['inputs_sha256'].items():
        assert digest(Path(path)) == expected
        frozen[path] = expected
trace_path = run / 'geomean-final-trace.jsonl'
trace = [json.loads(line) for line in trace_path.read_text().splitlines()]
assert len(trace) == 1225 and all(r['function'] == 'GEOMEAN' for r in trace)
failures = [r for r in trace if r['actual'] != r['expected']]
remaining = read(run / 'geomean-final-remaining.json')
assert remaining == dict(total=1225, mismatches=failures) and len(failures) == 1
assert hashlib.sha256(b''.join((r['actual'] + '\n').encode() for r in trace)).hexdigest() == after['actual_sha256_by_function']['GEOMEAN']
failure = failures[0]
assert failure['args'] == ['0x0010000000000000'] * 3
assert failure['expected'] == '0x0010000000000000' and failure['actual'] == 'fault:NotImplemented'
product = binary64('3ff0000000000000')
for arg in failure['args']:
    product = B.round(E.round(product * binary64(arg)))
assert product == 0
frozen[str(trace_path)] = digest(trace_path)
assert frozen == {path: digest(Path(path)) for path in frozen}
receipt = dict(complete=True, total=after['total'], before_pass=before['matched'], after_pass=after['matched'],
    inclusive_fail=after['failed'], unread=after['unread'], gains=57, losses=0,
    changed_functions=['GEOMEAN'], geomean=dict(total=1225, before_pass=1167, after_pass=1224, failure=failure),
    remaining_classification='Unchanged out-of-scope fallback: three minimum-normal inputs make the stored product zero, selecting the preexisting fallback. The typed refusal remains one inclusive failure under the parent package.',
    source=str(Path(__file__).resolve()), hashes=frozen | {str(Path(__file__)): digest(Path(__file__))})
with output.open('x') as stream:
    json.dump(receipt, stream, indent=2); stream.write('\n')
print(json.dumps({k: v for k, v in receipt.items() if k not in ('hashes', 'source', 'geomean')}))
