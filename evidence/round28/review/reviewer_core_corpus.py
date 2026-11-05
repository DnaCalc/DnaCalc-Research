#!/usr/bin/env python3
"""Independent complete corpus delta and exact original84 GEOMEAN scope audit."""
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import sys
import zipfile
import xml.etree.ElementTree as ET

run, output = map(Path, sys.argv[1:])
before, after = run / 'corpus-before-core-final', run / 'corpus-after-core'
delta_path = run / 'core-corpus-delta.json'
original_path = run.parent.parent / 'engine13/original-root-causes.json'
if not original_path.is_file():
    original_path = Path('/home/ubuntu/codex-runs/engine13/original-root-causes.json')
frozen = {}


def read(path):
    raw = path.read_bytes()
    frozen[str(path)] = hashlib.sha256(raw).hexdigest()
    return json.loads(raw)


summary_before, summary_after = read(before / 'summary.json'), read(after / 'summary.json')
assert summary_before['run']['corpus_inventory_sha256'] == summary_after['run']['corpus_inventory_sha256']
assert [(w['name'], w['sha256'], w['input_identity']) for w in summary_before['workbooks']] == [
       (w['name'], w['sha256'], w['input_identity']) for w in summary_after['workbooks']]
declared = read(delta_path)
reports = sorted(p for p in (before / 'workbooks').rglob('report.json')
                 if not p.parent.name.startswith('attempt-'))
expected_paths = {str(p.relative_to(before)) for p in reports}
assert expected_paths == {str(Path(p['before']).relative_to(before)) for p in declared['pairs']}
assert expected_paths == {str(p.relative_to(after)) for p in (after / 'workbooks').rglob('report.json')
                         if not p.parent.name.startswith('attempt-')}
key = lambda row: (row['sheet'], row['row'], row['column'])
counts, classes = Counter(), Counter()
scored_catalogue = defaultdict(list)
extras = []
scored_changes = []
for old_report in reports:
    relative = old_report.relative_to(before)
    new_report = after / relative
    old, new = read(old_report), read(new_report)
    old_eval, new_eval = read(old_report.with_name('evaluations.json')), read(new_report.with_name('evaluations.json'))
    old_triage, new_triage = read(old_report.with_name('triage.json')), read(new_report.with_name('triage.json'))
    a, b = [{key(r): r for r in d['cells']} for d in (old, new)]
    x, y = [{key(r): r for r in d['evaluations']} for d in (old_eval, new_eval)]
    t, u = [{key(r): r for r in d['cells']} for d in (old_triage, new_triage)]
    assert a.keys() == b.keys() and x.keys() == y.keys()
    capture = None
    if old_report.parent.name.startswith('session-'):
        cp = old_report.parent.parent / (old_report.parent.name + '.capture.json')
        np = new_report.parent.parent / (new_report.parent.name + '.capture.json')
        capture = read(cp)
        assert hashlib.sha256(np.read_bytes()).hexdigest() == frozen[str(cp)]
        frozen[str(np)] = frozen[str(cp)]
    counts['sessions'] += 1
    counts['formula_cells'] += len(a)
    counts['evaluations'] += len(x)
    for k in a:
        policy = bool(t[k].get('external_links') or t[k].get('circular_rule'))
        assert policy == bool(u[k].get('external_links') or u[k].get('circular_rule'))
        assert (a[k]['verdict'] == 'awaiting_capture') == (b[k]['verdict'] == 'awaiting_capture')
        if policy:
            counts['external_or_circular'] += 1
        else:
            was, now = a[k]['verdict'] == 'matched', b[k]['verdict'] == 'matched'
            counts['before_matched'] += was
            counts['after_matched'] += now
            counts['gains'] += now and not was
            counts['losses'] += was and not now
            if was != now:
                scored_changes.append(dict(workbook=str(relative), cell=list(k), before=a[k], after=b[k],
                                           original_evaluation=x.get(k), current_evaluation=y.get(k)))
        if 'calc-func_catalogue_oracle.xlsx-' in str(relative):
            scored_catalogue[k].append((relative, a[k], b[k], x.get(k), y.get(k), capture))
    for k in x:
        assert {j: v for j, v in x[k].items() if j != 'returned'} == {j: v for j, v in y[k].items() if j != 'returned'}
        if x[k].get('returned') == y[k].get('returned'):
            continue
        counts['changed_values'] += 1
        if k not in a:
            classification = 'outside_scored_session_inventory'
            extras.append((relative, k, x[k], y[k]))
        else:
            assert a[k]['verdict'] != 'matched' and b[k]['verdict'] == 'matched'
            classification = 'gained_exact_match'
        classes[classification] += 1
counts['awaiting_capture'] = 0
counts['workbooks'] = len(summary_before['workbooks'])
assert counts == declared['counts']
assert classes == declared['calculated_change_classes']
assert scored_changes == declared['changes']
assert counts['before_matched'] == summary_before['totals']['counts']['pass']
assert counts['after_matched'] == summary_after['totals']['counts']['pass']
original = read(original_path)
scope = [r for r in original['cells'] if r['group'] == 'geomean-accumulation']
assert len(scope) == original['groups']['geomean-accumulation']['count'] == 84
package = Path('/home/ubuntu/projects/DnaCalc/DnaCalc/crates/calc-func/tests/fixtures/catalogue/oracle.xlsx')
frozen[str(package)] = hashlib.sha256(package.read_bytes()).hexdigest()
source_formulas = {}
ns = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
with zipfile.ZipFile(package) as archive:
    rels = {r.attrib['Id']: r.attrib['Target'] for r in ET.fromstring(archive.read('xl/_rels/workbook.xml.rels'))}
    for sheet in ET.fromstring(archive.read('xl/workbook.xml')).find('m:sheets', ns):
        if sheet.attrib['name'] not in {r['sheet'] for r in scope}:
            continue
        target = rels[sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']]
        target = target.lstrip('/') if target.startswith('/') else 'xl/' + target
        for cell in ET.fromstring(archive.read(target)).findall('.//m:c', ns):
            formula = cell.find('m:f', ns)
            if formula is not None:
                source_formulas[sheet.attrib['name'], cell.attrib['r']] = formula.text
scope_answers = {}
scope_records = []
for row in scope:
    k = row['sheet'], row['row'], row['column']
    items = scored_catalogue[k]
    assert len(items) == 1
    relative, old, new, x, y, capture = items[0]
    want = row['tool_cell']['difference']['expected']['number']
    assert old['verdict'] == 'mismatched' and new['verdict'] == 'matched'
    assert y['returned']['value'] == dict(number=want)
    assert capture['source'] == row['source']
    formula = next(f for f in capture['formulas'] if key(f) == k)
    reading = next(o for o in capture['observations'] if key(o) == k)
    # Excel's returned Formula2 text normalizes prefixes, .0 and exponent signs.
    # The immutable source package retains the exact original source spelling.
    assert row['source'] == 'sha256:' + frozen[str(package)]
    assert source_formulas[row['sheet'], row['cell']] == row['formula']
    assert reading['excel_raw']['number'] == want
    scope_answers[k] = want
    scope_records.append(dict(cell=list(k), session=str(relative), matched_original_number=want))
extra_counts = Counter()
for relative, k, old, new in extras:
    assert 'calc-func_catalogue_oracle.xlsx-' in str(relative) and k in scope_answers
    assert new['returned']['value'] == dict(number=scope_answers[k])
    extra_counts[k] += 1
assert len(extras) == 336 and set(extra_counts) == set(scope_answers)
assert all(n == 4 for n in extra_counts.values())
assert frozen == {path: hashlib.sha256(Path(path).read_bytes()).hexdigest() for path in frozen}
receipt = dict(complete=True, counts=counts, change_classes=classes,
    original_bead_cells=84, original_bead_matched=84, original_scope=scope_records,
    outside_inventory_classification='Exactly4 additional catalogue calculation sessions for each original84 cell; values now equal their original recorded Excel bits, and are not counted as extra passes.',
    remaining_corpus_GEOMEAN_mismatch_classes=[r for r in summary_after['root_causes']['mismatches'] if r['function'] == 'GEOMEAN'],
    source=str(Path(__file__).resolve()), hashes=frozen | {str(Path(__file__)): hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
with output.open('x') as stream:
    json.dump(receipt, stream, indent=2)
    stream.write('\n')
print(json.dumps({k: v for k, v in receipt.items() if k not in ('hashes', 'original_scope', 'source')}))
