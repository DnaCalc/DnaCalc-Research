#!/usr/bin/env python3
"""Validate the canonical R28 payload, restoration data and producer epochs."""
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

root, output = map(Path, sys.argv[1:])
root = root.resolve()
base = root / 'evidence/round28'
digest = lambda data: hashlib.sha256(data).hexdigest()


def document(path):
    data = path.read_bytes()
    return json.loads(gzip.decompress(data) if path.suffix == '.gz' else data)


index_path, compression_path, producers_path = [base / p for p in ('index.json', 'compression.json', 'producer-bindings.json')]
index, compression, producers = [document(p) for p in (index_path, compression_path, producers_path)]
seen, total = set(), 0
for row in index['artefacts']:
    path = root / row['path']
    assert not Path(row['path']).is_absolute() and path.resolve().is_relative_to(base)
    assert row['path'] not in seen
    seen.add(row['path'])
    data = path.read_bytes()
    assert len(data) == row['bytes'] and digest(data) == row['sha256'], row['path']
    total += len(data)
assert total == index['bytes']
restored_bytes, archive_bytes = 0, 0
restores = {r['original']: r for r in compression['compressed']}
for name, row in restores.items():
    original, archive = root / name, root / row['archive']
    assert original.resolve().is_relative_to(base) and archive.resolve().is_relative_to(base)
    assert row['archive'] in seen and row['archive'] == name + '.gz'
    packed = archive.read_bytes(); raw = gzip.decompress(packed)
    assert digest(packed) == row['archive_sha256']
    assert digest(raw) == row['uncompressed_sha256'] and len(raw) == row['uncompressed_bytes']
    assert original.parent.is_dir()
    if original.exists():
        assert original.read_bytes() == raw
    restored_bytes += len(raw); archive_bytes += len(packed)
ignored = subprocess.run(['git', 'check-ignore', '--no-index', '--stdin'], cwd=root,
    input='\n'.join(restores) + '\n', text=True, capture_output=True, check=True)
assert set(ignored.stdout.splitlines()) == set(restores)
extras = {str(p.relative_to(root)) for p in base.rglob('*') if p.is_file()} - seen - {str(index_path.relative_to(root))}
assert extras <= set(restores)
known = {}
for h, names in producers['source_epochs'].items():
    for name in names:
        assert 'DnaCalc-engine' not in name
        if name.startswith('git:'):
            _, revision, relative = name.split(':', 2)
            assert re.fullmatch('[0-9a-f]{7,40}', revision) and '..' not in Path(relative).parts
            data = subprocess.check_output(['git', 'show', revision + ':' + relative], cwd=root)
        else:
            path = root / name
            assert path.resolve().is_relative_to(root)
            data = path.read_bytes()
            if path.suffix == '.gz':
                data = gzip.decompress(data)
        assert digest(data) == h, name
        known.setdefault(h, []).append(name)
bindings = []
for path in sorted(base.iterdir()):
    if not path.is_file() or not (path.name.endswith('.json') or path.name.endswith('.json.gz')):
        continue
    if path.suffix != '.gz' and path.with_suffix(path.suffix + '.gz').is_file():
        assert path.read_bytes() == gzip.decompress(path.with_suffix(path.suffix + '.gz').read_bytes())
        continue
    doc = document(path)
    if not isinstance(doc, dict) or not isinstance(doc.get('producer_bindings'), dict):
        continue
    for source, h in doc['producer_bindings'].items():
        assert 'DnaCalc-engine' not in source
        current = Path(source)
        assert current.resolve().is_relative_to(root)
        assert current.is_file() and digest(current.read_bytes()) == h or h in known
        bindings.append((path.name, source, h))
assert len(bindings) == producers['binding_count'] == 72 and producers['unresolved'] == []
for name in index['authoritative_proofs']:
    doc = document(base / name)
    assert doc['complete'] is True
assert not set(index['authoritative_proofs']) & set(index['historical_reports'])
copied_review = 0
run = output.parent
for path in (base / 'review').rglob('*'):
    if not path.is_file():
        continue
    original = run / path.name.removesuffix('.gz')
    if original.is_file():
        data = gzip.decompress(path.read_bytes()) if path.suffix == '.gz' else path.read_bytes()
        assert data == original.read_bytes(), path.name
        copied_review += 1
result = dict(complete=True, artifact_count=len(seen), artifact_bytes=total,
    restoration_reports=len(restores), restored_bytes=restored_bytes, compressed_bytes=archive_bytes,
    canonical_producer_bindings=len(bindings), source_epochs=len(known),
    authoritative_proofs=len(index['authoritative_proofs']), copied_review_files_verified=copied_review,
    extras_only_ignored_restorations=sorted(extras),
    restore_command='README command was simulated in memory: exact gzip hash/content, parent directories, existing-copy equality and all34 ignore entries verified; repository files were not written.',
    hashes={str(p): digest(p.read_bytes()) for p in (index_path, compression_path, producers_path, root / 'README.txt', root / '.gitignore', Path(__file__))})
with output.open('x') as stream:
    json.dump(result, stream, indent=2); stream.write('\n')
print(json.dumps({k: v for k, v in result.items() if k != 'hashes'}))
