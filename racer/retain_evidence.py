"""Generalised compact retention route used by rounds38/39. Original MIT.

Only explicit inputs are retained. Existing payloads are immutable; an index
can grow. Both original and compressed hashes are verified on every update.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def retain(run, destination, inputs):
    run=Path(run).resolve();destination=Path(destination).resolve()
    destination.mkdir(parents=True,exist_ok=True)
    index=destination/'files.json'
    records={r['original_path']:r for r in json.loads(index.read_text())}if index.exists()else{}
    for original in inputs:
        source=Path(original).resolve()
        relative=str(source.relative_to(run))
        raw=source.read_bytes();compressed=relative.endswith('.gz')
        payload=raw if compressed else gzip.compress(raw,mtime=0)
        path=Path('retained')/(relative if compressed else relative+'.gz')
        target=destination/path;target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists():
            if target.read_bytes()!=payload:raise ValueError('changed retained input: '+relative)
        else:
            with target.open('xb')as stream:stream.write(payload)
        records[relative]=dict(path=str(path),original_path=relative,
            source_sha256=hashlib.sha256(raw).hexdigest(),stored_sha256=hashlib.sha256(payload).hexdigest(),
            bytes=len(raw),stored_bytes=len(payload),already_compressed=compressed)
    for record in records.values():
        payload=(destination/record['path']).read_bytes()
        raw=payload if record['already_compressed']else gzip.decompress(payload)
        if hashlib.sha256(payload).hexdigest()!=record['stored_sha256']or hashlib.sha256(raw).hexdigest()!=record['source_sha256']:
            raise ValueError('retention hash differs: '+record['original_path'])
    index.write_text(json.dumps(sorted(records.values(),key=lambda r:r['original_path']),indent=2)+'\n')
    return dict(records=len(records),stored_bytes=sum(r['stored_bytes']for r in records.values()))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('run');parser.add_argument('destination');parser.add_argument('inputs',nargs='+')
    args=parser.parse_args()
    print(retain(args.run,args.destination,args.inputs))
