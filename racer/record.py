"""Run research commands with a new, never overwritten transcript."""
import json
from pathlib import Path
import resource
import subprocess
import sys
import time
import hashlib
import shutil

log = Path(sys.argv[1])
log.parent.mkdir(parents=True, exist_ok=True)
start = time.monotonic()
with log.open('x') as f:
    executable=Path(sys.argv[2])
    if not executable.is_file():executable=Path(shutil.which(sys.argv[2])or sys.argv[2])
    binary_hash=None
    if executable.is_file():
        hasher=hashlib.sha256()
        with executable.open('rb')as binary:
            for block in iter(lambda:binary.read(1024*1024),b''):hasher.update(block)
        binary_hash=hasher.hexdigest()
    f.write(json.dumps({'command': sys.argv[2:],'executable_sha256':binary_hash}) + '\n')
    f.flush()
    result = subprocess.run(sys.argv[2:], stdout=f, stderr=subprocess.STDOUT)
summary = dict(returncode=result.returncode, seconds=time.monotonic()-start,
               child_peak_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
               log=str(log),executable_sha256=binary_hash)
with log.with_suffix(log.suffix + '.resource.json').open('x') as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary))
print('\n'.join(log.read_text().splitlines()[-8:]))
sys.exit(result.returncode)
