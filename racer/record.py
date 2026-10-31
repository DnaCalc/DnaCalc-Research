"""Run research commands with a new, never overwritten transcript."""
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

log = Path(sys.argv[1])
log.parent.mkdir(parents=True, exist_ok=True)
start = time.monotonic()
with log.open('x') as f:
    f.write(json.dumps({'command': sys.argv[2:]}) + '\n')
    f.flush()
    result = subprocess.run(sys.argv[2:], stdout=f, stderr=subprocess.STDOUT)
summary = dict(returncode=result.returncode, seconds=time.monotonic()-start,
               child_peak_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
               log=str(log))
with log.with_suffix(log.suffix + '.resource.json').open('x') as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary))
print('\n'.join(log.read_text().splitlines()[-8:]))
sys.exit(result.returncode)
