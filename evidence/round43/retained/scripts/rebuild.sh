#!/bin/bash
# rebuild.sh <outname>: bind changed application files into round43-unit.json, build racer, copy binary
set -e
cd ~/projects/DnaCalc/DnaCalc-Research
python3 - <<"PY"
import json, subprocess
main = "/home/ubuntu/projects/DnaCalc/DnaCalc"
u = json.load(open("round43-unit.json"))
rev = u["application_revision"]
changed = subprocess.run(["git", "-C", main, "diff", "--name-only", rev, "--", "crates", "apps", "tools/calc-verdict/src"], capture_output=True, text=True).stdout.split()
if changed:
    u["application_working_blobs"] = {p: subprocess.run(["git", "-C", main, "hash-object", p], capture_output=True, text=True).stdout.strip() for p in changed}
else:
    u.pop("application_working_blobs", None)
json.dump(u, open("round43-unit.json", "w"), indent=2); open("round43-unit.json", "a").write("\n")
print(u.get("application_working_blobs"))
PY
export TMPDIR=~/codex-runs/r43/tmp PATH=~/.cargo/bin:$PATH
mkdir -p $TMPDIR
DNA_RESEARCH_UNIT=round43-unit.json cargo build --release -q -p symbolic-racer 2>&1 | grep -E "^error|panicked" -A7 | head -30 || true
cp target/release/symbolic-racer ~/codex-runs/r43/$1
