#!/bin/bash
# usage: measure.sh <label>
set -u
export TMPDIR=$HOME/codex-runs/r43/tmp PATH=$HOME/.cargo/bin:$PATH DNACALC_PAINT_PYTHON=python3
cd $HOME/projects/DnaCalc/DnaCalc
L=$1; O=$HOME/codex-runs/claude-r43/$L
mkdir -p $O
git rev-parse HEAD > $O/head.txt; git diff > $O/diff.patch
( time cargo xtask witness judge --bank full --report $O/bank.json ) > $O/bank.log 2>&1
( time cargo xtask corpus --oracle --no-capture --out $O/corpus ) > $O/corpus.log 2>&1
echo done > $O/DONE
