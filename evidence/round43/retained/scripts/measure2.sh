#!/bin/bash
set -u
cd ~/projects/DnaCalc/DnaCalc
git checkout -- crates/calc-func/src/ported/discrete_dist.rs crates/calc-func/src/ported/special_math_common.rs
~/codex-runs/r43/measure.sh head2
cp ~/codex-runs/r43/beta-change/discrete_dist.rs ~/codex-runs/r43/beta-change/special_math_common.rs crates/calc-func/src/ported/
touch crates/calc-func/src/ported/discrete_dist.rs crates/calc-func/src/ported/special_math_common.rs
~/codex-runs/r43/measure.sh after-beta
