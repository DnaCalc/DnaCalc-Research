R22 regression workbook questions:105 formulas,27 selected LINEST outputs and78 operand echoes.
These are discovery questions, not Excel answers. The original scalar R22 packet stays frozen.
The section belongs only in crates/calc-func/tests/fixtures/dcj9k_r22/oracle.xlsx.
It must never be appended to the captured catalogue. The existing
make-catalogue-fixtures.ps1 generates this separate workbook from the identical
committed sheet-dcj9k_r22.txt. From the published Windows checkout:
  cargo xtask capture --pending --list
  cargo xtask capture --pending
Preserve all original worker outputs, saved packages, settings, attempts and
canonical captures; return them with the unchanged packet under
~/codex-runs/excel-results-r22. Predictions are application calculations only.
The original3x3 center is an intentional control. Reversed rows, swapped columns
and a fixed dyadic column scale ask about arithmetic invariance; no Excel
invariance is assumed. Operand echoes measure literal-parser bits separately.
No coefficient, threshold, rank tolerance or per-row arithmetic choice is fitted.
