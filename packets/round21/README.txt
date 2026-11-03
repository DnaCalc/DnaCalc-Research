R21 hypergeometric questions: 20 scalar calls; no expected answers.
From the published DnaCalc Windows checkout, use the existing witness reader:
  cargo xtask witness read --func HYPGEOM.DIST --inputs <packet>/requests-HYPGEOM.DIST.json
Do not use --new-only: recorded center controls are intentional. Only requests-*.json
are reader inputs. Predictions are application calculations, never Excel readings.
Preserve the canonical WitnessSet in OxFunc-Evidence, requested/submitted/Value2
bit audits, worker settings, original requests/results, frozen files, packet hashes
and partial attempts. Return them under ~/codex-runs/excel-results-r21.
No second oracle or capture lineage is introduced. Missing/unread outcomes fail.

The same packet asks 204 scalar formula cells on catalogue sheet dcj9k_r21:
120 LINEST/LOGEST selected outputs and84 operand echoes. They extend the one
catalogue oracle; the existing generator registers sheet-dcj9k_r21.txt.
From the published DnaCalc checkout on Windows:
  cargo xtask capture --pending --list
  cargo xtask capture --pending
Keep sections/dcj9k_r21.capture.json and the original worker capture outputs.
Do not place application predictions into a capture or a WitnessSet. The operand
echoes record Excel's literal parser results; requested binary64 arithmetic and
parsed literals need not have identical bits. These are discovery questions,
not fitted rank rules or a coefficient program. Return all failures and attempts.
