R20 discrete distribution questions: 36 BINOM.DIST calls and8 HYPGEOM.DIST calls.
From the published DnaCalc Windows checkout, use the existing witness reader:
  cargo xtask witness read --func BINOM.DIST --inputs <packet>/requests-BINOM.DIST.json
  cargo xtask witness read --func HYPGEOM.DIST --inputs <packet>/requests-HYPGEOM.DIST.json
Do not use --new-only. Center points intentionally repeat recorded controls.
Only requests-*.json are reader inputs. Predictions and route audits are current
application calculations, never Excel answers or expected values.

Preserve the full requested/submitted/Value2 bit audits, worker settings,
canonical WitnessSets, original worker requests/results and frozen predictions.
Return those files, the unchanged packet, file hashes and partial attempts in
~/codex-runs/excel-results-r20. Preserve canonical readings in OxFunc-Evidence.
Input changes, missing/unread results and automation failures remain failures.
This packet uses the established witness reader and capture lineage only.
