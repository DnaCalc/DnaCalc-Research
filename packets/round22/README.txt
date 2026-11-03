R22 BESSELI scalar questions; no Excel answers included.
From the published Windows DnaCalc checkout:
  cargo xtask witness read --func BESSELI --inputs <packet>/requests-BESSELI.json
Do not use --new-only: independently recorded controls intentionally repeat.
Only requests-BESSELI.json is reader input. Application predictions and route
audits are calculations, never Excel answers or expected values.
Preserve canonical WitnessSets in OxFunc-Evidence, requested/submitted/Value2
bit audits, worker settings, original requests/results, file hashes and partial
attempts. Return those files and the unchanged packet under
~/codex-runs/excel-results-r22. Missing/unread outcomes remain failures.

Do not append research questions to the captured catalogue workbook. R21's
204 cells now live in crates/calc-func/tests/fixtures/dcj9k_r21/oracle.xlsx;
fn2 questions live in fixtures/fn2_questions/oracle.xlsx. The original catalogue
has been restored byte-for-byte and needs no recapture for those1633 old cells.
Use cargo xtask capture --pending for the separate pending workbooks.
This supersedes only R21's instruction to capture its questions in catalogue;
the historical R21 request/answer formats and hashes remain unchanged.
