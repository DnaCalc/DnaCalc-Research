Round 14 variance surface and component audit: 44 requests, zero Excel answers.

The first six distinct triples in the corpus selection supply 36 calls:
VAR.S, VAR.P, STDEV.S, STDEV.P, DEVSQ and AVERAGE on each tuple.
Four fixed one-ULP steps toward positive infinity add eight VAR.S/DEVSQ calls.
No complete program survives the preceding finite searches. These questions
audit transport and components; they do not validate a replacement kernel.

The source captures called functions with one array argument. These requests
use three direct numeric Value2 arguments. They preserve the numeric tuple,
and deliberately change the argument surface. Same typed-input identity and
unchanged private Excel operands are not asserted.

Use the existing reader on Windows from the DnaCalc checkout, for each of
VAR.S, VAR.P, STDEV.S, STDEV.P, DEVSQ and AVERAGE:
  cargo xtask witness read --func VAR.S --inputs ..\DnaCalc-Research\packets\round14-variance\requests-VAR.S.json
Substitute the function in both places. Do not pass --new-only.

Requests contain no answers. application-predictions-* contain current product
predictions only; mathematical-controls.json contains exact mathematical
values only. Neither is an Excel reading. Missing, refused or changed-input
responses remain failed obligations. No annuity question or reading is added.

Keep the reader's workbook compatibility/settings metadata with its answers.
Those settings may differ from the original array captures; a changed result
alone does not isolate an argument-surface effect.
