R16 Excel questions: 250 scalar calls and 810 catalogue formula cells

On Windows, with the published DnaCalc master containing sheet dcj9k_r16:
  powershell -NoProfile -ExecutionPolicy Bypass -File RUN-R16.ps1 -DnaCalcCheckout C:\Work\DnaCalc\DnaCalc-win
Supply -EvidenceCheckout if OxFunc-Evidence is not beside that checkout.
Supply a new -OutputDirectory for another attempt. Preserve partial attempts.

The runner verifies packet and section hashes. It uses the existing catalogue
capture and witness reader, never another workbook generator or Excel adapter.
capture --pending reads dcj9k_r16 and any other pending sections. Adding a sheet
also changes the workbook context of prefix_metadata and dcuar_cell, which call
SHEETS()/INFO("numfile"); their 1,633 old cells require a new context reading.
They are not 1,633 new research questions. No saved Excel answer is regenerated
here. The standard capture tool controls and records its own Excel processes.

The scalar calls use exact numeric Value2 inputs through the standard reader.
Do not use --new-only: this packet deliberately repeats earlier calls.
The catalogue compares separate references, ranges, horizontal/vertical array
literals, transposed ranges and separate literals. Operand echo cells retain
the worksheet parser's results. A changed echo is an observation of different
inputs; it cannot silently become an answer to the intended binary64 tuple.
The captures record their calculation and compatibility settings; differences
from the reader's new workbook are not attributed to the argument surface alone.

Bessel questions use fixed adjacent/binary grids, orders and parity controls,
plus independent EXP/SIN/COS and stored/expression phase calls. These standalone
calls are not observations of Bessel's private leaves.
Regression questions vary fixed binary perturbations and independent scales,
with exact dependency and zero-column controls. A surviving rank-only policy
is not a coefficient program and does not authorize a fitted tolerance.
GEOMEAN questions repeat the 1,5,3 product under permutations, counts and surfaces.

Only requests-*.json go to witness read. catalogue-questions.json documents the
registered section and identifies its cells; the section file in this packet
mirrors that text, not a second oracle workbook. application-predictions-* are
current application calculations, never Excel readings. No R15 annuity answer
is available or supplied. The new questions have no expected-answer fields.

Return the complete output directory, including packet hashes, canonical
WitnessSets, frozen worker predictions, raw request/result/scripts, section
captures and capture work directories. Unread, missing, changed-input and
failed runs remain failed obligations. Keep original files byte for byte.
