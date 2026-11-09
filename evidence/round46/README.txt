R46: ETS reset counts and short-history seasonality; YIELD/discrete refutations

Original MIT research by GPT-6.1 Sol (Codex), dna-firehorse, implementer.
Sources: original Excel20527 question-service answers, IEEE754, elementary
ETS/raw-moment and zero-coupon mathematics, current cleared repository code.
No tier3/4 implementation source or excluded checkout was read or changed.

Application changes: ETS uses raw trend moments with separate tail/divisor
counts, stored-PC64 operations, and no longer refuses the measured
minimum-alpha forecast region. SEASONALITY returns0 for2..4 prepared points.
No parameter point, fitted real constant or new arithmetic piece was added.
Full seasonal search, short/inexact initialization and longer detection stay
unidentified and retain counted refusals. See design/ets-search.txt.
An offline erratum corrects the unimplemented length4 seasonal candidate:
its original forecast distance added two extra trend steps. The corrected
distance is index-phase, proven with exact rational arithmetic and every
zero-horizon result equaling the last observation. It still contradicts
368/384 forecasts (original diagnostic incorrectly counted380); STAT560/672
and CI276/288 contradictions are unchanged. Original predictions and answers
remain immutable; short-seasonal-corrected-analysis.json gives the correction.

The 1,111 initial corpus refusals are partitioned by first-refusal tracing:
615 automatic detection/search,204 short initialization,180 seasonal
innovations,112 seasonal history/cycle limits. refusal-census.json retains
all overlapping history/option/missing-point/duplicate dimensions and IDs.
Tracing used temporary owner strings; the final source retains dc-j9k.
Original/candidate/final workbook outputs and their hashes are on R2.

Reset identification: repeated adjacent29-point histories, exact cash scales,
and trend-shift interval bisection identify the count discrepancy. The odd
width sqrt(182/183) is derived from integer counts, not a fitted constant.
The selected stored-PC64 program fits54 retained boundaries,74 earlier reset
decisions,512 fresh triples,20 targeted forecasts and460 independent forecasts
from115 qualifying histories. The latter stopped at its fixed20,000-history
generation limit. Individual stores are nonunique (1,228 surviving programs);
this validates the program without claiming unique private instructions.
Search loads and tie order remain undetermined; the product table is unchanged.

Final R46 registered replay:5,807/7,347 matched,1,540 refused,zero numeric
mismatches. This is mixed discovery/holdout coverage, not7,347 fresh inputs.
All6,735 original IDs/arguments/expectations reconcile:1,397 gains/zero losses;
1,338 previously refused readings now exact. All612 later readings match.
Application-published repeats every prediction after a unit-metadata refresh;
the final ETS source blob is d3260fae (full blob in the compiled provenance).

Wider R39/R40 replay:61,983/84,998 matched;22,547 refused and468 numeric
mismatches remain counted failures. Compared with final R40 receipts,159
gains/zero losses:155 prior refusals and4 STAT beta values become exact.
retained-transitions.json explicitly gives final-receipt precedence. The first
draft's160 count reused one older CI baseline and is superseded; no input or
expected answer changed. These older questions are discovery evidence.

YIELD:5,005 cash-exponent questions map small/large cash limiting answers.
All100 declared Newton/secant seed programs and150 bracket programs fail a
necessary iterate-membership test on12 dyadic controls. No stop/seed constant
was fitted. Analytic price/stores are conditional;16-ULP padding does not
bound another private price implementation. POISSON:516 adjacent controls
equal the public chi-square tail.
BINOM:221/774 public-beta coupling contradictions all have inexact complements;
520/520 exact-complement controls agree. Retaining both probability operands
in the current bratio body still gives only6/774 exact answers. No financial
solver or distribution kernel changed. See design/yield-stop.txt.

Main light check:707 passed,0 failed,3 ignored,1 not run; workspace check and
portable wasm compilation passed. Research:154 Rust passed,2 ignored;
165 Python passed. No Windows/browser sweep or CI dispatch is claimed.
UBS:0 critical;143 warnings were reviewed (bounded histories, deliberate exact
comparisons and test assertions/clones). New findings are test clones/asserts
and bounded casts. ubs-review.json retains the review, not a blanket clean claim.
Research's staged evidence has no supported active source files; the session
Python drivers were scanned separately:0 critical,35 reviewed warnings. A
fixed rclone copy/check argv was initially flagged as executable selection;
the reviewed annotation is in sources/archive-reviewed.py.gz. The original
archiving source epoch stays immutable in archive.py.gz. No behavior changed.
Early attempts failed on a stale research unit, a disabled route bridge and
the inherited full /tmp target. Their logs remain; explicit current-unit builds
and the clone's normal target succeeded. No failed attempt supplied an answer.

Service:155 audited replies,14,932 calls; all use interactive priority and retain
original frames, exact echoes, environments and replies. Every reply is checked
by the existing fingerprint.accept route. No wait exceeded3 seconds; maximum
0.863584s. Rejected initial shell truncation requests were replaced with file
patches/append-only logs; no approval override or service restart was used.

files.json uses the existing retain_evidence route and records original/stored
SHA256 for every small payload. Large raw outputs are in
dnacalc-evidence/research/round46/20261109T121141Z-dna-firehorse:
2 files,25,621,257 bytes;4,225 verified archive members,748,955,968 raw bytes.
Checksum copy and one-way check passed; r2-receipt.json gives both file hashes.
Final bank/corpus publication details are in measurement-bindings.json.

DnaCalc implementation commit: bb2873d321106ed5d765e5e2c7754994d8fe0c90.
