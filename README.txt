DNA Calc arithmetic program research

This separate native-only research workspace interprets candidate data. Its code
does not enter the application. The first unit is PMT's finite discount combine,
specified by DnaCalc/docs/plan/SYMBOLIC_SEARCH.md and pmt-r6-unit.json (dc-j9k).

The unit records the implementer role, permitted and banned sources, exact source
revision and numeric file hashes, bounded grammar, model semantics and resources.
The arithmetic subset adapts OxFunc's MIT racer at 05196c394f5400eb329c79b01b4a8d2bac7c228b:
smart-fuzzer/tools/calc_graph_racer/src/{dsl,enumerate,eval,scheduler}.rs.
The MIT notice is in LICENSE. No legacy numeric substrate is copied. build.rs
compiles DnaCalc's current cleared numeric operations with visibility/path changes
only, and refuses source changes relative to the unit's frozen application source.
Keep this repository beside the DnaCalc checkout, as DnaCalc-Research.

Build from this directory with Rust 1.94.1:
  CARGO_TARGET_DIR=target cargo build --release --offline -p symbolic-racer

Evidence comes exclusively from the existing DnaCalc witness loader and verdict
encodings. From the DnaCalc sibling directory:
  cargo xtask witness judge --func PMT --bank full --trace R6/application-before.jsonl

Replace R6 below by an outside-repository results directory. Prior gaps are the
round-3 pmt-quotient-gap-certificates.jsonl. From this research directory:
  target/release/symbolic-racer prepare R6/application-before.jsonl PRIOR_GAPS R6/frozen
  target/release/symbolic-racer self-test R6/frozen/scoped.jsonl
  target/release/symbolic-racer run R6/frozen/scoped.jsonl R6/native hard
  target/release/symbolic-racer run R6/frozen/scoped.jsonl R6/software soft
  python3 racer/certify.py R6/native
  python3 racer/certify.py R6/software
  python3 racer/audit.py R6

Prepare verifies the current combine through every scoped application trace row.
It preserves original args, selected result, expected answer, source and environment
inside each scoped/excluded record. Unknown scoped expected outcomes abort the run
as a failed load. No input, answer or source is synthesized. A missing result is
never a candidate pass. The original loader's unread-row accounting must accompany
the trace; PMT has zero unread rows in this unit.

The six structural templates each have eight arithmetic sites, with every one of
3^8 strict/stored-PC64/resident-PC64 assignments: exactly 39,366 serialized graphs.
Candidate hashes accelerate lookup; full graph equality checks duplicates. The
three combines are the explicit association choices in the unit. Only 0/1 constants
are admitted, including canonical payment timing. log1p and spilled expm1 leaves
remain frozen. Strict consumers store extended operands before arithmetic. Extended
values retain all ten bytes until a store. Final publication follows Site::number,
including its zero, subnormal refusal and nonfinite policies.

First-kills.jsonl records every candidate's first contradiction with expected and
calculated outcomes; programs.jsonl supplies each complete program. Witness-kills.jsonl
includes every scoped witness, including zero kills. The first 734 rows are the
previous quotient gaps; the first 32 are discovery, all previously inspected.
Elimination checkpoints retain next witness order and surviving candidate IDs.
Completed reachability records carry stable order/trace IDs; absent records mean
incomplete coverage, not a successful search. This initial runner can restart a
run; it does not yet resume an interrupted reachability worker.

Reachability uses every candidate, including previously killed ones. Exact tagged
value sets share evaluation at subtrees, without claiming global program equivalence.
There are no computed shared subtrees in these templates, so the Cartesian joins
are independent. The self-test verifies full graph enumeration against these sets
on 32 initial gaps plus 32 evenly spaced scope rows, on both backends. It also checks
360,000 arithmetic/continuation parity cases and three actual invoke-route publication
controls. Synthetic arithmetic controls are not Excel observations.

Unreachable-root records retain every possible final binary64/Ext80 value. certify.py
decodes each value as an exact rational and proves it lies outside the recorded
target's RN53 rounding interval, including tie endpoint parity. No approximate
inverse is used. These certificates constrain this finite family and frozen leaves;
they do not exclude other Excel programs or identify its internal helper precision.
The native and software interval certificates and all raw roots agree exactly.

R6's result is zero survivors, with 38,432 of 250,730 scoped targets unreachable.
An individually reachable answer is not an Excel pass or a selectable row-specific
program. No product kernel or recorded Excel answer changes. There are no frozen
questions because there is no surviving program set to distinguish. The next unit
should search the expm1 precision/store boundary, using the three FV=0 killers as
constraints and keeping the one-period identity as a separate mathematical control.

Complete external artifacts, checksums and work report:
  /home/ubuntu/codex-runs/dna-firehorse-20261031-round6/
  /home/ubuntu/codex-runs/work-report-6.md
Compact outcome/provenance is retained in pmt-r6-result.json and on dc-j9k.

Round 7 extends the same evaluator and exact interval code. PMT's final bounded
family has1980 helpers, all60 four-leaf multiplication/division trees and81
three-model combines per tree:9622800 programs, zero survivors. The independent
shared-helper projection rejects all4860 combine profiles even with nonunit
periods and timing fixed. The full cross-matrix attempts are retained as partial
diagnostics; the complete first-kill group proof is a separate result.

Poisson work starts from the defining PMF and uses typed count-one ASTs, exact
constant-load certificates, scalar/Gamma forms, reciprocal logs, normalization
loads, underflow guards and five arithmetic/store models. The largest mixed-build
family has26043360 programs and no survivor. Build20527's88 tiny-mean readings
alone leave52908 programs in the1585920-program refinement; this is a scoped
hypothesis set and does not identify a full Excel program or confirm build drift.
Questions use the existing xtask witness read --inputs schema; predictions stay
separate and are never ingested as answers. No application kernel changes.

From this directory, with R7 as the outside-repository round7 artifact directory:
  target/release/symbolic-racer boundary-wide R7/frozen/scoped.jsonl OUT hard
  python3 racer/group_eliminate.py R7/frozen/scoped.jsonl OUT PROOF
  target/release/symbolic-racer precision-test
  target/release/symbolic-racer mass R7/poisson-frozen/scoped.jsonl OUT hard
  target/release/symbolic-racer count-one-full R7/poisson-frozen/scoped.jsonl poisson-r7-wide-constants.json OUT hard
  target/release/symbolic-racer count-one-full-test R7/poisson-frozen/scoped.jsonl poisson-r7-wide-constants.json
  target/release/symbolic-racer count-one-builds R7/poisson-frozen/scoped.jsonl poisson-r7-wide-constants.json OUT hard
  target/release/symbolic-racer normalizer-matrix R7/poisson-frozen/scoped.jsonl poisson-r7-wide-constants.json OUT hard
  python3 racer/normalizer_intersect.py OUT PROOF

Use soft instead of hard for software replay. DNA_RESEARCH_WORKERS controls1..8
reachability workers. Exact per-input/backend tagged-request memoization preserves
the physical program description; uncached reference checks audit that optimization.
The numerator projection never selects its unknown denominator as a constant.
Canonical unit descriptions distinguish original frozen metadata from later prose
corrections, and name the actual three-constant or wide-TAU certificate used.
Full artifacts/report: /home/ubuntu/codex-runs/dna-firehorse-20261031-round7/
and /home/ubuntu/codex-runs/work-report-7.md.

Round7 publication correction: Poisson aliases use out_with_subnormal and preserve
finite subnormals. PMT keeps its distinct refusal policy. Initial Poisson output
sets and question-draft are superseded; final runs are poisson-final-*, mass-final-*
and builds-final-* in the artifact directory. Run certify.py DIR poisson for those
runs. Four actual invoke alias/false-flag controls exercise numeric subnormals.
Final questions add a count-one mean720 domain-control pair and use distinct final
IDs; requests contain no answers. Program predictions are separate files.

Final upstream master e663103 adds WorksheetErrorCode::excel_spelling, an unused
additive method in the guarded calc-value source directory. A fresh research build
correctly refuses that delta; the frozen-source guard was retained. The numerical
executables and question generator were built at ec05eae before that pull. Rebuild
this unit against the sibling DnaCalc at ec05eaeb3a30a56e492361d1c76fa17d0e649eaa
(or the original frozen revision with identical guarded source bytes). No numeric
source changed in the final upstream pull. The complete compiled binaries and
source hashes are in the external artifact directory.

Round 8 first scores the immutable Round7 predictions against the 38 new Excel
answers. Every one of 52,908 hypotheses has a separate exact pass/fail row in
R8/poisson-blind/per-hypothesis-blind-scores.csv. None scores 38/38; the highest
score is 30/38 (196 hypotheses). Only after this blind result is frozen does exact
elimination reject all 52,908, at the mean 1e-6 boolean-FALSE question. The corresponding
numeric-zero flag gives the same recorded answer. No Poisson program is promoted.

Annuity structure discovery precedes the new searches. It retains exact typed
rate/nper/environment groups, conditioned kind-transition brackets, the original
38,432/11,932 gap identities, and isolated growth-factor preimages. NPER/RATE
output-derived coordinates are separate projections; no internal RATE solver
trace is inferred from externally evaluated FV observations.

On 235 isolated FV growth observations, only the named binary-squaring diagnostics
remain compatible with every applicable exact preimage. This motivates joint
PV/FV programs rather than a new PMT combine ranking. The 354,294-program compact
integer family has 7,164 survivors on 1,321 original rows. Its 13-provider fractional
refinement has 1,824 survivors on all 3,785 recorded positive-base/nonzero-rate
PV/FV rows (including negative and fractional periods). Zero-rate/nonpositive-base
controls remain outside that grammar. Larger separately factored families constrain
the same four integer-growth IDs but remain broader hypotheses, not confirmations.

The immutable 30-question packet is at:
  /home/ubuntu/codex-runs/fresh-questions-r8/
It has 16 PV and 14 FV requests, 26 new typed identities and 4 intentional rechecks.
Requests contain no answers; native/software predictions and all 1,824 physical
program descriptions are separate. Question source is fc14027112ea72444b43c720218e05503a07d555.
The matching executable is retained outside this repository. Use on Windows:
  cargo xtask witness read --func PV --inputs requests-PV.json
  cargo xtask witness read --func FV --inputs requests-FV.json
Do not pass --new-only. Score these frozen predictions before any later elimination.

Three PMT-specific followups have zero survivors: positive shared PV/FV growth
(5,773,680 components), five independently specified log producers with coupled
growth/expm1 (36,085,500), and the shared interest/amortization decomposition
constrained by NPER (1,603,800). Exact per-build reverse projections also test
arbitrary common helpers, alternative numerator/division/multiplication forms,
binary64 FMA and five arithmetic models. Their completed sufficient-row certificates
are reviewed independently. These exclude specified structures, never all annuity
programs. The hardware FYL2XP1 leaf is measured natively in BOTH replays; the second
run is hybrid, not a software implementation of that instruction. The log series
is DLMF 4.6.E6, summed with PC64 until its value stops changing; no correct-rounding
claim is made. Sources: https://dlmf.nist.gov/4.6.E6 and Intel's public SDM.

Input preservation correction: an early research-only decimal f64 JSON path changed
150 arguments in 142 integer cases. All such joint runs are invalidated; only the
bits-* reruns and later u64-backed runs count. Excel observations were not changed.
Resident growth preimages also now include final RN53 publication; only corrected
growth-preimage/comparison directories count. The frozen packet uses both fixes.
The independent PC53 audit found a software-backend underflow discrepancy at the
extended exponent boundary; no global native/software parity claim follows from
ordinary-range controls or the matching recorded-input exclusions.

Use R8=/home/ubuntu/codex-runs/dna-firehorse-20261031-round8. Build this round with:
  DNA_RESEARCH_UNIT=pmt-prefix-r8-unit.json CARGO_TARGET_DIR=target cargo build --release --offline -p symbolic-racer
Source units, complete grammars, commands/resource logs, exact first kills,
survivors, supersession notes and independent reviews are retained there.
The compact final index is round8-result.json; the report is:
  /home/ubuntu/codex-runs/work-report-8.md
No application kernel, evidence pin or recorded Excel answer changes.

Later R8 structure checks keep the whole-function constraints separate. A single
shared annuity factor with one PV multiply and one PMT divide is rejected by all
25 five-model profiles on paired build20131 observations; this does not exclude
more complex consumers of a shared helper. Direct F2XM1 expm1 with coupled growth
rejects all 5,314,410 PMT programs on 220,888 bounded rows. Its raw native/software
components sometimes differ even though the first kills agree.

The one-loan backward balance check treats internal PMT as every one of the 2,049
Ext80 values publishing to its observed value, with the same unknown PMT across
seven periods. All 899,640 components / 2,698,920 partial programs reject in IPMT;
PPMT is not reached. Transferring the frozen PV/FV growth into that same balance
rejects 504 components / 1,378,944 joint partial assignments. This leaves the 1,824
standalone PV/FV hypotheses intact and does not fit or install a latent payment.

Five-model arbitrary-helper proofs reject the original 37,500 four-leaf profiles
separately on both builds for all/nonunit/timing cohorts. The stronger nonunit plus
timing build20131 proof also excludes all 37,500: an original singleton group
rejects the remaining 312 after the retained 37,188-profile partial passes.
The stopped sequential run and
the independent sufficient-subset closure have distinct provenance and order.

RATE discovery has 231 mismatching cases whose recorded external FV companions
all match. Only one of 60 one-period RATE answers is the RN53 exact linear root;
that is a mathematical diagnostic and supplies no internal stopping tolerance.
All 625 five-model one-step Newton consumers of published FV(guess) are rejected.
Transferring the actual PV/FV internal BG/paid/sum into four residual forms rejects
22,800,000 logical assignments on both backends; inactive site choices are retained
as duplicates. All 109,440 original FV publications are checked per backend across
the 60 cases. This rejects bounded one-step consumers, not all RATE algorithms,
and leaves the standalone PV/FV hypotheses intact.

The final upstream8311ace pull changes udf.rs in the guarded source tree. The old
R8 unit correctly refuses a fresh build there. annuity-r8-rebuild-unit.json pins
that revision explicitly for rebuild validation; all eight numerical source hashes
remain unchanged. Original units, run executables and the question generator retain
their frozen revisions and hashes. No measurement or packet is regenerated.

The late NPER-constrained I*G family forms I=PV*rate before consuming G. Its
48,114,000 difference/correction and 7,085,880 direct-F2XM1 profiles both reject
on the same 220,888 rows, on native/hybrid backends. This is a separate earlier
structure from PV*G first or I+I/em. The present-value RATE residual also rejects
22,800,000 bounded one-step logical assignments; original PV-mode metadata has an
explicit active-field erratum, and later source labels select the active mode.
These late followups pin e9cb285 in their source units; earlier units and the
separate8311ace rebuild validation keep their own original source boundaries.

Round 9 adopts the Pro review's first milestone. checker/exact.py is an independent
integer-coordinate/rational nearest-even implementation, with no racer or Fpu
imports. It checks inverse-rounding cells, open/tie endpoints, lattice absence,
original observation identities and full-width arithmetic traces. Checker receipts
freeze source hashes at startup and refuse success if the source changes mid-run.
Run its deliberate controls with:
  python3 -m unittest discover -s checker -v

The original 734 quotient gaps are reproduced over 250,730 rows, with 249,996 constructive
reachable quotients and every original quotient/product calculation checked. The
54-bit result is lattice reachability only. 17e42ce changed the growth arm's initial
quotient store; the gaps concern the discount final shell RN53(RN64(q*r)). They
refute that binary64-quotient shell universally, not the separate growth hypothesis.
The independently checked 4,860 shared-helper combines and the 37,500 five-model
profiles per build have no survivors. Even arbitrary common helpers cannot rescue
those declared families; neither the frozen helper nor Excel's source is identified.

The one widened axis is ordinary association: 4 numerators x 8 combines x 3^8 plus
the shorter (FV+(PV+FV)/em)*(r/tf) with 3^7 assignments, 212,139 programs in total.
Both backends reject every program within 8 existing rows. Programs and every
full-width first-kill trace are byte-identical between backends; the separate
checker verifies 2,967,759 node values. No programs survive to justify questions.
The concrete function claim and assumption audit are pmt-r9-claim.json; the frozen
unit is pmt-r9-unit.json and the result index is round9-result.json.

Replace R3/R7 below with the original outside-repository artifact directories and
R9 with a new output directory to reproduce the quotient and grouped checks:
  python3 checker/check.py quotient R3/logs/dna-round3-pmt-quotient-lattice-raw.tsv R3/pmt-quotient-gap-certificates.jsonl R9/quotient.json
  python3 checker/check.py groups R7/frozen/scoped.jsonl R7/wide-native R9/groups.json R7/shared-group-wide-first-kill/all-first-kills.jsonl R7/shared-group-wide-first-kill/nonunit_timing-first-kills.jsonl --count 4860

Build the frozen unit, then use the same evaluator and witness loader as before:
  DNA_RESEARCH_UNIT=pmt-r9-unit.json cargo build --release --offline -p symbolic-racer
  target/release/symbolic-racer arithmetic-reference-controls R9/controls.jsonl
  python3 checker/check.py controls R9/controls.jsonl R9/controls-summary.json
  target/release/symbolic-racer associations R7/frozen/scoped.jsonl R9/native hard
  target/release/symbolic-racer associations R7/frozen/scoped.jsonl R9/software soft
  python3 checker/check.py graphs R7/frozen/scoped.jsonl R9/native/programs.jsonl R9/native/first-kills.jsonl R9/checked-graphs.json

Format only research files, using rustfmt with skip_children=true; cargo fmt --all
traverses the path dependencies into DnaCalc. The source guard rejected that
unintended format diff; it was backed up and reversed before any numerical run.
All receipts, supersession notes, source/executable freezes and review are under
/home/ubuntu/codex-runs/dna-firehorse-20261101-round9/. No product code, recorded
Excel answer, evidence pin or application pass changed. PMT remains unidentified.

The later upstream UI split at 065687a adds unrelated clipboard/axis helpers in
calc-value/src/lib.rs. The original unit's source guard remains intact; select
DNA_RESEARCH_UNIT=pmt-r9-rebuild-unit.json to rebuild against that source instead.
All eight numeric source hashes are unchanged. Both rebuilt backend runs have
byte-identical programs, full first-kill traces and survivors; the 70 independent
reference comparisons pass. Original measurement units/binaries stay unchanged;
this rebuild is not fresh Excel validation and supersedes no recorded result.

Round 10 first freezes 64 PMT observation questions: 24 certified gaps and eight
passing controls re-read numerically, plus 16 literal, eight logical-timing and
eight text-rate audits. Requests have no answers; historical/application
predictions remain separate. The reader records requested/submitted/read-back
inputs and settings. Literal observations stay in counted transport audits
outside semantic banks. No build 20527 answers exist yet. Use the two commands in
/home/ubuntu/codex-runs/fresh-questions-r10/RUN-ROUND10.txt on Windows.

The next function is NORM.S.DIST CDF: 13,716 failures in 21,810 calls, one scalar
argument and published rational regions. GAMMA.DIST CDF has 50,562 failures but
also shape, scale and iterative structure. Published Cody/SPECFUN and BSD3
SciPy CDFLIB constants supply fixed sketches; the separate Hart body uses only
the shorter West 2004 Figure 2 table. Full-Hart higher precision tables remain
unsearched. No coefficient or helper value is fitted.

The original 18,000 Cody programs and its 72,000 reconstruction extension have
zero survivors separately on builds 20131 and 20430; the original is a subset,
not an extra distinct family count. The 7,200 Hart-body profiles and 18,000
central projections also reject. Central checks use no exponential leaves.
Eight final independent integer/rational checks verify 15,423,150 full-width
trace node values. Concrete exponential leaves remain conditional, never identified Excel
helpers. Native/software traces agree for the original and extended Cody families.

All 275 declared tiny-linear profiles reject even with an arbitrary common Ext80
coefficient per observation group. Of 660 sign/reconstruction relations, 655
reject and five remain feasible on 167 old sign pairs. These are necessary
constraints with an arbitrary latent value per pair, not complete numerical
programs. The 24 frozen normal questions include eight new positive companions
to historical negative anchors, and four new mathematical lattice sign pairs.
Exact conditional outcome sets and a fixed three-case probability interval table
are recorded separately from requests. MPFR supplies mathematical proposals,
never Excel answers or private helpers.

The claim, bounds and source audit are normal-cdf-r10-claim.json; receipts and
source hashes are indexed by round10-result.json. Superseded failed attempts
remain in the external artifact root, including an invalid add inverse rejected
by the checker and corrected with the existing inverse_add implementation.
No numeric kernel, evidence pin or recorded Excel answer changes. All 27,473
application trace rows are byte-identical before/after: 13,757 pass, 13,716 fail.

For reproduction, set R10 to the outside-repository round10 artifact directory.
Build against the guarded source revision in the selected unit, then use:
  DNA_RESEARCH_UNIT=normal-cdf-central-r10-unit.json cargo build --release --offline -p symbolic-racer
  target/release/symbolic-racer normal-cdf R10/application-NORM.S.DIST-before.jsonl normal-cdf-r10-constants.json R10/new-original 20131 hard
  python3 checker/normal.py R10/new-original normal-cdf-r10-constants.json R10/application-NORM.S.DIST-before.jsonl R10/new-original-checked.json
  python3 checker/normal.py linear R10/application-NORM.S.DIST-before.jsonl R10/normal-linear-projection-corrected R10/new-linear-checked.json
  python3 checker/normal.py symmetry R10/application-NORM.S.DIST-before.jsonl R10/normal-symmetry-complete R10/new-sign-checked.json
  python3 checker/normal.py questions R10/application-NORM.S.DIST-before.jsonl /home/ubuntu/codex-runs/fresh-questions-r10 R10/new-questions-checked.json
  python3 -m unittest discover -s checker -v

Primary numerical sources:
  https://netlib.org/specfun/erf
  https://raw.githubusercontent.com/scipy/scipy/v1.14.1/scipy/special/cdflib.c
  https://s2.smu.edu/~aleskovs/emis/sqc2/accuratecumnorm.pdf
SciPy's BSD3 notice is retained in notices/SciPy-cdflib.txt. The source audit
records incidental nonadmitted search previews that were not used or saved.

Round 11 prioritizes the workbook cells. The census accounts for all 1,184
baseline mismatches; its formulas/XML and classifications are retained in
/home/ubuntu/codex-runs/corpus-mismatch-r11.md and the round 11 artifact directory.
Original N, FILTERXML and exact-proportional-regression fixes in application
commit 4326fc06 (published through e8a15f90) give 59 corpus gains and zero losses.
All 3,874,816 decoded bank outputs remain byte-identical; 116 unread rows still fail.

The Round 10 answers in evidence 51bed087 were scored against frozen predictions
before the new search. PMT numeric rechecks reproduce 32/32 historical outcomes,
including 24 certified gaps; application predictions match 8/32. All 8 logical controls
reproduce their paired numeric answers. Fourteen of 16 literal calls change echo
inputs and remain failed transport audits. All 8 normal anchors reproduce, and
two of the five necessary sign relations pass all 12 pairs; no complete CDF follows.
The reviewed blind-v2 score and subsequent elimination are separate artifacts.

CUMPRINC has 157 corpus mismatches. Its 540 bank observations provide 60 loans, each
with 9 intervals and an independently recorded PMT companion on the same build.
The 90,000 logical integer-growth/payment/recurrence assignments all fail the first
singleton. There are 18,000 instantiated graphs, with 1,509,000 independently checked
node values; this supplies no evidence for an accumulation order. Both backends
have identical complete traces, kills and growth tables. Two separate 50-profile
singleton projections exclude their specified growth consumers even with an
arbitrary shared principal or an internal payment constrained by external PMT.
No coefficient is fitted. Numerical research changes no annuity or approximation
body or saved answer; the separate structural fixes above change the product.

Use cumulative-r11-unit.json for the guarded build. The claim, original bounds,
receipt/source hashes and remaining assumptions are cumulative-r11-claim.json;
round11-result.json indexes the measured outcomes. The 64 frozen coupled questions
are versioned in packets/round11 and outside the repository in fresh-questions-r11.
Requests and predictions remain separate. The exposure erratum adds the two
known PMT controls omitted from the frozen selection's partial exposure listing.
The external packet and its original generator snapshot remain immutable.

From the DnaCalc checkout on Windows after pulling both repositories:
  $packet = (Resolve-Path '..\DnaCalc-Research\packets\round11').Path
  cargo xtask witness read --func CUMPRINC --inputs "$packet\requests-CUMPRINC.json"
  cargo xtask witness read --func PMT --inputs "$packet\requests-PMT.json"
  cargo xtask witness read --func PPMT --inputs "$packet\requests-PPMT.json"
  cargo xtask witness read --func CUMIPMT --inputs "$packet\requests-CUMIPMT.json"
Do not use --new-only. Score every immutable prediction first when answers arrive.

The separate admission/text audit uses packets/structural-r11-v2 (52 calls).
Its numeric-only CONCAT/CONCATENATE requests replace the original packet's
unnecessary empty second argument, whose Value2 read-back might be blank.
The original packet is preserved; no returned answer motivated this revision.
Use the v2 directory for BESSELI, BESSELJ, CONCATENATE and CONCAT, with each
function's requests-<function>.json. These are boundary/route audits; their
kind-only hypotheses are never exact numerical function matches.

Upstream engine/value-model work at application16b00e66 changes the guarded
source tree. The numerical source hashes remain unchanged. The original unit
still pins e8a15f90; cumulative-r11-rebuild-unit.json explicitly pins the later
tree for rebuilds, without replacing any measurement, packet or original unit.
