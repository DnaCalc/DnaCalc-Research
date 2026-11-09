R48: exact alternating ETS diagnostics; broader search trees refuted

Original MIT research by GPT-6.1 Sol (Codex), dna-firehorse, implementer.
Sources are original Excel20527 question-service replies, IEEE754, published
Hyndman ETS mathematics, the Microsoft Power View2014 article, and cleared
MIT DnaCalc/Research code. References and ranked trees: design/ets-search.txt.
No tier3/4 implementation source or excluded checkout was read or changed.

Landed scope: exact zero-innovation alternating affine histories, length5..32,
including unfinished final cycles, now admit period detection, reported ETS
statistics and zero confidence widths through the existing initializer and
shared loop. Existing exact lattice/magnitude bounds remain. Forecasts in the
new region stay refused; the previous complete-cycle forecast region remains.
The original23-point parameter search and score-only loads are unchanged.

The initial wider trial fc7a0907 included forecasts. Reviewer-selected
DC/trend/amplitude lattice-edge histories expose four near-cancellation forecasts
where Excel returnszero and that trial returns+/-2^-400. Commit10798f3d fixes
forward by withdrawing the entire newly extended forecast consumer. No
near-zero threshold, horizon exception or arithmetic correction is fitted.
All183 original exact call expressions remain. One detector refusal becomes
an exact period2 pin; original six-point automatic forecast refusals remain.
Original trial sources/answers/commits are preserved, including the contrary
answers. The first title's accidental identifier spacing is corrected by
canonical references in the follow-up commit, without rewriting history.

Original corpus stays75209/76957:1748 inclusive failures, including952
refused. ETS refusals before->after: FORECAST256->256, CONFINT256->256,
SEASONALITY368->368, STAT72->72. First-refusal blocks remain automatic528,
short initialization204, seasonal innovations180, cycle limits40. These
original refused corpus cells are outside the newly proved diagnostic set.
The isolated ETS change gainszero original corpus or bank matches.

Incoming cumulative-finance work9871f8b1 changes three still-failing values;
it changes no match/refusal count. Incoming SWITCH3e8d2645 adds two bank
matches:3887881->3887883/3910773, with22890 remaining failures, including
4 unread and2064 excluded. Incoming J4 fixtures95c1b3bf add604 formula cells:
579 match,14 mismatch,11 VDB refusals. Expanded corpus75788/77561 has1773
inclusive failures:614 mismatches,963 total refusals,1 missing,73 external,
14 circular and108 pending. ETS still supplies952 of those963 refusals.
The added coverage is separate from product gains. J4 changes workbook
identity and session numbering, so the strict original-context audit rejects
it. j4-shared-readings.json instead matches sections on every capture field
except source, then verifies201691 original evaluation payloads and76849
observed verdicts unchanged;108 pending cells remain in the denominator.
That supplementary result does not assert identical workbook context.
The abandoned path-only j4-original-scope-mapping belongs to the rejected
strict comparison; final section identities are in j4-shared-readings.json.

Search:36 complete bounded-grid/refinement/simplex/score-window trees on96
histories give74 grid/simplex,10 error-model and59 window separators. All35
alternatives fail. The sparse-grid survivor already fails R47 fresh bits;
no search change lands. MSE and profiled additive likelihood are monotone on
the same fixed residual window; no real-minimizer separating history exists.
Detector:all30 simplified spectrum/ACF/error-minimizing trees fail220 fresh
histories. A fixed exact larger-period extension fails108/288. No cutoff,
peak correction or extra arithmetic piece follows these counterexamples.

Short initialization:all3072 joined raw-moment and2048 centered-covariance
programs fail318 retained conditional histories. Original trials mistakenly
coupled PC53 initialization to PC53 final publication. Review found this;
V2 fixes publication independently toPC64, withzero survivors again. Both
epochs remain. Reported reset parameters do not certify private raw beta.
Seasonal innovations:all12 conditional initializer trees fail96 fresh
positive-trend histories; existing forward/raw initialization fails57/384.
Cycle limits:128 four-point period0/2 pairs give3840 equal scalar readings,
which identifies routing but not the missing short arithmetic. Those regions
stop at these counterexamples and retain their typed failures.

Final selected route proofs:4352/5376 discovery,13032/15288 independent,
900/1260 lattice-edge match;1024/2256/360 forecast refusals respectively,
zero numeric discrepancies. Every diagnostic selection matches. The shared
152376-selection registered cohort gains8268 exact diagnostics from refusals
(STAT5770, CONFINT1640, SEASONALITY858), withzero lost matches or newly
calculated discrepancies; every shared forecast value remains unchanged.
Final mixed replay:108078/153925 matched,44557 refused,1290 mismatched.
These include discovery and conditional questions, not153925 fresh proofs.
All refusals and discrepancies remain failures, with the full denominator.

Canonical final bindings: closing-measurement-bindings.json and route-latest.
Final guarded Research unit is round48-final-unit.json at application6f0a01fb;
J4 changes no calculation sources and retains the exact calculation executable.
The isolated final-measurement-bindings and combined integrated/latest files
remain useful earlier source epochs. Earlier route-final/measurement files
refer to the wider trial; narrowed/published/latest supersede that scope.
source-epochs-corrected.json supersedes an earlier prose qualification of the
baseline executable; original compiled provenance was always retained.
The initial inherited /tmp target ran out of space; subsequent builds unset
both inherited target variables and use each clone's normal target.

Validation: final main light check749 passed,0 failed,3 ignored,1 not run;
Research154 Rust passed,2 ignored,165 Python passed. UBS haszero critical
findings;157 ETS Rust and67 Python warnings are reviewed in retained receipts.
The incoming keeper source scan adds57 reviewed warnings andzero critical.
No Windows/browser sweep or CI dispatch. The existing audit handles an empty
slow-wait list withzero rather than raising; all reply validation remains.
Service:288 unique audited replies,27137 numeric selections,zero service
refusals or unqualified answers. All use interactive priority and pass the
existing input/frame/provenance checks. No wait over3s; maximum0.462061s.

Storage: r2-complete-receipts.json records all checksum copies and one-way
checks under dnacalc-evidence/research/round48:5 files/85288911 bytes,
8511 verified archive members across four fresh dated prefixes. Bulk corpus
and registered replay payloads stay on R2. Small complete question replies,
frozen predictions/programs, counterexamples, epochs and summaries stay here.
files.json uses the single retain_evidence helper; every original/compressed
SHA256 is checked, and previously retained bytes are immutable. Immutable
historical draft commit/comment text is not proof of what was published;
Git history and the source/publication receipts identify the actual actions.
An incoming9b13be30 keeper tracker update changes only tools/xtask/src/keep.rs;
its merge is checked with both calc-func and xtask tests, then published
with both requested trailers. final-publication-bindings.json identifies it.
Git-created merge9fd77961 retains its automatic message without the requested
trailers; published history is preserved. Authored work commits carry both.

Reproduction uses the existing racer/ets_proof.py registered replay, the
existing audited question client and retained frozen sources/inputs. The
active bounded designs are racer/ets_round48.py; use a new OUT directory for
live questions because original output creation is exclusive. The existing
racer/corpus_cases.py problem-functions command regenerates the sole table
through cargo xtask corpus and the full witness judge. Research unit source
guards require the recorded application sources when rebuilding old epochs.

Read-only review is complete with no remaining material concern. Its receipt
is retained. Final index: 543 immutable payloads/6977724 stored bytes,
all verified by the existing retainer after adding the review and commit body.
