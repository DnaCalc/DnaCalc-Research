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
