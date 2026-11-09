Round 42 (Claude Opus, dna-firehorse): NORM.S.INV and LINEST/LOGEST QR

The Windows question service answered nothing this round: from the first
request at 02:13 UTC every call returned "worker-down: no worker is up" (all
three workers stopped at 02:16 after excel-timeout restarts; the broker kept
running). No Excel answer was taken, so no fresh proof exists and no product
code landed. Everything below uses answers retained by rounds 38-41, Excel
16.0 build 20527. The service was not restarted from this lane.

Measurements (DnaCalc c3646c36, clean): bank 3,887,837 / 3,910,773; corpus
58,451 / 60,296 (1,845 failures: 538 mismatched, 1,111 refused, 1 missing,
73 external, 14 circular, 108 awaiting capture). With the three candidate
changes (retained/scripts/r42-candidates.patch): bank 3,887,839 (+2, TREND
10 -> 12/12, no other function changes); corpus 58,454 (+3, -0: LOGINV x2 and
LOGNORM.INV in the catalogue's prefix_metadata sheet).

NORM.S.INV (design/normal-inverse.txt):
- Below p - 0.5 = -0.425 Excel publishes the AS241 PPND16 tail at p with
  every Horner step and the quotient stored from precision 64. R41 dense
  sweep: 81,987 -> 393,140 / 393,216 (+311,160, -7); R39/R40 retained answers
  1,646 -> 2,715 / 2,716, the 1,283 upper answers above 15/16 unchanged.
- The central branch is sqrt(2)*erfinv(2p - 1): near 0.5 the program
  y = RN53(RN64(2q*RN(sqrt(pi)/2))), x = RN53(RN64(y*RN(sqrt 2))) matches all
  16,383 retained answers; it is a function of the stored q (0 of 23,550
  sharing pairs differ). Its noise is that of a monomial rational in t^2 with
  large alternating coefficients, with a piece break between |t| = 0.673 and
  0.720. Pieces, degrees and coefficients are not identified; Pade
  approximants about 0 (m + n <= 25, binary64 or x87-stored) are refuted.
- (0.925, 15/16] is a separate approximation: not the AS241 tail (Excel
  minus it -1..+5 ulp), not Newton-refined, not the tail scaled by 1/sqrt 2.

LINEST/LOGEST (design/regression-qr.txt):
- Rank-revealing stop: when the largest remaining column fails the
  collinearity test, every unstepped column is left out with it.
- A diagonal entry that is exactly zero reflects as a negative one.
- Both against c3646c36: R41 LINEST 13,915 -> 14,864 / 15,150 (+1,057, -108),
  LOGEST 10,973 -> 11,861 / 12,140 (+996, -108); R38 designs (not used to find
  the rules) LINEST 1,002 -> 1,239 / 1,296 (+276, -39), LOGEST 937 -> 1,168
  (+270, -39). All 147 losses have a zero column beside a collinear pair (123
  with a constant); that ordering rule is still open.

Frozen questions for the next round (retained/questions, hashes in
freeze.sha256, predictions of both builds beside each): NORM.S.INV lower tail
2,050 points (1,607 separate the builds); LINEST/LOGEST 754 selections each
(372 separate stop from continue); 576 + 960 LINEST selections for the zero
diagonal rule (3 + 12 separate); 79 NORM.S.INV mapping windows (20,224
points: the central split, (0.925, 15/16], above 15/16). retained/scripts/
ask.py asks them all and score.py compares the replies with both builds.

Validation of the candidates (not committed to DnaCalc): calc-func unit tests
pass, including three new ones pinned to retained Excel readings. The light
check on 15098f52 plus the patch passes 693/0 (3 ignored, 1 not run) after
`cargo xtask witness refresh --func TREND` re-pinned the two TREND gains.
Every Excel failure is in retained/excel-waits.tsv (each a 0.2 s refusal; no
wait over 3 s). files.json maps every compact payload to its source and both
hashes. Sources: retained Excel answers, IEEE 754, AS241 tables, published
Householder QR; no tier 3/4 implementation was read.
