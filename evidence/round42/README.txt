Round 42 (Claude Opus, dna-firehorse): NORM.S.INV and LINEST/LOGEST QR

The question service answered no call from 02:13 to 04:36 UTC ("worker-down:
no worker is up"; the workers stopped at 02:16 after excel-timeout restarts).
The orchestrator restarted it at 04:37; every frozen question was then asked
(Excel 16.0 build 20527, env f3269992527c). Analysis before 04:37 used the
answers retained by rounds 38-41.

Landed in DnaCalc 296cedb6 (on dd35297b; light check 697 passed, 0 failed):
- NORM.S.INV/NORM.INV below p - 0.5 = -0.425 is the AS241 tail at p, Horner
  steps and quotient stored from precision 64. Fresh frozen lower-tail set:
  1,999 / 2,000 (old program 392). Retained R41 sweep 81,987 -> 393,140 /
  393,216; R39/R40 1,646 -> 2,715 / 2,716 (1,283 upper answers kept).
- LINEST/LOGEST/TREND/GROWTH: pivot on the running remainder (an exact tie
  goes to the later column with a constant, the earlier without), stop at
  the first pivot that fails the collinearity test, zero diagonal entry
  reflects as negative. Fresh set 3, frozen with that program: LINEST and
  LOGEST 1,303 / 1,344 each (old 903 and 932). Retained R41 LINEST 13,915 ->
  14,852 / 15,150 and LOGEST 10,973 -> 11,849 / 12,140, R38 1,002 -> 1,266
  and 937 -> 1,195 / 1,296, every one with zero losses.
Bank 3,887,837 -> 3,887,839 / 3,910,773 (TREND 10 -> 12/12); corpus 58,451 ->
58,468 / 60,296 (+17, -0), both measured at 15098f52 without and with the
change (measurements/). Corpus failures 1,845 -> 1,828: 521 mismatched,
1,111 refused, 1 missing, 73 external, 14 circular, 108 awaiting capture.

How the rules were reached (design/regression-qr.txt): fresh set 1 refuted
the plain stop with fresh-squares pivoting (D1 143/288) and showed the
running-remainder pivot (288/288); set 2 (four builds frozen) confirmed it
(447/448) and exposed the constant-dependent tie; set 3 proved the combined
program. Still open: 41 set-3 selections, all zero-column layouts beside a
collinear pair.

NORM.S.INV central (design/normal-inverse.txt): sqrt(2)*erfinv(2p - 1); near
0.5 y = RN53(RN64(2q*RN(sqrt(pi)/2))), x = RN53(RN64(y*RN(sqrt 2))) matches
all 16,383 retained answers. Fresh mapping (79 windows) puts the erfinv
piece split at t^2 = 1/2 (|q| = 0.35355); pieces, degrees and coefficients
are not identified; Pade approximants about 0 are refuted. (0.925, 15/16] is
a separate approximation with a flat +0.85 ulp bias, unidentified.

Retained: questions/ (every frozen set, freeze.sha256, predictions of each
build), fresh/ (the original replies), scripts/ (ask.py, score.py, the
analysis scripts, r42-final.patch = the landed diff), measurements/.
files.json maps each payload to its source and both hashes. Waits over 3 s:
none; the 9 lane probes and 8 minute-polls that failed (worker-down, or
no reply while the broker was stopped) before 04:37 are in excel-waits.tsv and poll.log, and were
re-asked after the restart. Sources: Excel answers, IEEE 754, AS241
tables, published Householder QR; no tier 3/4 implementation was read.
