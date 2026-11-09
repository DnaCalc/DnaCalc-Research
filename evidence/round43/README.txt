Round 43 (Claude Opus, dna-firehorse): LINEST/LOGEST QR order and test,
the regularized beta at a = 1 and BINOM.DIST one short of n; NORM.S.INV
central and incomplete-beta mapping

Every question went through the question service (Excel 16.0 build 20527,
env f3269992527c); no request took longer than 1.2 s and none waited over
3 s (fresh/waits.tsv is empty; the direct scripts record their own seconds).

Landed in DnaCalc:
- 6996cf6d (regression_forecast_family.rs): with a constant the ones column
  is stored after the predictors and swapped to the front; pivoting swaps,
  an exact tie goes to the earlier column, a negative running remainder is
  read as zero; while columns are still to come the collinearity test reads
  vv against 2^-41 of the initial sum of squares. Fresh proof frozen before
  asking (questions/proof.sha256): LINEST 2,225 -> 2,413 / 2,439, LOGEST
  2,255 -> 2,414 / 2,439 (25 unread per function, counted as failures).
  Retained cohorts +791 / -0. Bank and corpus unchanged (measurements/
  baseline and after-linest, both at 70ecfa0b).
- 431bcdc3 (special_math_common.rs, discrete_dist.rs): BETA.DIST(x, 1, b) for
  x < 0.1 is -expm1 of the x87 log1p_times chain; BINOM.DIST(n - 1, n, p, TRUE)
  for 1 - p < 0.1 is -expm1(n * log1p(-q)) with a binary64 product. Fresh
  proofs frozen before asking: BETA.DIST 1,997 / 2,000 (430), the doubles
  below 0.1 100 / 100, BETADIST 298 / 300; BINOM 998 / 1,000 (844), RANGE
  200 / 200. Bank 3,887,839 -> 3,887,866 / 3,910,773 (BETA.DIST +11, FDIST +8,
  NEGBINOM.DIST +6, F.TEST +2); corpus 75,042 / 76,957 unchanged
  (measurements/head2 and after-beta, both at 6996cf6d).

How: design/regression-qr.txt section 5 (x87 emulation scripts/qr2.py on
2,807 retained degrees of freedom; threshold maps questions/LINEST.thr*),
design/incomplete-beta.txt (shape map, error scale, the a = 1 map and leaf
ranking, the BINOM-beside-BETA map).

NORM.S.INV (design/normal-inverse.txt section 5): split at t^2 = 1/2 to
within four doubles (fresh/NORM.S.INV.split), piece C accurate to one ulp,
piece A a cancelling form (lattice map fresh/NORM.S.INV.lattice); every
published erfinv candidate in the brief refuted; no kernel changed.

files.json lists every retained payload with its size and SHA-256. Sources:
Excel answers, IEEE 754, AS241 tables, published Householder QR, the TOMS
708 bratio already in calc-func, mpmath; no tier 3/4 implementation was read.
