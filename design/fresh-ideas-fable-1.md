# Fresh ideas for the stuck cases (fable-1)

Claude Fable 5.1, Windows (GOVERT-OMEN), 2026-11-09. Implementer role; sources: the design notes in this
folder, reports 40-48 and ld1, the review answer on symbolic search, IEEE 754, public mathematics, and
fresh answers from Excel 16.0.20527 (env f3269992527c) through the question service only. No tier 3/4
source. Every number below that is new is in `evidence/fable-1/` (scripts `e*.py`, replies `e*.json`).

## Premises I doubt across the board

1. **"The stuck part is the function."** It is usually a small core behind cheap structural facts that
   nobody nailed because the rounds went straight to coefficient tables. GAMMALN below 8 is the example:
   two hours of structural questions shrank the unknown from [0, 8) to [0.7, 3) (section 1).
2. **"Store masks are the search space."** Excel here is 64-bit (`C:\Program Files\...\EXCEL.EXE`), and MSVC
   x64 never emits x87 for `double`. Kernels that show PC64 stores (AS241 tail, Stirling, QR, PMT) are
   therefore either assembly, or C with 80-bit `long double` under a compiler that supports it on x64
   (Intel). Either way the store points follow *source-level variables*, not arbitrary node masks, and a
   typed constant is RN64 of its decimal text, not binary64. Enumerate three store policies (every op,
   per assignment, none) over plausible source lines, and add "RN64 of the decimal literal" and "RN64 of
   the exact rational" to every constant's loading hypotheses (the plateau tool only steps binary64
   neighbours).
3. **"Match counts kill."** Qualitative regime signatures kill faster and cannot be fitted: which root a
   solver lands on, whether an error vanishes at a root, whether a mirror identity holds bit for bit,
   whether a decision threshold sits at 0. Sections 3 and 4 run on those.
4. **"Excel's other functions are only consumers."** They are free oracles for intermediates: GAMMALN vs
   LN(GAMMA), NORM.INV vs NORM.S.INV, BETA pdf vs cdf, ROUND vs ROUNDDOWN (all below). Ask the relative
   before the absolute.
5. **"A question is one input."** The most informative question is a bisection on a predicate between two
   public functions (`bisect_edge` exists): the GAMMALN 0.7 join took 60 questions.

## 1. GAMMALN below 8

**Rounds assumed** Cody-shaped bands with unknown tables, or shift-up-to-8 programs, and looked for
coefficients through spectral stage frequencies.

**New facts (e1*, 2,500 fresh answers).**
- The upward recurrence `GAMMALN(x) = RN(GAMMALN(x+1) - LN(x))` holds for all x < 0.7 and stops at the
  *literal* 0.7: last holds 0x3fe6666666666664, first fails 0x3fe6666666666666 = RN(0.7). Cody's
  0.6796875 is excluded. Someone typed `0.7`.
- **[3, 8) is identified:** `GAMMALN(x) = RN(GAMMALN(x - k) + LN(p))`, k = floor(x) - 2, p the binary64
  product `(x-k)(x-k+1)...(x-1)` accumulated from the *lowest* factor upward, LN the public LN. 600/600
  across five bands (e1b). The same product with the highest factor first fails 104-113/120; a sum of
  logs fails 61-99/120. GAMMA(x) is *not* `RN(GAMMA(x-k) * p)` (71-84/120): the two functions share
  the reduction idea, not the arithmetic.
- The core [0.7, 3) has a **flat absolute error** of about 2-5e-16 (max 1.2e-15), the same at the roots
  1 and 2 as elsewhere (e1c/e1e). So it is not `(x-1)*R(x)` or `(x-2)*S(x)` in binary64 (those vanish at
  the roots), and not Cody's shape. It is also **not `ln` of a stored double**: near the roots the output
  pins exp(GAMMALN) to 2^-57, finer than the binary64 grid, and 74/558 land on the grid against 75.5
  expected by chance (e1d). And it is not LN(GAMMA(x)) (1-5 of 60 per band).

**Doubt:** that the core is a rational at all. A flat absolute error of a few 1e-16 is the fingerprint of
a *difference of O(1..8) terms*, Lanczos style (`-tmp + ln(sqrt(2pi) * ser / x)`), or of `ln` applied to an
extended (not stored) gamma, or of `stirling(x + n) - ln(product)` kept in registers.

**Approaches, each with its probe.**
- (a) *Binade signature instead of coefficients.* For Lanczos, the error sd is set by `ulp(tmp)` with
  `tmp = (x + g + 0.5) - (x + 0.5) ln(x + g + 0.5)`; it jumps by 2 where `tmp` crosses a power of two, at an
  x that depends only on g. Ask 64-double windows at 40 x in [0.7, 3), compute sd(GAMMALN - truth) per
  window, locate factor-2 jumps, solve for g (candidates 5, 7, 607/128, 4.7421875, 9). No coefficient is
  fitted; a g with no jump at the right x is dead.
- (b) *The 1.5 join.* `bisect_edge` on the predicate "window sd changes" or on any identity around 1.5:
  is it literal 1.5, and does something structural change (direction of a recurrence, a different core)?
  The flat-error finding says the pieces do not vanish at the roots, so 1.5 is not Cody's seam.
- (c) *Extended-gamma test.* If the core is `RN(ln64(G_ext))`, GAMMALN and the public GAMMA differ as
  `RN53(G_ext)` vs `RN(exp(RN(ln G_ext)))`; the implied G from e1d lands on the 2^-64 grid by construction,
  but EXP(GAMMALN) vs GAMMA then has a predictable +-1 ulp pattern against the fractional position of
  G_ext within its binary64 interval. e1 already has 400 triples; compute the pattern offline first.
- (d) *Differential leaves.* GAMMA.DIST pdf = exp((a-1) ln x - x - lnGamma(a)): ask it beside GAMMALN(a)
  and LN(x) at tiny x and check whether the public GAMMALN is its leaf. If yes, GAMMA.DIST becomes a
  second window on the same core.
- (e) *The one miss in [3,4) (e1, 54/55)*: check whether it is the double-rounding case RN53(RN64(a+b))
  vs RN53(a+b). One answer tells whether the final sum is PC64-then-store.

## 2. Incomplete beta (BETA.DIST, BINOM.DIST, NEGBINOM.DIST, F.DIST)

**Rounds assumed** TOMS 708 bratio with Excel leaves, then attacked the tiny-x prefactor's exponent,
normalizer and scaling, with GAMMALN as the limiting leaf.

**New facts (e3, e3b, 1,120 answers).**
- At tiny x the cdf is **not one rounding of its own pdf * x / a** (17/160 exact, most beyond 3 ulp), while
  **GAMMA.DIST's cdf nearly is** (109/160 exact, all within 3 ulp). BETA.DIST's cdf and pdf have separate
  prefactor code; the gamma pair shares one. The pdf matches none of three public-leaf exp/GAMMALN
  compositions (16-26/160).
- Both pdf and cdf scale by 2^(a-1) and 2^a under x -> 2x only to 10-1,000 ulp for integer a (e3b): both
  are exponent calculations, confirming R44 for the cdf and extending it to the pdf.
- The a = 1 arm uses expm1 and log1p (R43). TOMS 708's helpers `rexp` and `alnrel` are exactly those;
  this supports the NSWC/TOMS 708 lineage more than any leaf test did.

**Doubts.**
- *That the tolerance is TOMS 708's.* bratio sets `eps = max(machine eps, 1e-15)` and derives every series
  and continued-fraction stop from it (`tol = eps/a` in bpser, `15 eps` in bfrac). A hand port with
  `DBL_EPSILON`, `1e-16` or 2^-63 changes the number of terms in **every** branch by one: a shared cause for
  "every branch differs 29-60%" that no store placement can repair (R43 found exactly that for bpser).
  Probe: choose (a, b, x) in bpser where the series stops after k terms at eps = 1e-15 and after k+1 at
  2.2e-16, with the extra term flipping the last bit of the sum; 200 such cells decide eps (emulate bpser
  in Python with exact rationals; the series is a few lines).
- *That the prefactor must be solved through the cdf.* GAMMA.DIST shares its prefactor between pdf and cdf
  and has one fewer log-gamma term: identify `exp((a-1) ln x - x - lgamma(a) - a ln beta)`'s store and
  leaf pattern there first, then transfer the pattern to the beta's `exp(a ln x + b ln(1-x) - lnB)`.
- *That BINOM's p < 0.5 route is a different route.* R46 showed every contradiction has an inexact 1 - p:
  that is bratio's (x, y) interface being fed both p and the rounded complement. Nothing to find there
  until the body is right; stop asking BINOM beside BETA.

**Approaches.** (a) eps probe above. (b) GAMMA.DIST prefactor map as template. (c) Isolate helpers by
shape: `algdiv(a, b)` (b >= 8, a small) dominates the tiny-x normalizer for (a, b) = (1.5, 20..100):
intersect exact preimages of the normalizer over 50 b per a; if the intersection is a single extended
value per (a, b), compare with algdiv's published expansion evaluated at PC64. (d) T.DIST(t, df) =
I_{df/(df+t^2)}(df/2, 1/2): its argument formation is public arithmetic, so it is a second, independent
caller of the same bratio with half-integer shapes; the mirror F.DIST(x, 1, df) must agree bit for bit if
both go through bratio the same way, and the first disagreement localizes a caller-side store.

## 3. NORM.S.INV central region

**Rounds assumed** `sqrt(2) * erfinv(2p - 1)` with pieces in s = t^2 and searched rational forms, Pade
tables and cancellation suffixes.

**New facts (e2, 2,800 answers).**
- `NORM.S.INV(1 - p) = -NORM.S.INV(p)` **bit for bit** in every region except p in (0.925, 15/16], where the
  lower mirror (1/16, 0.075) is the AS241 tail and the upper side is 1-6 ulp larger in magnitude (24/300
  equal). So the code tests `p < 0.075` first and `p > 15/16` second, and piece C is the central routine's
  outermost piece, whose lower mirror is unreachable. Two typed thresholds, not a design.
- NORM.INV(p, m, s) = RN(m + RN(s * NORM.S.INV(p))) on 400/400: the inverse is stored before any consumer.
  No extended leak through NORM.INV.

**Doubt:** that the lattice means "a large cancelled intermediate in the numerator". A stored-double
*denominator* cannot give a 14.6-unit lattice (its relative granularity is at most 2 units), but a
denominator formed as `1 - stored(R(s))` or `C - stored(R(s))` can: the relative granularity is
`ulp(R)/D(s)` and grows smoothly as D shrinks. R43's lattice map says the spacing grows from 2.5 units at
s = 0.25 to 14.6 at s = 0.50 with no binade step. If D is linear near its root, `(s0 - 0.25)/(s0 - 0.5) =
5.84` gives **s0 = 0.552 (t = 0.743)**: a denominator with a simple zero just past the piece's end. That is a
classic low-degree rational `t * N(s) / (1 - s/s0)` or a `[m/1]`-like form, not a 7/7 table.

**Approaches.** (a) *Offline first:* fit `1/L(s)` on the 154 retained lattice windows; if it is linear in
s, read s0 to four digits and D's degree. (b) Then recover N(s) by exact preimage intersection with D
fixed (the R44 inverse-rounding tool, with the division inserted). (c) Fresh windows at s = 0.46, 0.48,
0.495 where the lattice must reach 20-30 units if s0 = 0.552; a flat lattice kills it in 768 questions.
(d) Piece B (1/sqrt2 < |t| <= 0.85) has sd 1.4-3 ulp: ask whether it is the same rational continued (its
lattice would shrink again past s0 as D grows negative) or a second piece; one 1,024-window at t = 0.75.
(e) Piece C is accurate and tiny (0.85, 0.875]: test the tail family in `z = sqrt(-ln(1 - t))` *and* in
`w = -ln(1 - t^2)` (Giles) by their exact preimages, since the range is so short that a two-term form may
already be exact.

## 4. YIELD start and stop, ODDFYIELD start

**Rounds assumed** Newton or secant from fixed or cash-dependent seeds, bracket programs, a scale-free
or absolute stop, and checked iterate membership to 16 ulp with analytic prices.

**New facts (e4, and the round-46 exponent map re-read).** For scale-free zero-coupon designs:
- Below 2^-100 the answers saturate; pr = 25 (root 1) returns **exactly 1.0**; pr = 80 returns the
  **negative root -2.118** (so the solver is an unbracketed Newton and the small-cash chain starts near
  y0 = 1 for N = 2: from 1, Newton on 100/(1+y)^2 lands on exactly 1.0 for pr = 25, on the negative root
  for pr = 80 and on the positive roots for 30, 37, 50, which is what Excel does).
- From e = 8 to 30 the deficit obeys `L + c/pr` (halves per exponent); from e = -12 to -2 it obeys
  `A + k*pr` (doubles per exponent); in between the iteration count changes. A seed that is a bounded
  smooth function of pr *with an absolute constant* (par 100) explains all three regimes: the two
  plateaus are the chains from y0(0) and y0(infinity), the linear laws are first-order sensitivity.
- Textbook seeds with analytic Newton do not reproduce the large-cash deficits (e4b harness; it prints
  candidate/Excel deficits per exponent). The regime shape is right, the seed or the slope is not.

**Approaches.** (a) *Backward Newton on the plateaus.* On the small-cash plateau the seed is y0(0), a
function of N and f only. The two-period price is quadratic, so each backward step is algebraic with at
most two preimages; k <= 8 gives <= 256 candidate seeds per design, and N = 2, 4, 8 designs intersect to
one formula (2/N?). Repeat on the large-cash plateau for y0(infinity). The crossover exponent then fixes
the constant (100 or 1). (b) *Slope.* With y0 known, y1 is predictable; an exact match says analytic,
an O(h) offset gives the finite-difference step. (c) *Stop constant.* The jitter-zone edges are where
the untaken step crosses 1e-10; two designs pin it. (d) *ODDFYIELD shares or not:* an odd first period
equal to a regular period must give ODDFYIELD = YIELD bit for bit if the solver is shared; run the same
exponent scan for it (cheap, 1,000 calls). (e) *#NUM! edge:* for pr below the root's basin, the price
where the first Newton step crosses 1 + y/f <= 0 is a sharp edge pinning (y0, P'(y0)) jointly; bisect it.

## 5. FORECAST.ETS search and seasonality

**New fact (e5).** Reported alpha in {0.002, 0.126, 0.25, 0.5, 0.75, 0.9, 0.998}, beta in {0.001, 0.125,
0.249, 0.499, 0.749}, gamma = 2^-52 exactly on all 28 readings: a coarse grid with 0.001/0.002 boundary
offsets and DBL_EPSILON as the lower clamp. This is the rounds' 25-pair table seen from the bits; the 1%
residual is near-tie behaviour of the score, so the score is the question, not the table.

**Doubts and probes.** (a) *The metric.* Build histories where SSE and sum of absolute errors (or MASE, or
AIC with the parameter count) rank two grid candidates oppositely; Excel's pick names the metric in 20
questions. (b) *Accumulation.* Histories whose two candidate scores differ by less than one ulp of the
sum, with the sign flipping between sequential binary64, pairwise, and PC64 accumulation; each history
is one bit. (c) *Tie order.* Exact ties by construction (histories whose residuals are identical under
two candidates, e.g. constant series with a step) separate `<` from `<=` and the enumeration order. (d)
*Detection is not ACF.* 81 ACF programs fail; Hyndman-style code picks the period by refitting. Probe:
histories where the ACF peak lag and the best-fit period disagree (two sinusoids with tuned amplitudes),
and a 1-ulp change of the *last* point only: an ACF decision barely moves, a fit-score decision can flip.
Also add a strong linear trend: raw-ACF and detrended-ACF detectors answer differently.

## 6. BINOM.DIST shared producers

Follows section 2 entirely. The one independent check worth 300 calls: BINOM.DIST(k, n, p, TRUE) against
BETA.DIST(1 - p, n - k, k + 1) **restricted to p whose 1 - p is exact** (p = m/2^j): if those agree bit for
bit everywhere, the "other route" is only the (x, y) interface and needs no separate search.

## 7. TRUNC, ROUNDDOWN, ROUNDUP: the fifteen-digit reading

**New facts (e7, e7b, 166 answers).**
- The half decision is **sharp at distance 0 to within one unit of 2^-64** of the scaled value in all
  four scale classes (one-chunk n = 34 and -6, two-chunk n = 314 and -286). In the n = 314 class a point
  0.97 units *below* the half went up (Excel's product is biased up by about one unit there); in the
  n = -286 class points 0.30-0.31 *above* went down and 0.32 split (biased down by about 0.3 units).
  Opposite signs for positive and negative n: two code paths (multiply vs divide), as S4 assumes.
- Without scaling (16 significant digits, 16th exactly 5): ROUNDDOWN, TRUNC and ROUNDUP all read the
  half **up regardless of parity** (…344.5 -> …345, …346.5 -> …347). ROUND(x, 0) on the same inputs
  returns …344: **ROUND reads the sixteenth digit differently** (drops it, or rounds the 16-digit value by
  another rule). Two readers exist.

**Approach.** Per n, bisect the decision threshold in units of 2^-64 (band 1e-5 in the design scan gives
0.05-unit resolution); the bias b(n) is a deterministic function of the 10^n construction and the
mantissa. Compare against candidate constructions (chunk order, table rounding, 64- vs 53-bit products)
on three n per chunk boundary; most die on the sign alone. Then re-score the 52 S4 losses.

## 8. CUMPRINC's second pv-dependent rounding

**New result (e6, offline on ld1's 5,000 consecutive pv).** Exact interval intersection over a hidden
constant A refutes six two-stage shapes `P = g(RN(pv*A))` with g in {identity, /1.75, *RN(1/1.75),
*RN64(1/1.75), minus RN(T*RN(r/(1+r))), minus RN(T*r)/(1+r)}: every intersection is empty. Since equal
public PMT gives different P (ld1), pv enters twice. **Next:** the difference shape `RN(RN(pv*A) -
RN(pv*B))` with A one of the five doubles around the exact PMT factor (so B is one-dimensional), then
`RN(RN(pv*A) * RN(pv*B)/pv)`-type forms only if that fails. The intersection tool is in e6 and takes
seconds per shape.

## Cross-cutting tools to build (small)

- A regime-replay harness per solver (e4b is the seed): candidates must reproduce root choice, plateau
  values and crossover exponents before any bit comparison.
- An exact-preimage intersection library over one or two hidden constants (e6), used for normalizers,
  seeds and denominators alike.
- A "relative oracle" checklist per function: which other Excel functions must equal it (mirrors,
  complements, pdf/cdf ratios, LN/EXP round trips), asked first, bit for bit.

## Three for the next Codex round

1. GAMMALN core signature (1a): binade jumps of the error sd across [0.7, 3) to read g, after landing
   the [3, 8) reduction and the 0.7 join (both proved here, zero fitted constants).
2. NORM.S.INV denominator (3a-3c): fit 1/L(s) offline, predict the lattice at s = 0.46..0.495, ask.
3. YIELD backward Newton on the two plateaus (4a) with the quadratic designs, then ODDFYIELD's scan (4d).
