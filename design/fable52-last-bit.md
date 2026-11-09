# Last-bit functions: structure, proofs and recipes (fable-52)

Claude Fable 5.2 (claude-fable-5-2), Windows (GOVERT-OMEN), 2026-11-09. Implementer role. Sources: the design
notes in this folder, `fresh-ideas-fable-1.md`, calc-func at DnaCalc `origin/master` 384a90f1 (read only), IEEE 754,
and fresh answers from Excel 16.0.20527 through the question service only (env f3269992527c and 29bedb2b6722).
No tier 3/4 source. Scripts, frozen predictions and replies are in `evidence/fable-52/`.

The fable-2 lane was working POISSON.DIST, CUMPRINC and TRUNC/ROUNDDOWN/ROUNDUP when this was written (its
`evidence/fable-2/` had probes p1-p9, c1, t1-t4; its design note did not exist yet). I did not repeat those three.
This note covers NPV, IRR, MDURATION, BETADIST, F.DIST, GAMMA.DIST, GAMMA.INV and HYPGEOM.DIST.

Words: **proved** = a claim written down first (predictions saved to a `*_frozen.json` file before the question was
sent) that then predicted fresh Excel answers. **Measured** = an identity between Excel's own functions, fixed in the
script before asking. **Hypothesis** = everything else, with its killing probe. `st(x)` = the x87 store: round to a
64-bit significand, then to binary64. "Stored" means every operation's result goes through `st`.

## Summary

| Function (bank last-bit count) | Result | Status |
|---|---|---|
| NPV (94) | Different loop from calc-func's; every operation stored | proved, recipe below |
| MDURATION (147) | `DURATION / (1 + yld/freq)`, nothing of its own | measured 300/300; follows DURATION |
| BETADIST (149) | Bit-identical to BETA.DIST, with and without bounds | measured 300/300 each; follows the beta |
| F.DIST (189) | Same body as BETA.DIST; differs only in how the argument pair is formed | measured; interface is a hypothesis |
| GAMMA.DIST (172) | cdf scale is `st(x/beta)` then the unit-scale body; = CHISQ.DIST; pdf scale differs | measured |
| GAMMA.INV (139) | `st(beta * unit-scale inverse)`; CHISQ.INV(p, 2a) = GAMMA.INV(p, a, 2) | measured 300/300 |
| HYPGEOM.DIST (89) | Not symmetric in (sample, successes); single draw is not a quotient; an exp of logs | measured; body is a hypothesis |
| IRR (174) | Secant in v confirmed against rate; conversion confirmed; objective arithmetic not found | about 110 programs refuted |

One general fact from NPV that the other financial kernels should be re-read against: this kernel stores **every**
operation (multiply, divide and add each showed the two-step rounding on designed inputs, 200/200 against 0/200 for
one-step binary64). That is the pattern of a helper call per operation, not of expressions kept in registers, and it
comes with per-operation error checks (an infinite power gives `#DIV/0!` before it is used, section 1).

## 1. NPV: proved

**Shape.** calc-func today (`financial_time_value_family.rs:490`) runs `discount *= base; total += flow / discount`
in plain binary64. Excel does this:

```
b = st(1 + rate)
p = 1;  t = 0
for each value v, in order:
    t = st(t + st(v / p))
    p = st(p * b)
NPV = st(t / b)
```

The flows are discounted as if the first were at time 0, and the sum is divided by `1 + rate` once at the end.

**How it was found** (shape first, no ranking by match count).
- N1: 56 loop programs (divide by a growing power, multiply by a shrinking one, backward Horner; every variable
  stored or in the register) all contradicted on 300 random calls.
- N2/N3: by cashflow count, one flow is `RN(v / b)`; then the flow was isolated: `NPV(r, 0, ..., 0, v)` with `v` in
  position k = 2..8, 40 calls each, asked beside Excel's POWER.
- N4/N5: position 2 is exactly `RN(RN(v/b)/b)` (40/40) but position 3 is not the three-fold chain (20/40). The only
  form that holds for every k is `(v / b^(k-1)) / b` with the power by sequential stored products: 280/280.
  Binary powering fails from k = 5 (38, 28, 22, 25 of 40).

**Proving tests** (predictions frozen before asking).
- N6, seed 5206: 400 random calls plus 400 edge calls (rates from 1e-9 to 1e12, negative bases, 1 to 30 flows, zeros
  among the flows): 799/800 exact; the earlier 300 discovery calls replay 300/300. The one miss is an error rule
  (below), not arithmetic.
- N9: 80 designed calls where `b*b` rounds differently in one step and in two: stored 80/80, plain binary64 0/80.
- N10: the same design for the division (60 calls, `NPV(b-1, v)`) and the addition (60 calls, `NPV(0, v1, v2)`):
  stored 60/60 and 60/60, plain 0/60 and 0/60.

**Error rules (measured on 15 designed calls, N7 and one inline batch; not a full map).**
- `b = 2^100`, flows of 1 or 2^900: 11 flows give a number (the power overflows only after its last use); 12 or 13
  give `#DIV/0!` although `v / inf` would be 0. So an infinite `p` at the moment it is used is `#DIV/0!`.
- `b = 2^-30`: 35 flows of 1 give `#NUM!` (the final `t / b` overflows); 36 or more give `#DIV/0!` (the term
  `v / p` overflows at the subnormal power 2^-1050); 37 flows of 0 give `#DIV/0!` (`0 / 0` once `p` is 0).

**Recipe for a Codex lane.**
1. Replace the loop in `npv` by the program above using `excel_x87_add`, `excel_x87_div`, `excel_x87_mul`.
2. Inside the loop, before the term: `p` infinite or zero gives `#DIV/0!`; a non-finite term gives `#DIV/0!`. After
   the loop a non-finite result gives `#NUM!` (as today). `b == 0` stays `#DIV/0!`.
3. Tests: the 1,100 calls of `n1_npv_discover_5201.json` and `n6_npv_proof_5206.json` (inputs and Excel bits as hex),
   the 200 designed calls of N9/N10, the edge calls of N7. Then the usual no-loss corpus and bank check; NPV's bank
   rows that match today do so because typical rates and few flows rarely separate the two loops, so a loss there
   would be a finding worth reporting, not a reason to keep the old loop.
4. Open, each one question: a subnormal `p` with a finite term (`b = 2^-30`, 36 flows of 2^-100; the value rule says
   a number, a "power must be normal" rule says `#DIV/0!`); the first operation `1 + rate` on a rate where one-step
   and two-step rounding differ (assumed stored like the others, not separately probed).

## 2. IRR: what holds, what is refuted

calc-func's `irr_discount_secant` (`cashflow_rate_family.rs:460`): secant in `v = 1/(1 + rate)` from `1/(1 + guess)`,
second point 0.001 away on the side of the NPV's sign, stop when the step and the NPV are both under 1e-7, answer
`1/v - 1`, objective `sum flow / POWER(1/v, k)`.

Fresh evidence (I1: 160 random problems of 3 to 9 flows with guesses; I7: 200 two-flow problems; all frozen first):
- **Solver variable is v, proved against the alternative**: the same solver on the rate matches 1 of 153 numeric
  answers and is off by up to 1e-1 relative; every v variant matches 82 to 98. All seven error cells agree.
- **Final conversion** `st(st(1/v) - 1)`: the forms `(1 - v)/v` (30/153) and "quotient kept in the register"
  (30 to 38/153) are refuted relative to it. Not proved exact, because no objective is.
- **Objective arithmetic: not found.** Refuted on the 153 answers, each with at least 55 contradictions: NPV's loop
  with base `1/v` (the shared-helper hypothesis from section 1), binary powering both ways, `t += flow*q; q *= v`,
  Horner in v, `p /= v`, and all 16 stored/register masks of each (I5: best 98/153). I4 shows why match counts
  cannot steer here: every candidate lands on Excel's iterate to within a few units of the rate; the last secant step
  is driven by the rounding residue of the objective, so any wrong operation costs the same third of the cases.
- **Two flows already fail** (I7-I9): `P(v) = f0 + f1 v` has no power at all, yet `f0 + f1*v` gives 130/200,
  `f0 + f1/(1/v)` 149/200, and 63 programs over three objectives, seven update expressions and three stop rules
  top out at 152/200. In 32 of 200 cases no iterate of the best program equals Excel's answer at all.

So the unidentified piece is in the two-flow program: the objective's linear term, the update expression, or the
second point. Hypotheses and killing probes, in the order I would run them:
- H1 *the second point is not `v0 +- 0.001` in v* (for example `1/(1 + guess +- 0.001)`): two-flow problems are
  linear, so the start barely moves the answer, but it moves the residue. Probe: two-flow problems with f0 = -1 and
  f1 = 2^j (root v = 2^-j exact, every product exact): the objective is then exact and the only rounding left is in
  the update, so the answers enumerate the update expression directly. Then repeat with non-dyadic f1.
- H2 *the objective is evaluated on the rate*, `b = st(1 + st(st(1/v) - 1))`: identical to `st(1/v)` only while the
  subtraction is exact; choose guesses and flows with `1/v` above 2 and below 0.5, where it is not.
- H3 *the stop compares other quantities* (step on the rate, or relative): I4's ranges for "Excel one iterate later"
  overlap on every quantity I logged, so no single threshold exists on those; design problems whose step crosses
  1e-7 between two iterates with the NPV already below 1e-12, 20 problems.
The harness for all three is `i9_irr_two_programs.py` (offline replay of saved answers, any number of programs).

## 3. MDURATION, BETADIST: nothing of their own

- **MDURATION** (R1, 300 random bonds, formula fixed before asking): `MDURATION = st(DURATION / st(1 + st(yld/freq)))`
  with Excel's own DURATION as the leaf, 300/300. The three ways to write the divisor agree for every valid
  frequency (1, 2, 4: scaling by a power of two is exact; R2 found no separating input in 400 bonds), so there is
  nothing to choose. MDURATION's bank failures are DURATION's. Recipe: `mduration_kernel` must take the published
  DURATION double (not an internal extended value) and do one stored divide; if it already does, close it under
  DURATION.
- **BETADIST** (R1): `BETADIST(x,a,b) = BETA.DIST(x,a,b,TRUE)` 300/300; with bounds 300/300; and
  `BETA.DIST(x,a,b,TRUE,A,B) = BETA.DIST(st(st(x-A)/st(B-A)),a,b,TRUE)` 300/300. Recipe: one body; the bounded form
  is the stored quotient in front of it. BETADIST's 149 are the beta's.

## 4. F.DIST: same body as the beta, different argument pair

Measured (R1, R2, R3):
- `FDIST = F.DIST.RT` 300/300.
- With `z = d1 x/(d1 x + d2)` dyadic so that `z` and `1 - z` are both exact (120 designed calls):
  `F.DIST(x,d1,d2,TRUE) = BETA.DIST(z, d1/2, d2/2, TRUE)` 119/120 and
  `F.DIST.RT = BETA.DIST(1 - z, d2/2, d1/2, TRUE)` 119/120.
- On general inputs the same identities hold on 204/300 and 202/300 only (up to 32 and 67 units), and
  `F.DIST = 1 - F.DIST.RT` fails on a third of them with differences of up to 5e14 units: the two tails are computed
  separately, not as complements.
- R3 (400 calls): four worksheet compositions per tail (`B(z)`, `1 - B(w)`, `B(1 - w)`, `1 - B(1 - z)`, with
  `w = d2/(d2 + d1 x)`) leave 12 cumulative and 23 right-tail answers matching none. `w = st(1 - z)` on 174/400.

Hypothesis F1: F.DIST calls the beta body with the pair `(x, y) = (z, w)`, both formed by their own quotient, where
BETA.DIST can only pass `(x, st(1 - x))`. It explains exact agreement when the complement is exact, partial
agreement otherwise, and the unmatched cases (neither worksheet pair is `(z, w)`). It cannot be tested through the
worksheet; the killing test is in calc-func once the beta body is exact: run bratio with `(z, w)` on R3's 400
inputs. A cheaper partial test today: the 174 cases of R3 with `w = st(1 - z)` must all equal `BETA.DIST(z, ...)`
if F1 is right (the script's table has the counts by side; not yet cut by that condition).
Consequence for ranking: F.DIST (189), FDIST and BETADIST (149) wait on the beta body; none needs its own search.
The one miss in each 119/120 is retained in the replies and should be looked at first when the body is exact.

## 5. GAMMA.DIST and GAMMA.INV: scale is outside the body

Measured (R1, 300 points each, shapes 0.5 to 37):
- `GAMMADIST = GAMMA.DIST`; `CHISQ.DIST(2x, 2a, TRUE) = GAMMA.DIST(x, a, 1, TRUE)` 300/300;
  `GAMMA.DIST(y, a, 3, TRUE) = GAMMA.DIST(st(y/3), a, 1, TRUE)` 300/300. The cdf divides by the scale once, stored,
  and runs a unit-scale body shared with the chi-square.
- The **pdf does not factor that way**: `GAMMA.DIST(y, a, 3, FALSE) = st(GAMMA.DIST(st(y/3), a, 1, FALSE) / 3)` on
  177/300 only, within 3 units. The scale enters the pdf's expression (inside the exponent or as a different
  operation order). Killing probes for the three natural forms: scale 2^j (all forms agree, control), then scales
  3, 10, 0.1 at shapes 1 (pdf = exp(-y/b)/b, no power) and 2 (one factor y/b): shape 1 alone separates
  `exp(-(y/b))/b`, `exp(-y/b - ln b)` and `(1/b) exp(-(y/b))`.
- `GAMMA.INV(p, a, 3) = st(3 * GAMMA.INV(p, a, 1))` 300/300; `CHISQ.INV(p, 2a) = GAMMA.INV(p, a, 2)` 300/300.
  The inverse is a unit-scale search followed by one stored multiply.
- `GAMMA.DIST(GAMMA.INV(p, a, 1), a, 1, TRUE) = p` on 132/300 (within 7 units): the search does not finish on the
  public cdf's exact preimage; its stop rule is open.

Recipes: (i) make CHISQ.DIST, GAMMADIST and GAMMA.DIST's cdf one body behind `st(x/beta)`, and CHISQ.INV, GAMMAINV,
GAMMA.INV one search behind `st(beta * q)`; any last-bit count that differs between those names in the bank is then
a routing defect, visible without research. (ii) Research then needs only the unit-scale cdf and the unit-scale
inverse. For the inverse, ask first whether the answer is a fixed point of a Newton or bisection step on the
public cdf: `GAMMA.DIST` at the answer's two neighbours tells which neighbour a search stopping on "cdf brackets p"
would return.

## 6. HYPGEOM.DIST: an asymmetric exp-of-logs

Measured (R1, R2):
- `HYPGEOMDIST = HYPGEOM.DIST(.., FALSE)` 300/300.
- Swapping sample size and successes (mathematically equal) agrees on **7/300**, up to 533 units. Swapping
  successes and failures agrees on 195/300, always within 1 unit. So the body treats the sample size differently
  from the success count, and nearly symmetrically in successes against failures.
- `cdf(k) = st(cdf(k-1) + pmf(k))` on 116/300, up to 268 units: the cdf is not a forward running sum of the pmf.
- Single draw, `HYPGEOM.DIST(1, 1, K, N, FALSE) = K/N` mathematically (200 calls, N up to 5,000): equal to `RN(K/N)`
  on 37, errors from -9 to +10 units; equal to Excel's `EXP(LN(K) - LN(N))` on 33; equal to an exp of
  log-factorial sums built from Excel's GAMMALN on 0 (two groupings).

Reading: an error of about ten units on a value below 1 is what an exponential of a stored double of magnitude 4 to
16 gives (half a unit of 8 is 9e-16). Sums of GAMMALN of the population (magnitude 4e4) would be thousands of units
off, and are. So the body is `exp` of a short sum of *log-binomial-sized* terms, from a helper that is not the public
LN of the count.
Hypothesis G1: `pmf = exp(lc(K, k) + lc(N-K, n-k) - lc(N, n))` with a log-choose helper `lc` evaluated by a product
or series over `min(k, n-k)`-like counts. The 533-unit asymmetry then comes from `lc(N, n)` against `lc(N, K)`
taking different numbers of steps. Killing probes:
- error scale against N at fixed K/N (N = 10, 100, 1e3, 1e5, 1e7, 50 calls each): flat for a log helper of
  magnitude ln N, growing like N ln N for log-factorials;
- separability: the implied `ln(pmf)` for single draws must be `a(K) - a(N)` with one `a` per integer, within the
  exp's rounding interval; intersect over a 40 x 40 grid of (K, N). An empty intersection kills every "difference of
  two per-integer logs" form at once; a non-empty one gives `a(K)` to compare with LN, GAMMALN differences and
  LN of COMBIN;
- asymmetry against steps: `pmf(k; n, K, N)` minus its swap for n = 1..50 at fixed K: linear growth in n says a
  product of n factors, a jump at a threshold says two routes.
The earlier rounds' "single-draw quotient fails 59 of 86" and "thirty quotient programs" agree with this: no
quotient program can match, because the value is not a quotient.

## 7. Not worked here

- POISSON.DIST, CUMPRINC, TRUNC/ROUNDDOWN/ROUNDUP: fable-2's lane (see its note when it lands). One cross-check for
  it from section 5: `CHISQ.DIST = GAMMA.DIST` bit for bit, and round 46 found every Poisson window equal to the
  public chi-square tail, so the Poisson cdf, the chi-square and the gamma cdf are one unit-scale body with three
  callers; a recipe that fixes one must be run against the other two banks.
- DURATION: only the relation in section 3. The design notes' open item (off-coupon discount arithmetic) stands.
  NPV's result suggests re-reading it with "every operation stored, powers by sequential stored products" before
  any POWER-based form: DURATION's coupon sum is the same kind of loop.

## 8. Method notes for the next lane

- Isolate one term by zeroing the others (`NPV(r, 0, ..., 0, v)`): the sum disappears and the power's construction
  is read directly. The same works for any cashflow or coupon loop.
- Separate one-step from two-step rounding by construction, not by sampling: 839 random NPV calls contained no
  separating case; 200 designed ones settled it. Search for operands with `RN53(x op y) != RN53(RN64(x op y))`,
  then build the call around them.
- Ask the relatives first: seven of the twelve functions on the list reduce to another function plus one or two
  stored operations. Those reductions were each one batch of 300 calls.
- Freeze predictions in a file before the question goes out; every `*_frozen.json` here was written that way.
