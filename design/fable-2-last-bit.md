# Last-bit research, round fable-2

Claude Fable 5.1, Windows (GOVERT-OMEN), 2026-11-09. Implementer; sources: this folder's notes, reports 45-49,
ld-report-1, calc-func at DnaCalc origin/master (read only, never built), IEEE 754, the published methods already
cited there, and fresh answers from Excel 16.0.20527 (env f3269992527c) through the question service. Every
"proved" claim below predicted fresh Excel answers asked after the claim was fixed; emulations of calc-func use
Excel's own public LN, EXP, GAMMALN and POISSON.DIST as leaves, never a library. Numbers: `evidence/fable-2/`.

Order worked: POISSON.DIST, TRUNC family, CUMPRINC (as far as refutations), then out of time: F.DIST / GAMMA.DIST /
GAMMA.INV and DURATION / IRR / NPV / PRICE were not probed (section 5 says what transfers and what to ask first).

## 1. POISSON.DIST

### 1.1 Where the bank's 597 differences are (measured, p10/p11)

The 45,997 cache rows (Excel 20131) are 34,000 cumulative rows at k = 0 (= EXP(-m)) and 11,997 masses at
x in {1, 2, 3}, means 1e-6 to about 500. The current calc-func mass, emulated with Excel's leaves, reproduces
**11,996 of the 11,997 masses** (one miss: x = 1, m ~ 1e-6, 1 ulp). 312 re-asked cache rows equal today's bits
(no build drift). So the bank's failures are in the 3,584 W111 rows (Excel 20430): nearly all have means
700-760 (the normal/subnormal publication probes) and are cumulative rows of the incomplete gamma's asymptotic
branch (467 of 740 differ from the emulation, 309 of them within 3 ulp) plus mass rows whose published value is
subnormal (the emulation lacks calc-func's flush, so those are not real failures) and 488 Temme rows (a >= 20,
not emulated). **The mass is not the problem; the cumulative beyond m = a is.**

### 1.2 Mass (proved on fresh windows, nothing to land)

- x = 1 (first route): `exp_compensated(LN m, [m, GAMMALN(2) = 0])` with the register difference kept at 64 bits:
  **1,536/1,536** on twelve 128-windows from m = 0.3 to 300 (p9). An emulation that stores the difference to 53
  bits loses from m ~ 1 on (p3/p4): the lesson for anyone re-deriving the chain in Python. Excel's error against
  the exact m e^-m is the stored LN's (p5/p7: exact exponent after the public LN matches 85-98% alone).
- x = 2, 3, m < 1.1: the same chain with GAMMALN(3) = 3fe62e42fefa39f0, GAMMALN(4) = 3ffcab0bfa2a2002: 128, 121,
  128 / 128, 121, 106 of 128 at m ~ 0.3, 0.9, 1.05 (p8). The 1.05 window at x = 3 is the only soft spot.
- x = 2, 3, m >= 1.1 (Loader route, stored series deviance or stored `x LN(x/m) + (m - x)`): **1,023/1,024** in the
  series range (p9) and 768/768 outside (p8). Excel's mass there is inaccurate by tens of ulp (mean -25.7 ulp at
  x = 2, m = 40), and the chain reproduces that inaccuracy bit for bit.

### 1.3 Cumulative, 1 <= k <= 18 (a = k + 1 < 20)

Identities (p1, 462 rows): POISSON.DIST(k, m, TRUE) = CHIDIST(2m, 2k + 2) = CHISQ.DIST.RT(2m, 2k + 2) bit for bit in
every branch; it is RN(1 - GAMMA.DIST(m, k + 1, 1, TRUE)) only in part (89/108 Taylor, 7/98 between a and 31,
0/88 asymptotic, 112/168 Temme): P and Q are published separately, so GAMMA.DIST is a separate caller.

- **m < a (Taylor): proved.** P = st(M * ans), M the public mass at a (= POISSON.DIST(k + 1, m, FALSE)), ans the
  stored forward series (`c = st(c * st(m / rr)); ans = st(ans + c)` until `st(c / ans) <= 2^-53`), Q = 0.5 + (0.5 - P):
  **69/69 fresh** (k = 1..18, m from 0.3a to 0.95a, p16) and 16/16 W111. This is the current code.
- **a <= m < 31: open.** Excel's error against the exact Q: sd 2.9 ulp, max 17.5, uncorrelated with Q's size or k
  (p18). Refuted (p12, p16, p17, p18, p19; fresh rows): the binary64 finite sum e^-m sum m^j/j! forward (40/94),
  in PC64 (41/66), Horner (47/66); the mass at a times the finite backward series 1 + k/m + k(k-1)/m^2 + ... in
  thirteen store arrangements (<= 17/94); the mass at k times that series (0/252); the Cephes continued fraction
  with the mass (21/94); Morris's rcomp prefactor exp(st(a LN m) - m)/(a-1)! and five store variants (<= 10);
  the mass's first route (compensated or plain) as prefactor with asymptotic, fraction or /m tails (0); the
  Taylor series continued past a and complemented (<= 4/94). Note the W111 "closed-int" rows match the finite
  sum 5/8 and the old W109 note's 91% was on the 20131 cache rows (none in this region).
- **m >= 31: open.** Excel's error grows with m (sd 12 ulp below 100, 62 up to 450, 295 up to 700), the shape of
  a stored exponent of size m; yet every prefactor that carries such an error fails: the public mass at a with
  the NSWC asymptotic tail (23/110 fresh, 61/252 in p13), thirteen arrangements of mass x backward series
  (<= 29/110), the continued fraction (23/110), rcomp-like and first-route prefactors (<= 2), the finite sum
  (3/140), the mass at k (0/252). At means >= 700 (W111) the same chain differs by 1-3 ulp on 309 rows; a
  subnormal-grid store does not repair it (p13).

Killing probes for the next lane (a <= m): (a) hold m fixed and step k = 1..18: if Excel's Q(k) / Q(k-1) - 1
equals m^k/k! / (sum to k-1) to a rounding, the terms are summed forward; if Q(k) is one rounding of Q(k-1) + T_k
for a stored T_k, the sum is a stored running sum of stored terms (both are cheap to check with exact rationals);
(b) 256-windows at m = 12, 25, 40, 90 for k = 3: the sawtooth periods name the stored quantities (the method of
`x87_chains.rs`); (c) ask the same rows of GAMMA.DIST(m, a, 1, TRUE) and CHISQ.DIST(2m, 2a, TRUE): where P + Q = 1
exactly, P is RN(1 - Q_register) and the register Q is a leaf for the preimage tool.

### 1.4 Recipe (Codex lane)

Nothing to change for the mass. The cumulative's Taylor region is confirmed. For the bank, the honest expectation
is that the 597 POISSON.DIST differences are mostly the W111 asymptotic rows at means >= 700 and the Temme rows;
a lane should first run `cargo xtask` witness pins on `w111_poisson_rows.json` to confirm that census before
any kernel change, then take 1.3(a)-(c).

## 2. TRUNC, ROUNDDOWN, ROUNDUP: the fifteen-digit reading at |n| > 27 (proved)

Structure. For x with decimal exponent e10 and n = 14 - e10, |n| > 27, Excel's scaled value V is built from x's
exact binary value by the nibble chain, low chunk first: for each nonzero part c = (|n| >> 4i & 15) << 4i,
`V = RN64(V * K_c)` with `K_c = RN64(10^c)` when n > 0 and `K_c = RN64(10^-c)` when n < 0 (a multiplication, not a
division). The reading's 15-digit integer is floor(V) + [frac(V) >= 1/2] for the directed functions (ties up at
the 64-bit value). Everything else in `ExactDecimal::fifteen_x87` stays.

Evidence, all Excel 20527 through the service:
- fable-1's 160 half decisions (n = 34, 314, -286, -6): K2 160/160; the current one-constant K1 156/160 (t1).
- 81 designed x where the constructions disagree (n = 314, -286, 34, 154, -126; t2): **K2 81/81**; K1 54, K7
  (constants multiplied first) 59, high-chunks-first 55, divide-for-negative 64, truncating products 34.
- Tie classes (t4, t5): K2's V one 2^-64 unit below the half -> Excel down **105/105**; exactly the half -> up
  **108/108**; one unit above -> up **106/106** (n = 314, -286, 154, -126, 34).
- Fresh random x, 1e-300..1e-30 and 1e30..1e300, 240 per function (t6): the K2 reading equals Excel's for
  ROUNDDOWN, ROUNDUP and TRUNC **240/240** each (K1 too: random x never sits near a half).

The ld1 note's "S4 gains 396, loses 52 near exact halves" is not reproduced by any question here: K2 is S4 as the
research patch reads (low nibble first, RN64 table constants, RN64 products, 10^-c factors for n < 0), and Excel
goes up on every exact half of K2's V. If the lane sees losses on the bank, those rows contradict t5 and should
be sent back with their bits; a tie rule other than "half up at 64 bits" is excluded for ROUNDDOWN/TRUNC/ROUNDUP.
ROUND's reading (ties down in the current code) was not probed at these scales: keep its rule, switch only the
construction, and watch the 70 ROUND tie rows.

Recipe: in `fifteen_x87`, replace the `shift.abs() > 27` arm (one `round_to_64_bits` of the mantissa by the nearest
64-bit 10^|shift|) by the loop above over the nibbles of |shift| (the decomposition `scaled_pc64` already uses),
each factor `Self::unit(+-c).pc64_parts()`, each product `round_to_64_bits`; keep the existing half decision and
`ties_up`. Expected: the 200 larger differences per function (worst 86 ulp) and most of the near-half last-bit
rows; zero losses expected on the directed functions. Reference emulation: `t1_trunc_constructions.py::chain`.

## 3. CUMPRINC (and CUMIPMT through it): refutations and a fingerprint, nothing proved

Setting: the last period's principal CUMPRINC(r, n, pv, n, n, 0) with n = 1000, where expm1(-n L) is exactly -1,
so q = pv / em = -pv exactly and the exact answer is pv r / (1 + r). 300 consecutive pv at nine rates (c4), plus
ld1's 2,000-row scan at r = 0.75.

Measured fingerprint (c4): the effective constant -P/pv deviates from r/(1+r) by, in units of 2^-53 relative,
r = 0.75: -2.26, 0.5: -0.40, 0.25: +0.18, 0.1: -0.96, 0.05: +1.16, 0.3: +0.10, 1.5: +1.78, and exactly 0 at
r = 1 and 3 (1 + r a power of two). Excel's per-row error sd is 0.3-0.5 ulp. The deviation is neither the stored
LN(1 + r)'s rounding (c7: +0.24, -0.03, -0.08, -0.63, -0.39, -0.15, -0.37) nor EXP(-LN(1+r))'s (-0.5, -0.5, +0.5,
-0.26, -0.48, +0.42, -0.75), so a constant of the chain is off by one to two units more than any public leaf.

Refuted by exact preimage intersection over a free real constant (c5, c6; empty at every rate where 1 + r is not
a power of two, non-empty at r = 1 and 3): P = st(pv K); st(st(pv r) K); st(st(pv/(1+r)) K); st(st(pv r) - st(pv B));
st(st(pv r) - st(st(pv B1) r)); st(st(pv r) - st(st(pv r) B2)); st(T - st(T - st(pv B3))) with T = st(pv r). With fixed
public constants u = EXP(-LN(1+r)) and r (c3, 2,000 rows): nineteen two-stage products at best 163; the balance
difference st(st(q u_e) - st(q u_s)) 0 (c2). Together with ld1 (one rounding of pv x any real: 82%) and fable-1's
e6, every one- and two-rounding shape over pv with the loan's public leaves is dead; the second pv-dependent
rounding either involves a third constant or a constant two units from its public leaf.

Next probes (cheap, designed): (a) the same 300-pv scans at n = 2 and 3 (em far from -1): every intermediate is
then enumerable and the preimage tool (c6) can run three-stage shapes; (b) step the period k at fixed pv and
rate: if the deviation grows linearly with m_e = n - k + 1 the logarithm L is the off constant (then bisect it
with CUMIPMT, which shares the chain); (c) rates 2^-j: with 1 + r not a power of two but r exact, the deviation
table above gains points at which u, v = 1/(1+r) and r/(1+r) separate in their last bits.

## 4. What was not reached

F.DIST, GAMMA.DIST, GAMMA.INV; DURATION, MDURATION, IRR, NPV, PRICE: no question asked. Two things transfer:
GAMMA.DIST(x, a, 1, TRUE) with integer a shares the incomplete gamma of 1.3 (its Taylor region should be the
current code; a <= x is open and publishes P separately), and the mass's Loader route (1.2) is GAMMA.DIST's
density prefactor, so GAMMA.DIST's 171 last-bit rows are worth splitting by branch with `p16`'s method before
touching bratio for F.DIST. For the bond and rate solvers, start with the oracle relations (DURATION vs PRICE's
discount factors; IRR vs NPV at the returned rate) on 64-windows, as done here for POISSON vs CHIDIST.

## 5. For the Codex lane, in order

1. TRUNC/ROUNDDOWN/ROUNDUP nibble chain (section 2): proved, a recipe, a no-loss expectation.
2. POISSON.DIST census of the 597 (section 1.4), then the cumulative probes 1.3(a)-(c).
3. CUMPRINC probes 3(a)-(c) with the preimage tool from `c6`.
