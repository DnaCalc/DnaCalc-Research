fable-2: last-bit research (Claude Fable 5.1, Windows GOVERT-OMEN, 2026-11-09)

Excel 16.0.20527 (env f3269992527c) through the question service only, interactive
priority; the bank rows cited are the OxFunc cache (Excel 20131) and the W111 answer
files (Excel 20430) in CalcFunc-Evidence. calc-func was read at origin/master, never
built or run; every "current model" result below is a Python emulation of its chain with
Excel's own public LN/EXP/GAMMALN/POISSON.DIST leaves (common.py: st = RN64 then RN53).
Implementer role; no tier 3/4 source. Design note: design/fable-2-last-bit.md.

Scripts and what each showed (the .json beside a script holds the asked inputs, Excel's
bits and the build):
  p1_poisson_identities      POISSON.DIST(k,m,TRUE) == CHIDIST(2m,2k+2) == CHISQ.DIST.RT 462/462 in every
                             branch; == 1-GAMMA.DIST(m,k+1,1,TRUE) only 89/108, 7/98, 0/88, 112/168.
  p2_poisson_mass_oracles    the mass is not GAMMA.DIST's density nor a stored product of it.
  p3/p4/p5/p7                x = 1 windows: public-leaf programs, error map; the stored LN is the error source.
  p6_poisson_bank_drift      312 cache rows (20131) re-asked: 312 equal today's bits.
  p8_poisson_x23_windows     x = 2, 3: first route below 1.1, Loader route above (stored deviance); error maps.
  p9_poisson_series_and_fixed_d  x = 1 chain with the register difference at 64 bits: 1,536/1,536; series-range
                             Loader route x = 2, 3: 1,023/1,024.
  p10_poisson_bank_emulation current mass model over the 11,997 cache masses: 11,996 (one miss, x = 1, m ~ 1e-6).
  p11_poisson_w111_emulation the W111 rows: failures sit in the cumulative's asymptotic branch at means >= 700
                             and in subnormal flushes the emulation lacks.
  p12..p19                   the cumulative for 1 <= k <= 18: Taylor region proved (69/69 fresh, 16/16 W111);
                             a <= m regions open, 40+ refuted shapes (finite sums, backward series on the mass at a
                             or at k, continued fraction, rcomp-like and first-route prefactors, complements).
  t1_trunc_constructions     fable-1's 160 half decisions: K2 (nibble chain) 160/160; K1 (current) 156/160.
  t2_trunc_designed          81 points where constructions disagree: K2 81/81, K1 54, K7 59, K4 55, K3 64, K2z 34.
  t3 (no usable rows in the 20430 ROUNDDOWN cache); t4/t5 tie classes: K2 one unit below the half -> down
                             105/105, exactly the half -> up 108/108, one unit above -> up 106/106.
  t6_trunc_fresh_random      240 random x per function at |n| > 27: K2 reading 240/240 for ROUNDDOWN, ROUNDUP, TRUNC.
  c1 (stopped, too slow), c2, c3, c5, c6, c7  CUMPRINC singleton shapes by exact preimage intersection: all empty
                             or below 10% where 1 + r is not a power of two.
  c4_cumprinc_last_period_bias  9 rates x 300 pv: the effective constant's deviation per rate (fingerprint).
  w111_poisson_rows.json     the W111 POISSON.DIST rows (source, x, m, cumulative, Excel bits, id).
