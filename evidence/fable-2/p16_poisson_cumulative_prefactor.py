"""P16: POISSON.DIST(k, m, TRUE), 1 <= k <= 18, fresh rows in the three incomplete-gamma regions (m < a; a <= m < 31; m >= 31),
Excel's error against the exact Q, and prefactor hypotheses for the tail formulas:
 PF1 the public mass at a = k+1 (accurate);  PF2 Morris's rcomp: R = EXP(E) / (a-1)!, E = st(st(a LN m) - m), divided by a
 PF3 R with E = st(a LN m - m) formed in the registers;  PF4 R = st(EXP(E) * (1/(a-1)!)).
 Tails: Taylor (m <= max(a, ln 10)): P = st(PF ans), ans the stored forward series, Q = 0.5 + (0.5 - P);
        asymptotic (m >= 31): Q = st(st(PF st(a/m)) (1 + summ));  between: the closed-int finite sum from EXP(-m) (no prefactor)
        and, as alternatives everywhere, the continued fraction Q = st(st(PF a) CF)."""
from common import *
from decimal import Decimal, getcontext
import json, math, collections, statistics, random
from p15_poisson_cumulative_from_mass_k import exact_q, err_ulps, nswc_summ
from p14_poisson_cf import gamma_fraction
getcontext().prec = 80
random.seed(21)


def taylor_ans(a, x):
    rr = float(a); c = 1.0; ans = 1.0
    while True:
        rr += 1.0
        c = d(st(F(c) * st(F(x) / F(rr))))
        ans = d(st(F(ans) + F(c)))
        if d(st(F(c) / F(ans))) <= 1.1102230246251565e-16:
            return ans


rows = []
for k in [1, 2, 3, 5, 8, 12, 18]:
    a = k + 1
    for m in [0.3 * a, 0.7 * a, 0.95 * a, 1.2 * a, 2.0 * a, 12.0, 25.0, 30.5, 31.5, 45.0, 90.0, 200.0, 450.0]:
        for _ in range(3):
            mm = m * (1 + random.uniform(-0.04, 0.04))
            if mm > 0.05:
                rows.append((k, mm))
Q = vals(ask("POISSON.DIST", [[float(k), m, boolean(True)] for k, m in rows]))
MASS = vals(ask("POISSON.DIST", [[float(k + 1), m, boolean(False)] for k, m in rows]))
LN = vals(ask("LN", [[m] for k, m in rows]))
EXN = vals(ask("EXP", [[-m] for k, m in rows]))
E2 = [d(st(st(F(k + 1) * F(l)) - F(m))) for (k, m), l in zip(rows, LN)]
E3 = [d(st(rn(F(k + 1) * F(l), 64) - F(m))) for (k, m), l in zip(rows, LN)]
EX2 = vals(ask("EXP", [[e] for e in E2]))
EX3 = vals(ask("EXP", [[e] for e in E3]))
by = collections.defaultdict(collections.Counter); tot = collections.Counter(); errs = collections.defaultdict(list)
out = []
for (k, m), q, mass, l, en, e2, e3, x2, x3 in zip(rows, Q, MASS, LN, EXN, E2, E3, EX2, EX3):
    if not all(isinstance(v, float) for v in (q, mass, l, en, x2, x3)):
        continue
    a = k + 1
    reg = "m<a" if m < a else "a<=m<31" if m < 31 else "m>=31"
    tot[reg] += 1
    errs[reg].append(err_ulps(q, exact_q(a, m)))
    g = float(math.factorial(a - 1))
    pf = {"PF1": F(mass), "PF2": st(st(F(x2) / F(g)) / a), "PF3": st(st(F(x3) / F(g)) / a), "PF4": st(st(F(x2) * F(1.0 / g)) / a),
          "PF5": st(F(x2) / F(g * a)), "PF6": st(F(x2) / st(F(g) * a))}
    res = {}
    if m <= max(a, 2.302585092994046):
        ans = taylor_ans(a, m)
        for n_, p in pf.items():
            P = st(p * F(ans)); res["tay " + n_] = 0.5 + (0.5 - d(P))
    if a <= m < 31:
        s = en; t = en; c = 0.0
        for n_ in range(1, a):
            c += 1.0; t = (m * t) / c; s += t
        res["sum f64"] = s
        hh = 1.0
        for j in range(k, 0, -1):
            hh = 1.0 + hh * j / m
        res["sum horner"] = en * hh
    if m >= 31:
        summ = nswc_summ(a, m)
        for n_, p in pf.items():
            res["asy " + n_] = d(st(st(p * st(F(a) / F(m))) * F(1.0 + summ)))
    if m > a:
        cf = gamma_fraction(a, m)
        for n_, p in pf.items():
            res["cf " + n_] = d(st(st(p * a) * cf))
    for n_, v in res.items():
        by[reg][n_] += v == q
    out.append({"k": k, "m": h(m), "cum": h(q), "mass_a": h(mass), "ln": h(l), "exp_neg_m": h(en), "E2": h(e2), "exp_E2": h(x2), "E3": h(e3), "exp_E3": h(x3)})
for reg in sorted(tot):
    e = errs[reg]
    print(f"{reg:8} n={tot[reg]:3} excel err mean {statistics.mean(e):+7.2f} sd {statistics.pstdev(e):7.2f} max {max(abs(v) for v in e):7.2f}")
    for n_, v in sorted(by[reg].items(), key=lambda t: -t[1]):
        print(f"    {n_:12} {v:4}")
save("p16_poisson_cumulative_prefactor.json", {"rows": out, "build": meta(ask("LN", [[2.0]]))})
