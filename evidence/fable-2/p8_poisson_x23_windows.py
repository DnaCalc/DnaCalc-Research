"""P8: POISSON.DIST(x, m, FALSE) at x = 2, 3 on windows of 128 adjacent m. Leaves asked of Excel: LN(m), LN(st(x/m)), EXP of
stored exponents. Candidates (G = Excel's GAMMALN(x+1); s_x the Stirling-error table; exact exp/ln to 60 digits):
 R1 st(exp(st(x LN m) - m - G))        R2 st(exp(x LN m - m - G))   R3 calc-func compensated public-EXP chain
 R4 st(exp(st(st(x LN m) - m) - G))    (binary64 stores after each subtraction, accurate exp)
 L1 calc-func Loader: e = st(-s_x - st(st(x LN(st(x/m))) + st(m - x))), st(EXP(e)/sqrt(2 pi x))
 L2 Loader, exponent exact after the public LN(st(x/m)): st(exp(-s_x - (x LN(x/m) + (m - x)))/sqrt(2 pi x))
 L3 L2 with the exact ln(x/m)."""
from common import *
from decimal import Decimal, getcontext
import math, statistics
getcontext().prec = 60
S = {2: 0.0413406959554092, 3: 0.0276779256849983}
G = {2: fb("3fe62e42fefa39f0"), 3: fb("3ffcab0bfa2a2002")}


def dexp(q):
    return F((Decimal(q.numerator) / Decimal(q.denominator)).exp())


def dln(q):
    return F((Decimal(q.numerator) / Decimal(q.denominator)).ln())


def exact_mass(x, m):
    v = Decimal(m)
    return F(v ** x * (-v).exp() / math.factorial(x))


def err_ulps(answer, truth):
    e = math.frexp(answer)[1]
    return float((F(answer) - truth) / F(2) ** (e - 53))


def sq2pix(x):
    return float(F(2) * F(3.141592653589793) * x) ** 0.5


centres = [0.3, 0.9, 1.05, 1.3, 2.6, 5.0, 40.0]
N = 128
rows = [(x, nxt(c, i)) for x in (2, 3) for c in centres for i in range(N)]
mass = vals(ask("POISSON.DIST", [[float(x), m, boolean(False)] for x, m in rows]))
ln = vals(ask("LN", [[m] for x, m in rows]))
ratio = [d(st(F(x) / F(m))) for x, m in rows]
lnr = vals(ask("LN", [[q] for q in ratio]))
# stored exponents for the public-EXP programs
e_first = []; e_loader = []
for (x, m), l in zip(rows, ln):
    t = d(st(F(x) * F(l)))
    e_first.append((t, d(st(F(t) - F(m)))))
for (x, m), lr in zip(rows, lnr):
    dev = d(st(st(F(x) * F(lr)) + st(F(m) - x)))
    e_loader.append(d(st(-F(S[x]) - F(dev))))
ex_first = vals(ask("EXP", [[d(st(F(n1) - F(G[x])))] for (x, m), (t, n1) in zip(rows, e_first)]))
ex_loader = vals(ask("EXP", [[e] for e in e_loader]))
res = {}; out = []
for (x, m), p, l, lr, (t, n1), ef, el, exl in zip(rows, mass, ln, lnr, e_first, ex_first, e_loader, ex_loader):
    key = (x, [c for c in centres if c <= m < nxt(c, N)][0])
    g = F(G[x])
    R1 = d(st(dexp(F(t) - F(m) - g)))
    R2 = d(st(dexp(F(x) * F(l) - F(m) - g)))
    # R3: compensated chain as calc-func: terms [m, G]
    stored = F(t); factors = []
    for term in (F(m), g):
        nxt_ = st(stored - term)
        dropped = st(st(stored - nxt_) - term)
        factors.append(st(1 + dropped)); stored = nxt_
    v = F(ef)
    for f_ in factors:
        v = st(v * f_)
    R3 = d(v)
    n2 = st(F(n1) - g)
    R4 = d(st(dexp(n2)))
    L1 = d(st(F(exl) / F(sq2pix(x))))
    dev_ext = F(x) * F(lr) + (F(m) - x)
    L2 = d(st(dexp(-F(S[x]) - dev_ext) / F(sq2pix(x))))
    L3 = d(st(dexp(-F(S[x]) - (F(x) * dln(F(x) / F(m)) + (F(m) - x))) / F(sq2pix(x))))
    r = res.setdefault(key, {"n": 0, "R1": 0, "R2": 0, "R3": 0, "R4": 0, "L1": 0, "L2": 0, "L3": 0, "err": []})
    r["n"] += 1
    for name, val in (("R1", R1), ("R2", R2), ("R3", R3), ("R4", R4), ("L1", L1), ("L2", L2), ("L3", L3)):
        r[name] += val == p
    r["err"].append(err_ulps(p, exact_mass(x, m)))
    out.append({"x": x, "m": h(m), "mass": h(p), "ln": h(l), "ln_ratio": h(lr), "exp_first": h(ef), "exp_loader": h(exl),
                "R1": R1 == p, "R2": R2 == p, "R3": R3 == p, "R4": R4 == p, "L1": L1 == p, "L2": L2 == p, "L3": L3 == p})
print("x centre   n   R1   R2   R3   R4   L1   L2   L3 | excel err mean sd min max")
for (x, c), r in res.items():
    es = r["err"]
    print(f"{x} {c:6} {r['n']:3} " + " ".join(f"{r[k]:4}" for k in ("R1", "R2", "R3", "R4", "L1", "L2", "L3")) +
          f" | {statistics.mean(es):+5.2f} {statistics.pstdev(es):4.2f} {min(es):+5.2f} {max(es):+5.2f}")
save("p8_poisson_x23_windows.json", {"rows": out, "build": meta(ask("POISSON.DIST", [[1.0, 1.0, boolean(False)]]))})
