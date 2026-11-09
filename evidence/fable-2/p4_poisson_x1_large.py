"""P4: POISSON.DIST(1, m, FALSE) where the compensated exp(ln m - m) fails (m above about 1). Candidates:
 H1  st(m * v64), v64 = RN64(exp(-m)) exact (an extended exp kept in the registers, one store)
 H1b st(m * EXP(-m)) (public, stored)   [= A of p3]
 H2  Loader: e = st(-s1 - dev), dev = st(st(LN(st(1/m))) + st(m - 1)) outside the series range; mass = st(EXP(e) / sqrt(2 pi))
 H3  Loader with EXP(e) kept extended: st(RN64(exp(e)) / sqrt(2 pi))
 H4  st(RN64(exp(st(LN m - m))) * st(1 + dropped))  (compensated D with an extended exp)
 H5  RN64 exp of the exact ln m - m (no stored argument): st(m * exp(-m)) computed as RN53(RN64(exact m e^-m))
"""
from common import *
from decimal import Decimal, getcontext
getcontext().prec = 60
S1 = 0.0810614667953272
SQ2PI = float(F(2) * F(3.141592653589793)) ** 0.5   # Rust: (2.0 * PI * x).sqrt() in binary64


def dexp(q):
    """exp of a Fraction as an exact Fraction to 60 digits (then rounded by the caller)."""
    v = Decimal(q.numerator) / Decimal(q.denominator)
    return F(v.exp())


def rn64(q):
    return rn(q, 64)


centres = [1.05, 1.3, 2.0, 3.3, 8.0, 40.0, 300.0]
N = 128
ms = [nxt(c, i) for c in centres for i in range(N)]
mass = vals(ask("POISSON.DIST", [[1.0, m, boolean(False)] for m in ms]))
ln = vals(ask("LN", [[m] for m in ms]))
exn = vals(ask("EXP", [[-m] for m in ms]))
inv = [d(st(F(1) / F(m))) for m in ms]
lninv = vals(ask("LN", [[v] for v in inv]))
es = []
for m, li in zip(ms, lninv):
    dev = d(st(F(li) + st(F(m) - 1)))
    es.append(d(st(-F(S1) - F(dev))))
exe = vals(ask("EXP", [[e] for e in es]))
res = {}
rows = []
for m, p, l, e, li, ee, exee in zip(ms, mass, ln, exn, lninv, es, exe):
    c = [c for c in centres if c <= m < nxt(c, N)][0]
    v64 = rn64(dexp(-F(m)))
    H1 = d(st(F(m) * v64))
    H1b = d(st(F(m) * F(e)))
    H2 = d(st(F(exee) / F(SQ2PI)))
    H3 = d(st(rn64(dexp(F(ee))) / F(SQ2PI)))
    t = d(st(F(l) - F(m)))
    dropped = d(st(st(F(l) - F(t)) - F(m)))
    H4 = d(st(rn64(dexp(F(t))) * st(1 + F(dropped))))
    H5 = d(st(F(m) * dexp(-F(m))))
    r = res.setdefault(c, [0] * 7)
    r[0] += 1
    for i, v in enumerate((H1, H1b, H2, H3, H4, H5)):
        r[i + 1] += (v == p)
    rows.append({"m": h(m), "mass": h(p), "ln": h(l), "exp_neg_m": h(e), "ln_inv": h(li), "loader_e": h(ee), "exp_loader_e": h(exee),
                 "H1": H1 == p, "H1b": H1b == p, "H2": H2 == p, "H3": H3 == p, "H4": H4 == p, "H5": H5 == p})
print("centre   n   H1  H1b   H2   H3   H4   H5")
for c, r in res.items():
    print(f"{c:7} {r[0]:3} " + " ".join(f"{v:4}" for v in r[1:]))
save("p4_poisson_x1_large.json", {"rows": rows, "build": meta(ask("POISSON.DIST", [[1.0, 1.0, boolean(False)]]))})
