"""P9: (a) x = 1: the compensated chain with the register difference kept at 64 bits (as calc-func does; p3/p4's D stored it to
53 bits, which loses the log's low bits once m exceeds a few units); (b) x = 2, 3 inside the deviance's series range
(|x - m| < 0.2 (x + m)): the Loader route with the stored atanh series (ported from the design note's description), EXP public.
Fresh windows at x = 2: m = 1.5, 2.2, 2.9; x = 3: m = 2.1, 3.6, 4.4; and p8's m = 2.6 windows."""
from common import *
import json, math, statistics
S = {2: 0.0413406959554092, 3: 0.0276779256849983}


def r64(q):
    return rn(q, 64)


def comp_chain(t, terms, exp_of):
    """calc-func exp_compensated: next = st(stored - term); dropped = st(rn64(rn64(stored - next) - term)); factor = st(1 + dropped)."""
    stored = F(t); factors = []
    for term in terms:
        nx = st(stored - F(term))
        dropped = st(r64(r64(stored - nx) - F(term)))
        factors.append(st(1 + dropped)); stored = nx
    v = F(exp_of(d(stored)))
    for f_ in factors:
        v = st(v * f_)
    return d(v)


def atanh_tail(v):
    stored = st(v)
    w = st(stored * stored)
    tail = F(0); power = stored; j = 1
    while True:
        power = st(power * w)
        nx = st(tail + r64(power / (2 * j + 1)))
        if nx == tail:
            return tail
        tail = nx; j += 1


def deviance_series(x, m):
    fx, fm = F(x), F(m)
    dd = st(fx - fm)
    v = r64(dd / st(fx + fm))
    tail = atanh_tail(v)
    quotient = st(v)
    s = st(quotient + tail)
    twice_x = F(2 * x)
    gap = st(fm - fx)
    near = r64(gap + st(twice_x * s))
    near = r64(near - st(twice_x * quotient))
    near = r64(near + dd)
    product = st(dd * quotient)
    return st(near + product)


def in_series(x, m):
    dd = abs(r64(F(x) - F(m)))
    return r64(F(0.2) * r64(F(x) + F(m))) - dd > 0


def sq2pix(x):
    return float(F(2) * F(3.141592653589793) * x) ** 0.5


# (a) x = 1 on the retained p3/p4 answers: need EXP(next) which p3/p4 recorded as exp_next (p3) / implicit (p4: ask again)
rows1 = json.load(open("p3_poisson_x1_windows.json"))["rows"] + json.load(open("p4_poisson_x1_large.json"))["rows"]
nexts = [d(st(F(fb(r["ln"])) - F(fb(r["m"])))) for r in rows1]
ex = vals(ask("EXP", [[t] for t in nexts]))
by = {}
for r, e in zip(rows1, ex):
    m = fb(r["m"]); p = fb(r["mass"]); l = fb(r["ln"])
    key = round(m, 1) if m < 10 else round(m)
    D2 = comp_chain(l, [m, 0.0], lambda t: e)
    c = by.setdefault(key, [0, 0]); c[0] += 1; c[1] += D2 == p
print("x=1 compensated chain, register difference at 64 bits:")
print("  " + "  ".join(f"{k}:{v[1]}/{v[0]}" for k, v in sorted(by.items())))

# (b) series windows
centres = {2: [1.5, 2.2, 2.6, 2.9], 3: [2.1, 2.6, 3.6, 4.4]}
N = 128
rows = [(x, nxt(c, i)) for x in (2, 3) for c in centres[x] for i in range(N)]
mass = vals(ask("POISSON.DIST", [[float(x), m, boolean(False)] for x, m in rows]))
exps = []
for x, m in rows:
    assert in_series(x, m), (x, m)
    exps.append(d(st(-F(S[x]) - deviance_series(x, m))))
ex2 = vals(ask("EXP", [[e] for e in exps]))
res = {}; out = []
for (x, m), p, e, ee in zip(rows, mass, exps, ex2):
    key = (x, [c for c in centres[x] if c <= m < nxt(c, N)][0])
    L1s = d(st(F(ee) / F(sq2pix(x))))
    r = res.setdefault(key, [0, 0]); r[0] += 1; r[1] += L1s == p
    out.append({"x": x, "m": h(m), "mass": h(p), "exponent": h(e), "exp": h(ee), "L1_series": L1s == p})
print("x=2,3 Loader route with the stored series deviance, public EXP, stored division by sqrt(2 pi x):")
for (x, c), (n, ok) in res.items():
    print(f"  x={x} m~{c}: {ok}/{n}")
save("p9_poisson_series_and_fixed_d.json", {"x1_fixed_D": {str(k): v for k, v in by.items()}, "series_rows": out,
                                             "build": meta(ask("POISSON.DIST", [[1.0, 1.0, boolean(False)]]))})
