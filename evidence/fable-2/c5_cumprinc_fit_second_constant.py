"""C5 (offline on c4 and the ld1 big-loan scan): with q = -pv exactly, fit P = st(T * K) where T is a known pv-dependent
rounding, T = st(pv r) (the payment) or T = st(pv / (1 + r)), and K a free real: intersect K's preimage intervals over the
rows. A non-empty intersection names the shape and pins K; then compare K with the loan's natural constants."""
from common import *
import json, struct


def pre(v):
    """Reals whose RN53 is the double v (closed interval); the RN64-then-RN53 preimage differs only at the endpoints."""
    v = abs(v)
    lo = nxt(v, -1); hi = nxt(v, 1)
    return (F(v) + F(lo)) / 2, (F(v) + F(hi)) / 2


def fit(samples, inner):
    cur = None
    for pv, P in samples:
        T = inner(F(pv))
        lo, hi = pre(P)
        iv = (lo / T, hi / T)
        cur = iv if cur is None else (max(cur[0], iv[0]), min(cur[1], iv[1]))
        if cur[0] > cur[1]:
            return None
    return cur


def rel(K, C):
    return float((K - C) / C / F(2) ** -53)


data = json.load(open("c4_cumprinc_last_period_bias.json"))
u = {float(k): fb(v) for k, v in data["exp_neg_ln"].items()}
by = {}
for r in data["rows"]:
    by.setdefault(r["rate"], []).append((fb(r["pv"]), abs(fb(r["cumprinc"]))))
big = [(fb(a), abs(fb(c))) for a, b, c in (l.split() for l in open("../ld1/pvscan-big1.txt") if l.strip())]
by[0.75] = by[0.75] + big
for r, samples in sorted(by.items()):
    C = F(r) / (1 + F(r))
    print(f"rate {r} ({len(samples)} rows)")
    for name, inner, Cexp in (("T = st(pv r)      K ~ 1/(1+r)", lambda pv: st(pv * F(r)), 1 / (1 + F(r))),
                              ("T = st(pv/(1+r))  K ~ r", lambda pv: st(pv / (1 + F(r))), F(r)),
                              ("T = pv            K ~ r/(1+r)", lambda pv: pv, C)):
        iv = fit(samples, inner)
        if iv is None:
            print(f"   {name}: empty")
        else:
            lo, hi = iv
            print(f"   {name}: K in [{float(lo)!r}, {float(hi)!r}]  rel. to exact: [{rel(lo, Cexp):+.3f}, {rel(hi, Cexp):+.3f}] units of 2^-53;"
                  f" RN53 {rel(rn(Cexp), Cexp):+.3f}, RN64 {rel(rn(Cexp, 64), Cexp):+.3f}, EXP(-LN(1+r)) {rel(F(u[r]), 1 / (1 + F(r))):+.3f}")
