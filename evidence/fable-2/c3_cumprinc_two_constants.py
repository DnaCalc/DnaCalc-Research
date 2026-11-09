"""C3 (offline): the big loan singleton (r = 0.75, n = 1000, k = 1000; q = -pv) as two-stage products of pv with the loan's own
constants u = EXP(-L) (public, or kept at 64 bits) and r, in every order and store pattern."""
from common import *
from decimal import Decimal, getcontext
getcontext().prec = 60
L = fb("3fe1e85f5e7040d0"); u = fb("3fe2492492492492"); r = F(0.75)
u64 = rn(F((-Decimal(L)).exp()), 64)          # an extended exp of the stored L
u64b = rn(F((Decimal(-0.75).ln().copy_negate()).exp()) if False else F((-Decimal(F(7, 4).numerator) / Decimal(4)).exp()), 64)  # placeholder, unused
rows = [l.split() for l in open("../ld1/pvscan-big1.txt") if l.strip()]
forms = {
    "st(st(q u) r)": lambda q: st(st(q * F(u)) * r),
    "st(st(q r) u)": lambda q: st(st(q * r) * F(u)),
    "st(q st(u r))": lambda q: st(q * st(F(u) * r)),
    "st(q rn64(u r))": lambda q: st(q * rn(F(u) * r, 64)),
    "st(rn64(q u) r)": lambda q: st(rn(q * F(u), 64) * r),
    "st(st(q u64) r)": lambda q: st(st(q * u64) * r),
    "st(q rn64(u64 r))": lambda q: st(q * rn(u64 * r, 64)),
    "st(rn64(q u64) r)": lambda q: st(rn(q * u64, 64) * r),
    "st(q st(u64 r))": lambda q: st(q * st(u64 * r)),
    "st(st(q r) u64)": lambda q: st(st(q * r) * u64),
    "st(q (1-u))": lambda q: st(q * st(1 - F(u))),
    "st(q - st(q u))": lambda q: st(q - st(q * F(u))),
    "st(st(q r)/(1+r))": lambda q: st(st(q * r) / F(7, 4)),
    "st(st(q/(1+r)) r)": lambda q: st(st(q / F(7, 4)) * r),
    "st(q st(r/(1+r)))": lambda q: st(q * st(r / F(7, 4))),
    "st(q rn64(r/(1+r)))": lambda q: st(q * rn(r / F(7, 4), 64)),
    "st(rn64(q r)/(1+r))": lambda q: st(rn(q * r, 64) / F(7, 4)),
    "st(st(q/(1+r)) - q)?": lambda q: st(st(q / F(7, 4)) - q),
    "st(q - st(q/(1+r)))": lambda q: st(q - st(q / F(7, 4))),
}
score = {k: 0 for k in forms}
bias = {k: 0 for k in forms}
for a, b, c in rows:
    pv, P = fb(a), fb(c)
    q = F(-pv)
    for k, f in forms.items():
        v = d(f(q))
        score[k] += v == P
        bias[k] += ulps(P, v)
print(len(rows), "rows")
for k in sorted(score, key=lambda k: -score[k]):
    print(f"  {k:24} {score[k]:5}   mean ulps(model - excel) {bias[k] / len(rows):+.2f}")
