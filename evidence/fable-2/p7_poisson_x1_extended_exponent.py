"""P7 (offline on p3/p4 answers): x = 1 with the public LN stored and everything after it accurate:
 H6a st(exp(LN m - m))  exact argument, exact exp, RN64 then RN53;  H6b RN53 only;
 H6c argument rounded to 64 bits first (x87 subtraction), then exact exp, RN64 then RN53;  H6d RN53 only.
 H7  the same with LN m replaced by the exact ln m (no stored log): tells whether the stored LN is the error source."""
from common import *
from decimal import Decimal, getcontext
import json
getcontext().prec = 60


def dexp(q):
    return F((Decimal(q.numerator) / Decimal(q.denominator)).exp())


def dln(m):
    return F(Decimal(m).ln())


rows = json.load(open("p3_poisson_x1_windows.json"))["rows"] + json.load(open("p4_poisson_x1_large.json"))["rows"]
by = {}
for r in rows:
    m = fb(r["m"]); p = fb(r["mass"]); l = fb(r["ln"])
    key = round(m, 1) if m < 10 else round(m)
    arg = F(l) - F(m)
    H6a = d(st(dexp(arg))); H6b = d(rn(dexp(arg)))
    arg64 = rn(arg, 64)
    H6c = d(st(dexp(arg64))); H6d = d(rn(dexp(arg64)))
    H7 = d(st(dexp(dln(m) - F(m))))
    c = by.setdefault(key, [0] * 6)
    c[0] += 1
    for i, v in enumerate((H6a, H6b, H6c, H6d, H7)):
        c[i + 1] += v == p
print("window   n  H6a  H6b  H6c  H6d   H7")
for k in sorted(by):
    print(f"{k:7} {by[k][0]:3} " + " ".join(f"{v:4}" for v in by[k][1:]))
