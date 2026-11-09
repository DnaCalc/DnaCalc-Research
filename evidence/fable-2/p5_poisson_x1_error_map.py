"""P5 (offline on p3/p4 answers): Excel's error at x = 1 in units of the answer's ulp against the exact m e^-m, per window."""
from common import *
from decimal import Decimal, getcontext
import json, math, statistics
getcontext().prec = 60


def exact_mass(m):
    v = Decimal(m)
    return F(v * (-v).exp())


def err_ulps(answer, truth):
    e = math.frexp(answer)[1]
    ulp = F(2) ** (e - 53)
    return float((F(answer) - truth) / ulp)


rows = json.load(open("p3_poisson_x1_windows.json"))["rows"] + json.load(open("p4_poisson_x1_large.json"))["rows"]
by = {}
for r in rows:
    m = fb(r["m"]); p = fb(r["mass"])
    by.setdefault(round(m, 1) if m < 10 else round(m), []).append((m, err_ulps(p, exact_mass(m))))
print("window   n   mean    sd     min    max   | first 24 errors")
for c in sorted(by):
    es = [e for m, e in by[c]]
    print(f"{c:7} {len(es):3} {statistics.mean(es):+6.2f} {statistics.pstdev(es):5.2f} {min(es):+6.2f} {max(es):+6.2f} |",
          " ".join(f"{e:+.2f}" for e in es[:24]))
