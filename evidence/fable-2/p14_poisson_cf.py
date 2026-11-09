"""P14 (offline on p13 answers): the cumulative at a = k + 1 < 20, m >= 31 as the upper tail's continued fraction (Cephes igamc
shape as calc-func keeps it for m < 31: 64-bit registers, stored once), Q = st(st(mass a) * CF) with Excel's public mass."""
from common import *
import json, collections


def r64(q):
    return rn(q, 64)


def gamma_fraction(a, x):
    BIG = 4503599627370496.0; BIGINV = 2.220446049250313e-16
    y = 1.0 - a; c = 0.0; z = x + y + 1.0
    pkm2, qkm2 = F(1), F(x)
    pkm1 = r64(F(x) + 1)
    qkm1 = r64(F(z) * F(x))
    ans = r64(pkm1 / qkm1)
    while True:
        c += 1.0; y += 1.0; z += 2.0
        yc = r64(F(y) * F(c))
        pk = r64(r64(pkm1 * F(z)) - r64(pkm2 * yc))
        qk = r64(r64(qkm1 * F(z)) - r64(qkm2 * yc))
        t = 1.0
        if qk != 0:
            nx = r64(pk / qk)
            t = abs(d(st(r64(ans - nx) / nx)))
            ans = nx
        pkm2, qkm2 = pkm1, qkm1
        pkm1, qkm1 = pk, qk
        if abs(d(st(pk))) > BIG:
            pkm2, pkm1 = r64(pkm2 * F(BIGINV)), r64(pkm1 * F(BIGINV))
            qkm2, qkm1 = r64(qkm2 * F(BIGINV)), r64(qkm1 * F(BIGINV))
        if t <= 1.1102230246251565e-16:
            return st(ans)


def variants(a, x, mass):
    cf = gamma_fraction(a, x)
    r = st(F(mass) * a)
    return {"CF st(r*cf)": d(st(r * cf)), "CF st(mass*a*cf) regs": d(st(r64(F(mass) * a) * cf)), "CF st(mass*st(a*cf))": d(st(F(mass) * st(a * cf)))}


p13 = json.load(open("p13_poisson_asymptotic_fresh.json"))
sc = collections.Counter(); n = 0; by = collections.defaultdict(collections.Counter); tot = collections.Counter()
for r in p13["fresh"]:
    if not r["cum"].startswith("3") and not r["cum"][0] in "0123": pass
    try:
        a_ = fb(r["cum"]); mass = fb(r["mass_k1"])
    except Exception:
        continue
    k = r["k"]; m = fb(r["m"])
    reg = "m<100" if m < 100 else "m<400" if m < 400 else "m<700"
    tot[reg] += 1; n += 1
    for name, v in variants(k + 1, m, mass).items():
        by[reg][name] += v == a_
print(f"fresh rows {n}")
for reg in sorted(tot):
    print(f"  {reg:6} {tot[reg]:4}: " + "  ".join(f"{k} {v}" for k, v in sorted(by[reg].items())))
rows = json.load(open("w111_poisson_rows.json"))
import math
w = [(math.floor(x), m, fb(bits[2:])) for src, x, m, cum, bits, rid in rows if cum and bits.startswith("0x") and 1 <= math.floor(x) <= 18 and m >= 31]
sc2 = collections.Counter(); n2 = 0
for (k, m, a_), mh in zip(w, p13["w111_mass_k1"]):
    try:
        mass = fb(mh)
    except Exception:
        continue
    if mass == 0.0:
        continue
    n2 += 1
    for name, v in variants(k + 1, m, mass).items():
        sc2[name] += v == a_
print(f"W111 rows with a nonzero public mass {n2}: " + "  ".join(f"{k} {v}" for k, v in sorted(sc2.items())))
