"""P15: the cumulative at a = k + 1 < 20, m >= 31 as the public mass AT k times the finite backward series
 1 + k/m + k(k-1)/m^2 + ... (the NSWC asymptotic form r (a/x) (1 + summ) has r a/x = mass(k, m)):
 W1 st(mass_k (1 + summ_nswc))  W2 st(mass_k h), h = Horner 1 + k/m (1 + (k-1)/m (...)) in binary64   W3 Horner with PC64 stores
 W4 one rounding of mass_k h    W5 forward binary64 sum of the terms   W6 forward with PC64 stores   W7 W1 with summ in PC64.
Also Excel's error against the exact Q in ulps."""
from common import *
from decimal import Decimal, getcontext
import json, math, collections, statistics
getcontext().prec = 80


def exact_q(a, x):
    X = Decimal(x); s = Decimal(0); t = Decimal(1)
    for j in range(a):
        s += t; t = t * X / (j + 1)
    return F(s * (-X).exp())


def err_ulps(ans, truth):
    e = math.frexp(ans)[1]
    return float((F(ans) - truth) / F(2) ** (e - 53))


def nswc_summ(a, x):
    amn = a - 1.0; t = amn / x; wk = [0.0] * 21; wk[1] = t; n = 20; broke = False
    for n_ in range(2, 21):
        amn -= 1.0; t *= amn / x
        if abs(t) <= 1e-3:
            n = n_; broke = True; break
        wk[n_] = t
    if not broke:
        n = 20
    summ = t
    while abs(t) > 5e-15:
        amn -= 1.0; t *= amn / x; summ += t
    for _ in range(n - 1):
        n -= 1; summ += wk[n]
    return summ


def variants(k, m, mk):
    out = {}
    s = nswc_summ(k + 1, m)
    out["W1"] = d(st(F(mk) * F(1.0 + s)))
    hh = 1.0
    for j in range(k, 0, -1):
        hh = 1.0 + hh * j / m
    out["W2"] = d(st(F(mk) * F(hh)))
    out["W4"] = d(rn(F(mk) * F(hh)))
    hx = F(1)
    for j in range(k, 0, -1):
        hx = st(1 + st(st(hx * j) / F(m)))
    out["W3"] = d(st(F(mk) * hx))
    t = 1.0; s5 = 1.0
    for i in range(1, k + 1):
        t = t * (k - i + 1) / m; s5 += t
    out["W5"] = d(st(F(mk) * F(s5)))
    t = F(1); s6 = F(1)
    for i in range(1, k + 1):
        t = st(st(t * (k - i + 1)) / F(m)); s6 = st(s6 + t)
    out["W6"] = d(st(F(mk) * s6))
    return out


p13 = json.load(open("p13_poisson_asymptotic_fresh.json"))
fresh = [(r["k"], fb(r["m"]), fb(r["cum"])) for r in p13["fresh"] if r["cum"][0] in "0123456789abcdef" and len(r["cum"]) == 16]
mk = vals(ask("POISSON.DIST", [[float(k), m, boolean(False)] for k, m, a in fresh]))
by = collections.defaultdict(collections.Counter); tot = collections.Counter(); errs = collections.defaultdict(list)
for (k, m, a), mass in zip(fresh, mk):
    if not isinstance(mass, float):
        continue
    reg = ("k<=3" if k <= 3 else "k<=8" if k <= 8 else "k<=18") + (" m<100" if m < 100 else " m<400" if m < 400 else " m<700")
    tot[reg] += 1
    errs[reg].append(err_ulps(a, exact_q(k + 1, m)))
    for name, v in variants(k, m, mass).items():
        by[reg][name] += v == a
print("region          n   W1   W2   W3   W4   W5   W6 | excel err mean  sd   max|.|")
for reg in sorted(tot):
    e = errs[reg]
    print(f"{reg:14} {tot[reg]:3} " + " ".join(f"{by[reg][n]:4}" for n in ("W1", "W2", "W3", "W4", "W5", "W6")) +
          f" | {statistics.mean(e):+6.2f} {statistics.pstdev(e):5.2f} {max(abs(v) for v in e):6.2f}")
rows = json.load(open("w111_poisson_rows.json"))
w = [(math.floor(x), m, fb(bits[2:])) for src, x, m, cum, bits, rid in rows if cum and bits.startswith("0x") and 1 <= math.floor(x) <= 18 and m >= 31]
wmk = vals(ask("POISSON.DIST", [[float(k), m, boolean(False)] for k, m, a in w]))
sc = collections.Counter(); n = 0
for (k, m, a), mass in zip(w, wmk):
    if not isinstance(mass, float) or mass == 0.0:
        continue
    n += 1
    for name, v in variants(k, m, mass).items():
        sc[name] += v == a
print(f"W111 rows (nonzero public mass at k) {n}: " + "  ".join(f"{k} {v}" for k, v in sorted(sc.items())))
save("p15_poisson_cumulative_from_mass_k.json", {"fresh_mass_k": [h(v) if isinstance(v, float) else str(v) for v in mk],
                                                  "w111_mass_k": [h(v) if isinstance(v, float) else str(v) for v in wmk], "build": meta(ask("LN", [[2.0]]))})
