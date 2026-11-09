"""P13: the cumulative's asymptotic branch (a = k + 1 < 20, m >= 31) on fresh ordinary rows, with Excel's own POISSON.DIST(k+1, m,
FALSE) as the prefactor: current chain qans = st(st(r st(a/m)) (1 + summ)), and variants: V2 one store of r (a/m) (1+summ) from
the registers; V3 st(r st(st(a/m) (1+summ))); V4 st(st(r/m) a (1+summ)) stores each. Also the W111 rows at m >= 700 with the
product chain's stores denormalized (binary64 subnormal grid) where they fall below DBL_MIN."""
from common import *
import json, math, random, collections
random.seed(13)


def st_sub(q):
    """x87 store with the binary64 subnormal grid: below 2^-1022 the result is a multiple of 2^-1074."""
    v = rn(q, 64)
    if abs(v) < F(2) ** -1022:
        g = F(2) ** -1074
        n_ = v / g
        k = n_.numerator // n_.denominator
        rem = n_ - k
        if rem > F(1, 2) or (rem == F(1, 2) and k % 2 == 1):
            k += 1
        return k * g
    return rn(v, 53)


def asym(a, x, r, store):
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
    one = F(1.0 + summ)
    return {
        "V1": d(store(store(F(r) * store(F(a) / F(x))) * one)),
        "V2": d(store(rn(F(r) * rn(F(a) / F(x), 64), 64) * one)),
        "V3": d(store(F(r) * store(store(F(a) / F(x)) * one))),
        "V4": d(store(store(store(F(r) / F(x)) * a) * one)),
        "V5": d(store(F(r) * rn(rn(F(a) / F(x), 64) * one, 64))),
    }


fresh = [(k, m * (1 + random.uniform(-0.05, 0.05))) for k in [1, 2, 3, 5, 8, 12, 18] for m in [35.0, 60.0, 120.0, 300.0, 600.0, 690.0] for _ in range(6)]
fa = vals(ask("POISSON.DIST", [[float(k), m, boolean(True)] for k, m in fresh]))
fm = vals(ask("POISSON.DIST", [[float(k + 1), m, boolean(False)] for k, m in fresh]))
sc = collections.Counter(); n = 0
for (k, m), a, r in zip(fresh, fa, fm):
    if not isinstance(a, float) or not isinstance(r, float):
        continue
    n += 1
    for name, v in asym(k + 1, m, r, st).items():
        sc[name] += v == a
print(f"fresh asymptotic rows {n}: " + "  ".join(f"{k} {v}" for k, v in sorted(sc.items())))
rows = json.load(open("w111_poisson_rows.json"))
w = [(math.floor(x), m, fb(bits[2:])) for src, x, m, cum, bits, rid in rows if cum and bits.startswith("0x") and 1 <= math.floor(x) <= 18 and m >= 31]
wm = vals(ask("POISSON.DIST", [[float(k + 1), m, boolean(False)] for k, m, a in w]))
sc2 = collections.Counter(); sc3 = collections.Counter(); n2 = 0; zero_mass = 0
for (k, m, a), r in zip(w, wm):
    if not isinstance(r, float):
        continue
    if r == 0.0:
        zero_mass += 1
    n2 += 1
    for name, v in asym(k + 1, m, r, st).items():
        sc2[name] += v == a
    for name, v in asym(k + 1, m, r, st_sub).items():
        sc3[name] += v == a
print(f"W111 asymptotic rows {n2} (public mass zero on {zero_mass}): plain stores " + "  ".join(f"{k} {v}" for k, v in sorted(sc2.items())))
print(f"   subnormal-grid stores " + "  ".join(f"{k} {v}" for k, v in sorted(sc3.items())))
save("p13_poisson_asymptotic_fresh.json", {"fresh": [{"k": k, "m": h(m), "cum": h(a) if isinstance(a, float) else str(a), "mass_k1": h(r) if isinstance(r, float) else str(r)} for (k, m), a, r in zip(fresh, fa, fm)],
                                           "w111_mass_k1": [h(r) if isinstance(r, float) else str(r) for r in wm], "build": meta(ask("LN", [[2.0]]))})
