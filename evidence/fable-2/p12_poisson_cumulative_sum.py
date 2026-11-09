"""P12: POISSON.DIST(k, m, TRUE) for 1 <= k <= 18 on the W111 rows (means mostly 700-760, results near the subnormal
boundary) and on fresh ordinary rows: is the cumulative the finite sum e^-m sum_{j<=k} m^j/j! for every mean (integer shape
parameter), rather than the asymptotic expansion from m >= 31? Variants of the sum's arithmetic:
 S1 binary64 forward: t = (m t)/c; s += t          S2 the same in 64-bit registers, stored once
 S3 binary64 forward with t = t (m/c)              S4 every operation stored from PC64: t = st(st(m t)/c), s = st(s + t)
 S5 Horner from the top: s = 1; for j = k..1: s = 1 + s m/j; Q = E s (binary64)   S6 S5 with PC64 stores
 E = the public EXP(-m)."""
from common import *
import json, math, collections, random
rows = json.load(open("w111_poisson_rows.json"))
cum = [(math.floor(x), m, fb(bits[2:])) for src, x, m, cum_, bits, rid in rows if cum_ and bits.startswith("0x") and 1 <= math.floor(x) <= 18 and m > 0]
random.seed(8)
fresh = []
for k in [1, 2, 3, 5, 8, 12, 18]:
    for m in [0.7, 2.5, 6.0, 12.0, 25.0, 35.0, 60.0, 120.0, 300.0, 600.0]:
        for _ in range(4):
            fresh.append((k, m * (1 + random.uniform(-0.05, 0.05))))
fresh_ans = vals(ask("POISSON.DIST", [[float(k), m, boolean(True)] for k, m in fresh]))
fresh_rows = [(k, m, a) for (k, m), a in zip(fresh, fresh_ans) if isinstance(a, float)]
allrows = [("w111", k, m, a) for k, m, a in cum] + [("fresh", k, m, a) for k, m, a in fresh_rows]
ms = sorted({m for _, k, m, a in allrows})
EXN = dict(zip(ms, vals(ask("EXP", [[-m] for m in ms]))))


def r64(q):
    return rn(q, 64)


def flush(v):
    v = float(v)
    return 0.0 if (v != 0.0 and abs(v) < 2.2250738585072014e-308) else v


def variants(k, m, E):
    out = {}
    # S1
    s = E; t = E; c = 0.0
    for n in range(1, k + 1):
        c += 1.0; t = (m * t) / c; s += t
    out["S1"] = s
    # S2
    s = F(E); t = F(E)
    for n in range(1, k + 1):
        t = r64(r64(F(m) * t) / n); s = r64(s + t)
    out["S2"] = d(rn(s))
    # S3
    s = E; t = E; c = 0.0
    for n in range(1, k + 1):
        c += 1.0; t = t * (m / c); s += t
    out["S3"] = s
    # S4
    s = F(E); t = F(E)
    for n in range(1, k + 1):
        t = st(st(F(m) * t) / n); s = st(s + t)
    out["S4"] = d(s)
    # S5
    h_ = 1.0
    for j in range(k, 0, -1):
        h_ = 1.0 + h_ * m / j
    out["S5"] = E * h_
    # S6
    h_ = F(1)
    for j in range(k, 0, -1):
        h_ = st(1 + st(st(h_ * F(m)) / j))
    out["S6"] = d(st(F(E) * h_))
    return out


score = collections.defaultdict(collections.Counter); tot = collections.Counter()
for src, k, m, a in allrows:
    E = EXN[m]
    if not isinstance(E, float):
        continue
    reg = src + (" m<a" if m < k + 1 else " m<31" if m < 31 else " m<700" if m < 700 else " m>=700")
    tot[reg] += 1
    for name, v in variants(k, m, E).items():
        score[reg][name] += flush(v) == a
print("region            rows   S1   S2   S3   S4   S5   S6")
for reg in sorted(tot):
    print(f"{reg:17} {tot[reg]:5} " + " ".join(f"{score[reg][n]:4}" for n in ("S1", "S2", "S3", "S4", "S5", "S6")))
save("p12_poisson_cumulative_sum.json", {"fresh": [{"k": k, "m": h(m), "poisson": h(a)} for k, m, a in fresh_rows],
                                         "exp_neg_m": {h(m): h(e) if isinstance(e, float) else str(e) for m, e in EXN.items()},
                                         "build": meta(ask("POISSON.DIST", [[1.0, 1.0, boolean(True)]]))})
