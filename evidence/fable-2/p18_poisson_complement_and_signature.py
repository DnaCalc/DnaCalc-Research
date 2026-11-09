"""P18 (offline on p16 answers): for a <= m, (i) is the cumulative the complement of the Taylor series continued past m = a
(T1 Q = 0.5 + (0.5 - st(M ans)); T2 Q = st(1 - M ans) with the product in the registers; T3 Q = st(1 - st(M ans)));
(ii) the error signature: Excel's error in ulps against the exact Q versus m and versus Q's size."""
from common import *
import json, math, statistics
from p15_poisson_cumulative_from_mass_k import exact_q, err_ulps
from p16_poisson_cumulative_prefactor import taylor_ans


def taylor_ans_regs(a, x):
    rr = F(a); c = F(1); ans = F(1)
    while True:
        rr += 1
        c = rn(c * rn(F(x) / rr, 64), 64); ans = rn(ans + c, 64)
        if d(st(c / ans)) <= 1.1102230246251565e-16:
            return ans


rows = json.load(open("p16_poisson_cumulative_prefactor.json"))["rows"]
sc = {"T1": 0, "T2": 0, "T3": 0, "T4": 0}; n = 0
sig = []
for r in rows:
    k = r["k"]; a = k + 1; m = fb(r["m"]); q = fb(r["cum"]); M = F(fb(r["mass_a"]))
    if m < a or m >= 31:
        continue
    n += 1
    ans = taylor_ans(a, m)
    T1 = 0.5 + (0.5 - d(st(M * F(ans))))
    T2 = d(st(1 - rn(M * F(ans), 64)))
    T3 = d(st(1 - st(M * F(ans))))
    T4 = d(st(1 - rn(M * taylor_ans_regs(a, m), 64)))
    for name, v in (("T1", T1), ("T2", T2), ("T3", T3), ("T4", T4)):
        sc[name] += v == q
    sig.append((m, q, err_ulps(q, exact_q(a, m)), k))
print(f"a<=m<31 rows {n}: " + "  ".join(f"{k} {v}" for k, v in sc.items()))
sig.sort(key=lambda t: t[1])
print("error vs size of Q (a<=m<31): bins of -log10 Q")
import collections
bins = collections.defaultdict(list)
for m, q, e, k in sig:
    bins[int(-math.log10(q))].append(e)
for b in sorted(bins):
    e = bins[b]
    print(f"   Q~1e-{b:<2} n={len(e):3} mean {statistics.mean(e):+6.2f} sd {statistics.pstdev(e):5.2f} max|.| {max(abs(v) for v in e):5.2f}")
bins = collections.defaultdict(list)
for m, q, e, k in sig:
    bins[k].append(e)
print("error by k:", {k: f"sd {statistics.pstdev(v):.2f} n={len(v)}" for k, v in sorted(bins.items())})
