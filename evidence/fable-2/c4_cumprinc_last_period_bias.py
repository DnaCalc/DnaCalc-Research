"""C4: the last period's principal CUMPRINC(r, n, pv, n, n, 0) at several rates with n large enough that expm1(-n L) is exactly -1
(so q = -pv and the answer is pv times a rate constant r/(1+r) in exact arithmetic). For 300 consecutive pv per rate: Excel's
error against the exact pv r/(1+r) in ulps (mean, sd), and how many rows equal one x87 rounding st(pv * C) for C = RN64(r/(1+r)),
C = RN53(r/(1+r)), and the public st(pv * st(EXP(-LN(1+r)) * r)). Rates where 1 + r is a power of two are the controls."""
from common import *
import math, statistics
rates = [0.75, 0.5, 0.25, 0.1, 0.05, 1.0, 3.0, 1.5, 0.3]
n = 1000
pvs = [nxt(3.3, i) for i in range(300)]
rows = [(r, pv) for r in rates for pv in pvs]
ans = vals(ask("CUMPRINC", [[r, float(n), pv, float(n), float(n), 0.0] for r, pv in rows]))
ln1r = dict(zip(rates, vals(ask("LN", [[1.0 + r] for r in rates]))))
u = dict(zip(rates, vals(ask("EXP", [[-ln1r[r]] for r in rates]))))
pmt = vals(ask("PMT", [[r, float(n), pv, 0.0, 0.0] for r, pv in rows]))
out = []
print("rate   mean_err  sd    st(pv*RN64(C)) st(pv*RN53(C)) st(pv*st(u r))  pmt==-r*pv  |err| of C_eff (units of 2^-53 rel)")
for r in rates:
    sub = [(pv, a, p) for (rr, pv), a, p in zip(rows, ans, pmt) if rr == r and isinstance(a, float)]
    C = F(r) / (1 + F(r))
    errs = [float((F(a) + F(pv) * C) / (F(2) ** (math.frexp(a)[1] - 53))) for pv, a, p in sub]
    m1 = sum(d(st(F(pv) * rn(C, 64))) == -a for pv, a, p in sub)
    m2 = sum(d(st(F(pv) * rn(C, 53))) == -a for pv, a, p in sub)
    m3 = sum(d(st(F(pv) * st(F(u[r]) * F(r)))) == -a for pv, a, p in sub)
    pm = sum(p == -r * pv for pv, a, p in sub)
    ceff = statistics.mean(float((-F(a) / F(pv) - C) / C) for pv, a, p in sub) / 2 ** -53
    print(f"{r:5} {statistics.mean(errs):+7.2f} {statistics.pstdev(errs):5.2f}   {m1:4}/{len(sub)}      {m2:4}/{len(sub)}      {m3:4}/{len(sub)}      {pm:4}      {ceff:+.2f}")
    out += [{"rate": r, "pv": h(pv), "cumprinc": h(a), "pmt": h(p) if isinstance(p, float) else str(p)} for pv, a, p in sub]
save("c4_cumprinc_last_period_bias.json", {"rows": out, "ln_1_plus_r": {str(r): h(v) for r, v in ln1r.items()}, "exp_neg_ln": {str(r): h(v) for r, v in u.items()},
                                           "build": meta(ask("LN", [[2.0]]))})
