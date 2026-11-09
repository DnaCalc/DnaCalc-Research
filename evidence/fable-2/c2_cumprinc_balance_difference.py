"""C2: CUMPRINC(start..end) as a difference of two discounted balances, each product stored:
   P = st( st(q * u_e) - st(q * u_s) ),  q = st(pv / em_n), u_e = EXP(-(n - end) L), u_s = EXP(-(n - start + 1) L), L = ln(1 + r)
 (the closed form sum_{k=s}^{e} pmt v^(n-k+1) = (pmt/r)(v^(n-e) - v^(n-s+1))). Tested first on the ld1 scan r = 0.75, n = 1000,
 k = 1000 (em_n = -1 exactly, so q = -pv), where every one-rounding form failed; then on pvscan-40 (n = 40, k = 40/20/1) with
 em_n = expm1(-40 L) taken as the correctly rounded value (a stand-in for Excel's internal expm1)."""
from common import *
from decimal import Decimal, getcontext
getcontext().prec = 60
r = 0.75
L = vals(ask("LN", [[1.75]]))[0]
print("L = LN(1.75) =", h(L))


def leaves(ms):
    ex = vals(ask("EXP", [[d(-F(m) * F(L))] if True else None for m in ms]))
    return dict(zip(ms, ex))


def shapes(pv, q, ue, us):
    q, ue, us = F(q), F(ue), F(us)
    return {
        "st(st(q ue) - st(q us))": d(st(st(q * ue) - st(q * us))),
        "st(q ue - q us) ext": d(st(rn(q * ue, 64) - rn(q * us, 64))),
        "st(q * st(ue - us))": d(st(q * st(ue - us))),
        "st(q * (ue - us)) ext": d(st(q * rn(ue - us, 64))),
        "rn53 all": d(rn(rn(q * ue) - rn(q * us))),
    }


# 1. the big loan
rows = [l.split() for l in open("../ld1/pvscan-big1.txt") if l.strip()]
n, k = 1000, 1000
lv = leaves([n - k, n - k + 1])
ue, us = lv[n - k], lv[n - k + 1]
print("u_e = EXP(0) =", ue, " u_s = EXP(-L) =", h(us))
score = {}
for a, b, c in rows:
    pv, P = fb(a), fb(c)
    q = d(st(F(pv) / F(-1)))
    for name, v in shapes(pv, q, ue, us).items():
        score[name] = score.get(name, 0) + (v == P)
print(f"big loan r=0.75 n=1000 k=1000, {len(rows)} rows:")
for name, s in score.items():
    print(f"   {name:28} {s}")

# 2. the n = 40 scans
em = d(rn(F((Decimal(-40) * Decimal(L)).exp()) - 1))
print("em_40 (correctly rounded expm1 stand-in) =", h(em))
for k in (40, 20, 1):
    rows = [l.split() for l in open(f"../ld1/pvscan-{k}.txt") if l.strip()]
    lv = leaves([40 - k, 40 - k + 1])
    ue, us = lv[40 - k], lv[40 - k + 1]
    score = {}
    for a, b, c in rows:
        pv, P = fb(a), fb(c)
        q = d(st(F(pv) / F(em)))
        for name, v in shapes(pv, q, ue, us).items():
            score[name] = score.get(name, 0) + (v == P)
    print(f"loan r=0.75 n=40 k={k}, {len(rows)} rows:")
    for name, s in score.items():
        print(f"   {name:28} {s}")
save("c2_cumprinc_balance_difference.json", {"L": h(L), "em_40_stand_in": h(em), "build": meta(ask("LN", [[1.75]]))})
