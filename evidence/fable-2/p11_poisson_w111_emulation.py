"""P11: emulate calc-func's POISSON.DIST on the W111 bank rows (Excel 20430) that the OxFunc cache lacks: masses at x >= 4
(Loader route; Stirling error from the table to 6 and the ten-term 64-bit series from 7), and the cumulative's Taylor,
closed-int and asymptotic branches of the incomplete gamma (a = k + 1 < 20). Temme (a >= 20) is not emulated here.
Leaves from Excel today: LN(st(x/m)), EXP(e), POISSON.DIST(k+1, m, FALSE) as the gratio prefactor, EXP(-m)."""
from common import *
from p9_poisson_series_and_fixed_d import deviance_series, in_series, sq2pix, r64
import json, math, collections
STIRLING_TABLE = [0.0, 0.0810614667953272, 0.0413406959554092, 0.0276779256849983, 0.0207906721037650, 0.0166446911898211, 0.0138761288230707]
TERMS = [F(1, 12), F(-1, 360), F(1, 1260), F(-1, 1680), F(1, 1188), F(-691, 360360), F(1, 156), F(-3617, 122400), F(43867, 244188), F(-174611, 125400)]
TERMS = [F(float(t)) for t in TERMS]     # binary64 literals, as Rust writes them


def stirling_series(n):
    t = r64(F(1) / r64(F(n) * F(n)))
    w = TERMS[9]
    for term in TERMS[:9][::-1]:
        w = r64(r64(w * t) + term)
    return st(w / F(n))


def stirling_error(n):
    return F(STIRLING_TABLE[int(n)]) if n < 7 else stirling_series(n)


rows = json.load(open("w111_poisson_rows.json"))
rows = [(src, x, m, cum, fb(bits[2:]), rid) for src, x, m, cum, bits, rid in rows if bits.startswith("0x")]
# --- masses, x >= 4, m >= 1.1 (Loader) and x <= 3 for completeness through the same code
mass_rows = [(i, x, m) for i, (src, x, m, cum, ans, rid) in enumerate(rows) if not cum and x >= 2 and m >= 1.1 and x == int(x)]
ratios = sorted({d(st(F(x) / F(m))) for i, x, m in mass_rows if not in_series(x, m)})
LNR = dict(zip(ratios, vals(ask("LN", [[q] for q in ratios]))))
exps = {}
for i, x, m in mass_rows:
    dev = deviance_series(x, m) if in_series(x, m) else st(st(F(x) * F(LNR[d(st(F(x) / F(m)))])) + st(F(m) - F(x)))
    exps[i] = d(st(-stirling_error(x) - dev))
keys = sorted(set(exps.values()))
EXPE = dict(zip(keys, vals(ask("EXP", [[e] for e in keys]))))
pred = {}
for i, x, m in mass_rows:
    pred[i] = d(st(F(EXPE[exps[i]]) / F(sq2pix(x))))
# --- cumulative, k >= 1, a = k + 1 < 20
cum_rows = [(i, math.floor(x), m) for i, (src, x, m, cum, ans, rid) in enumerate(rows) if cum and 1 <= math.floor(x) and math.floor(x) + 1 < 20 and m > 0]
masses_needed = sorted({(k + 1, m) for i, k, m in cum_rows if not (k + 1 <= m < 31)})
PM = dict(zip(masses_needed, vals(ask("POISSON.DIST", [[float(a), m, boolean(False)] for a, m in masses_needed]))))
expneg = sorted({m for i, k, m in cum_rows if k + 1 <= m < 31})
EXN = dict(zip(expneg, vals(ask("EXP", [[-m] for m in expneg]))))
branch = {}
for i, k, m in cum_rows:
    a = k + 1; x = m
    if a <= x < 31:
        # statement 210: integer a finite sum in binary64
        summ = EXN[m]; t = summ; c = 0.0
        for n in range(1, a):
            c += 1.0
            t = (x * t) / c
            summ += t
        pred[i] = summ; branch[i] = "closed-int"
    else:
        r_series = PM[(a, m)]
        if not isinstance(r_series, float):
            continue
        if x <= max(a, 2.30258509299404590109):
            rr = float(a); c = 1.0; ans = 1.0
            while True:
                rr += 1.0
                c = d(st(F(c) * st(F(x) / F(rr))))
                ans = d(st(F(ans) + F(c)))
                if d(st(F(c) / F(ans))) <= 1.1102230246251565e-16:
                    break
            p = d(st(F(r_series) * F(ans)))
            pred[i] = 0.5 + (0.5 - p); branch[i] = "taylor"
        elif x < 31:
            branch[i] = "cf"      # not emulated
        else:
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
            qans = d(st(st(F(r_series) * st(F(a) / F(x))) * F(1.0 + summ)))
            pred[i] = qans; branch[i] = "asymp"
tot = collections.Counter(); bad = collections.Counter(); fails = []
for i, (src, x, m, cum, ans, rid) in enumerate(rows):
    if i not in pred:
        continue
    cls = ("mass x<=3" if x <= 3 else "mass x<=6" if x <= 6 else "mass x<=20" if x <= 20 else "mass x>20") if not cum else "cum " + branch[i]
    tot[cls] += 1
    if pred[i] != ans:
        bad[cls] += 1
        fails.append({"class": cls, "x": x, "m": h(m), "excel": h(ans), "model": h(pred[i]), "ulps": ulps(ans, pred[i]) if pred[i] == pred[i] else None, "id": rid, "src": src})
print("class         rows fails")
for cls in sorted(tot):
    print(f"{cls:13} {tot[cls]:5} {bad[cls]:5}")
print("ulps histogram of fails:", sorted(collections.Counter(f["ulps"] for f in fails).items(), key=lambda t: (t[0] is None, t[0]))[:30])
for f in fails[:12]:
    print("  ", f)
save("p11_poisson_w111_emulation.json", {"fails": fails, "totals": dict(tot), "bad": dict(bad), "build": meta(ask("LN", [[2.0]]))})
