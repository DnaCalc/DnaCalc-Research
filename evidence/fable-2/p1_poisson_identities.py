"""P1: is POISSON.DIST(k, m, TRUE) the same bits as CHIDIST(2m, 2k+2) and as 1 - GAMMA.DIST(m, k+1, 1, TRUE) in every
branch of the incomplete gamma (Taylor m < k+1; closed-int k+1 <= m < 31; asymptotic m >= 31; Temme k+1 >= 20)?"""
from common import *
import random
random.seed(2)
rows = []
for k in [1, 2, 3, 5, 8, 12, 18, 19, 25, 40, 80]:
    for m in [0.3, 0.9, 1.7, 2.6, 4.4, 7.9, 12.3, 19.6, 27.5, 30.9, 31.2, 45.0, 77.0, 120.0]:
        for _ in range(3):
            mm = m * (1 + random.uniform(-0.02, 0.02))
            rows.append((k, mm))
pois = vals(ask("POISSON.DIST", [[k, m, boolean(True)] for k, m in rows]))
chi = vals(ask("CHIDIST", [[2 * m, 2 * k + 2] for k, m in rows]))
gam = vals(ask("GAMMA.DIST", [[m, k + 1, 1.0, boolean(True)] for k, m in rows]))
chi2 = vals(ask("CHISQ.DIST.RT", [[2 * m, 2 * k + 2] for k, m in rows]))
by = {}
out = []
for (k, m), p, c, g, c2 in zip(rows, pois, chi, gam, chi2):
    a = k + 1
    reg = "big" if a >= 20 else "taylor" if m < a else "cint" if m < 31 else "asymp"
    ok = isinstance(p, float) and isinstance(c, float) and isinstance(g, float)
    e1 = ok and p == c
    e2 = ok and p == 1 - g
    e3 = ok and c == c2
    d = by.setdefault(reg, [0, 0, 0, 0])
    d[0] += 1; d[1] += e1; d[2] += e2; d[3] += e3
    out.append({"k": k, "m": h(m), "pois": h(p) if isinstance(p, float) else str(p), "chidist": h(c) if isinstance(c, float) else str(c),
                "gamma_p": h(g) if isinstance(g, float) else str(g), "region": reg, "pois_eq_chi": e1, "pois_eq_1_minus_p": e2})
for reg, (n, e1, e2, e3) in by.items():
    print(f"{reg:7} n={n:3} POISSON==CHIDIST {e1:3}  POISSON==1-GAMMA.DIST {e2:3}  CHIDIST==CHISQ.DIST.RT {e3:3}")
save("p1_poisson_identities.json", {"rows": out, "build": meta(ask("POISSON.DIST", [[1.0, 1.0, boolean(True)]]))})
