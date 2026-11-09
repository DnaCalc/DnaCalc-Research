"""P2: which other Excel functions share the Poisson mass's arithmetic? POISSON.DIST(x, m, FALSE) against
GAMMA.DIST(m, x+1, 1, FALSE) (the same quantity m^x e^-m / x!), and against stored products of GAMMA.DIST(m, x, 1, FALSE)
(mass at x-1) with m/x; and for x = 1 against m * EXP(-m)."""
from common import *
import random
random.seed(3)
rows = []
for x in [1, 2, 3, 5, 8, 13, 21, 34, 60, 100]:
    for m in [0.05, 0.4, 0.9, 1.05, 1.15, 1.7, 2.6, 4.4, 7.9, 12.3, 19.6, 30.9, 45.0, 77.0, 120.0, 300.0]:
        for _ in range(3):
            rows.append((x, m * (1 + random.uniform(-0.02, 0.02))))
mass = vals(ask("POISSON.DIST", [[x, m, boolean(False)] for x, m in rows]))
pdf1 = vals(ask("GAMMA.DIST", [[m, x + 1, 1.0, boolean(False)] for x, m in rows]))
pdf0 = vals(ask("GAMMA.DIST", [[m, x, 1.0, boolean(False)] for x, m in rows]))
ex = vals(ask("EXP", [[-m] for x, m in rows]))
by = {}
out = []
for (x, m), p, g1, g0, e in zip(rows, mass, pdf1, pdf0, ex):
    reg = ("m<1.1" if m < 1.1 else "m<x" if m < x else "m<3x" if m < 3 * x else "m>=3x") + ("/x=1" if x == 1 else "/x>1")
    ok = all(isinstance(v, float) for v in (p, g1, g0, e))
    e1 = ok and p == g1
    # mass(x) = mass(x-1) * m / x ; GAMMA pdf alpha=x at m is m^(x-1) e^-m/(x-1)! = mass(x-1)
    e2 = ok and p == d(st(st(F(g0) * F(m)) / F(x)))
    e3 = ok and p == d(st(F(g0) * st(F(m) / F(x))))
    e4 = ok and x == 1 and p == d(st(F(m) * F(e)))
    e5 = ok and x == 1 and p == d(rn(F(m) * F(e)))
    dd = by.setdefault(reg, [0, 0, 0, 0, 0, 0])
    dd[0] += 1; dd[1] += e1; dd[2] += e2; dd[3] += e3; dd[4] += e4; dd[5] += e5
    out.append({"x": x, "m": h(m), "mass": h(p) if isinstance(p, float) else str(p), "gamma_pdf_x1": h(g1) if isinstance(g1, float) else str(g1),
                "gamma_pdf_x": h(g0) if isinstance(g0, float) else str(g0), "exp_neg_m": h(e) if isinstance(e, float) else str(e), "ulps_mass_vs_pdf": ulps(p, g1) if ok else None})
print("region        n  mass==GAMMA.DIST(m,x+1,1,F)  ==st(st(pdf0*m)/x)  ==st(pdf0*st(m/x))  x=1:==st(m*EXP)  ==rn(m*EXP)")
for reg, (n, e1, e2, e3, e4, e5) in sorted(by.items()):
    print(f"{reg:12} {n:3}  {e1:3}  {e2:3}  {e3:3}  {e4:3}  {e5:3}")
save("p2_poisson_mass_oracles.json", {"rows": out, "build": meta(ask("POISSON.DIST", [[1.0, 1.0, boolean(False)]]))})
