"""P3: POISSON.DIST(1, m, FALSE) = m e^-m on windows of 128 adjacent m. Candidate programs with public LN and EXP leaves:
 A  st(m * EXP(-m))                      B  rn(m * EXP(-m))
 C  EXP(st(LN m - m))                    D  compensated: EXP(next) * st(1 + dropped), next = st(LN m - m), dropped = st(st(LN m - next) - m)
 E  st(EXP(st(LN m - m)) * st(1 + dropped)) with dropped = st(LN m - next) - m in binary64 (not stored twice)
"""
from common import *
centres = [0.3, 0.9, 1.05, 3.3, 40.0, 300.0]
N = 128
ms = [nxt(c, i) for c in centres for i in range(N)]
mass = vals(ask("POISSON.DIST", [[1.0, m, boolean(False)] for m in ms]))
ln = vals(ask("LN", [[m] for m in ms]))
exn = vals(ask("EXP", [[-m] for m in ms]))
nexts = [d(st(F(l) - F(m))) for l, m in zip(ln, ms)]
ex_next = vals(ask("EXP", [[t] for t in nexts]))
res = {}
rows = []
for m, p, l, e, t, et in zip(ms, mass, ln, exn, nexts, ex_next):
    c = [c for c in centres if c <= m < nxt(c, N)][0]
    A = d(st(F(m) * F(e))); B = d(rn(F(m) * F(e))); C = et
    dropped = d(st(st(F(l) - F(t)) - F(m)))
    D = d(st(F(et) * st(1 + F(dropped))))
    dropped2 = d(rn(st(F(l) - F(t)) - F(m)))
    E = d(st(F(et) * st(1 + F(dropped2))))
    r = res.setdefault(c, [0, 0, 0, 0, 0, 0])
    r[0] += 1
    for i, v in enumerate((A, B, C, D, E)):
        r[i + 1] += (v == p)
    rows.append({"m": h(m), "mass": h(p), "ln": h(l), "exp_neg_m": h(e), "exp_next": h(et), "A": A == p, "B": B == p, "C": C == p, "D": D == p, "E": E == p})
print("centre   n    A    B    C    D    E")
for c, r in res.items():
    print(f"{c:7} {r[0]:3} " + " ".join(f"{v:4}" for v in r[1:]))
save("p3_poisson_x1_windows.json", {"rows": rows, "build": meta(ask("POISSON.DIST", [[1.0, 1.0, boolean(False)]]))})
