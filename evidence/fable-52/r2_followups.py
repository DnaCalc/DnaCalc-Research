"""R2: three follow-ups to R1.
(a) MDURATION: inputs where the three division forms differ, to pick one.
(b) F.DIST against BETA.DIST where the beta argument z and 1 - z are both exact (z dyadic): interface or body?
(c) HYPGEOM.DIST single draw (n = 1, k = 1, value K/N) against a quotient and against exp-of-log-factorial sums built from
    Excel's own GAMMALN and EXP."""
from common import *
rng = random.Random(int(sys.argv[1])); T, Fa = boolean(True), boolean(False)
q = lambda fn, rows: vals(ask(fn, rows))
D = lambda a, b: float(st(F(a) / F(b))); A_ = lambda a, b: float(st(F(a) + F(b))); M_ = lambda a, b: float(st(F(a) * F(b))); S_ = lambda a, b: float(st(F(a) - F(b)))
eq = lambda a, b: isinstance(a, float) and isinstance(b, float) and bits(a) == bits(b)
# (a)
bonds = []
while len(bonds) < 400:
    s = float(rng.randint(40000, 44000)); m = s + rng.randint(200, 9000); f = float(rng.choice([1, 2, 4]))
    bonds.append([s, m, round(rng.uniform(0.005, 0.12), 4), rng.uniform(0.005, 0.15), f, float(rng.randint(0, 4))])
du, md = q("DURATION", bonds), q("MDURATION", bonds)
forms = {"d/(1+y/f)": lambda d, y, f: D(d, A_(1.0, D(y, f))), "d/((f+y)/f)": lambda d, y, f: D(d, D(A_(f, y), f)), "d*f/(f+y)": lambda d, y, f: D(M_(d, f), A_(f, y)),
         "d/(1+y/f) one rounding of the denominator": lambda d, y, f: D(d, float(rn(1 + F(y) / F(f))))}
sep = [i for i, r in enumerate(bonds) if isinstance(du[i], float) and len({bits(fn(du[i], r[3], r[4])) for fn in forms.values()}) > 1]
for name, fn in forms.items():
    print("MDURATION", name, sum(eq(md[i], fn(du[i], bonds[i][3], bonds[i][4])) for i in sep), "/", len(sep), "separating inputs")
# (b)
rows = []; zs = []
for _ in range(300):
    k = rng.randint(1, 6); t = F(2 ** k - 1) if rng.random() < 0.5 else F(1, 2 ** k - 1)
    d1 = float(rng.choice([1, 2, 3, 4, 5, 8, 16, 7, 12])); d2 = float(rng.choice([1, 2, 4, 6, 8, 10, 16, 32, 5]))
    x = t * F(d2) / F(d1)
    if F(float(x)) != x: continue
    rows.append([float(x), d1, d2]); zs.append(float(t / (1 + t)))
fd = q("F.DIST", [r + [T] for r in rows]); bd = q("BETA.DIST", [[z, r[1] / 2, r[2] / 2, T] for z, r in zip(zs, rows)])
fr = q("F.DIST.RT", rows); br = q("BETA.DIST", [[1 - z, r[2] / 2, r[1] / 2, T] for z, r in zip(zs, rows)])
print("F.DIST = BETA.DIST(z, d1/2, d2/2) with z and 1-z exact:", sum(map(eq, fd, bd)), "/", len(rows))
print("F.DIST.RT = BETA.DIST(1-z, d2/2, d1/2) with z and 1-z exact:", sum(map(eq, fr, br)), "/", len(rows))
print("F.DIST = 1 - F.DIST.RT there:", sum(eq(a, S_(1.0, b)) for a, b in zip(fd, fr) if isinstance(b, float)), " F.DIST.RT = 1 - F.DIST:", sum(eq(b, S_(1.0, a)) for a, b in zip(fd, fr) if isinstance(a, float)))
# (c)
kn = []
while len(kn) < 200:
    Np = rng.randint(3, 5000); K = rng.randint(1, Np - 1)
    kn.append((float(K), float(Np)))
hy = q("HYPGEOM.DIST", [[1.0, 1.0, K, Np, Fa] for K, Np in kn])
err = [ulps(float(rn(F(K) / F(Np))), v) for (K, Np), v in zip(kn, hy)]
print("HYPGEOM single draw: equal to RN(K/N):", sum(e == 0 for e in err), "/", len(kn), " error range in ulp:", min(err), max(err))
lnK = q("LN", [[K] for K, _ in kn]); lnN = q("LN", [[Np] for _, Np in kn])
e1 = q("EXP", [[S_(a, b)] for a, b in zip(lnK, lnN)])
print("  = EXP(LN(K) - LN(N)):", sum(map(eq, hy, e1)), "/", len(kn))
G = lambda rows_: q("GAMMALN", [[x] for x in rows_])
gK1, gK, gNK1, gNK, gN1, gN = G([K + 1 for K, _ in kn]), G([K for K, _ in kn]), G([Np - K + 1 for K, Np in kn]), G([Np - K for K, Np in kn]), G([Np + 1 for _, Np in kn]), G([Np for _, Np in kn])
# nine log-factorials for k = n = 1: lf(K) + lf(N-K) + lf(1) + lf(N-1) - lf(N) - lf(1) - lf(K-1) - lf(0) - lf(N-K)
s1 = [S_(S_(A_(A_(a, b), d), e), c) for a, b, c, d, e in zip(gK1, gNK1, gK, gN, gN1)]         # ((lfK + lf(N-K)) + lf(N-1)) - lfN - lf(K-1); lf(N-K) once
s2 = [S_(S_(a, c), S_(e, d)) for a, c, d, e in zip(gK1, gK, gN, gN1)]                          # lnC(K,1) - lnC(N,1)
for name, s in (("EXP(sum of log-factorials, numerators first)", s1), ("EXP(lnC(K,1) - lnC(N,1)) from GAMMALN", s2)):
    ev = q("EXP", [[x] for x in s]); print("  =", name + ":", sum(map(eq, hy, ev)), "/", len(kn))
save(f"r2_followups_{sys.argv[1]}.json", {"hyp_err_ulps": err})
