"""R3: F.DIST's interface to the beta. z = d1x/(d1x+d2), w = d2/(d2+d1x) (each op stored), B1 = BETA.DIST(z,d1/2,d2/2,TRUE),
B2 = BETA.DIST(w,d2/2,d1/2,TRUE). Which of B1, 1-B2 is F.DIST, which of B2, 1-B1 is F.DIST.RT, split by the side of z."""
from common import *
from collections import Counter
rng = random.Random(int(sys.argv[1])); T = boolean(True); q = lambda fn, rows: vals(ask(fn, rows))
D = lambda a, b: float(st(F(a) / F(b))); A_ = lambda a, b: float(st(F(a) + F(b))); M_ = lambda a, b: float(st(F(a) * F(b))); C1 = lambda a: float(st(1 - F(a)))
eq = lambda a, b: isinstance(a, float) and isinstance(b, float) and bits(a) == bits(b)
fp = [(rng.choice([0.05, 0.3, 0.9, 1.0, 1.7, 4.0, 12.0]) * rng.uniform(0.9, 1.1), float(rng.choice([1, 2, 3, 5, 8, 13, 40])), float(rng.choice([1, 2, 4, 7, 12, 30, 100]))) for _ in range(400)]
z = [D(M_(d1, x), A_(M_(d1, x), d2)) for x, d1, d2 in fp]; w = [D(d2, A_(d2, M_(d1, x))) for x, d1, d2 in fp]
fd = q("F.DIST", [[x, d1, d2, T] for x, d1, d2 in fp]); fr = q("F.DIST.RT", [[x, d1, d2] for x, d1, d2 in fp])
b1 = q("BETA.DIST", [[zz, d1 / 2, d2 / 2, T] for zz, (x, d1, d2) in zip(z, fp)]); b2 = q("BETA.DIST", [[ww, d2 / 2, d1 / 2, T] for ww, (x, d1, d2) in zip(w, fp)])
b3 = q("BETA.DIST", [[C1(ww), d1 / 2, d2 / 2, T] for ww, (x, d1, d2) in zip(w, fp)]); b4 = q("BETA.DIST", [[C1(zz), d2 / 2, d1 / 2, T] for zz, (x, d1, d2) in zip(z, fp)])
t = Counter()
for i in range(len(fp)):
    side = "z<.5" if z[i] < 0.5 else "z>=.5"; comp = "1-z exact" if F(z[i]) + F(C1(z[i])) == 1 and F(w[i]) == F(C1(z[i])) else "inexact"
    t[(side, "F.DIST", "B1" * eq(fd[i], b1[i]) + " 1-B2" * eq(fd[i], C1(b2[i])) + " B(1-w)" * eq(fd[i], b3[i]) + " 1-B(1-z)" * eq(fd[i], C1(b4[i])))] += 1
    t[(side, "F.DIST.RT", "B2" * eq(fr[i], b2[i]) + " 1-B1" * eq(fr[i], C1(b1[i])) + " B(1-z)" * eq(fr[i], b4[i]) + " 1-B(1-w)" * eq(fr[i], C1(b3[i])))] += 1
for k, v in sorted(t.items()): print(k, v)
print("w == stored 1 - z on", sum(bits(a) == bits(C1(b)) for a, b in zip(w, z)), "/", len(fp))
