"""R1: relative oracles, bit for bit, between Excel's own functions (one batch per function, 300 random points each).
Each line is 'identity: equal / asked'. No DNA Calc code is involved: these are facts about Excel 16.0.20527 only."""
from common import *
rng = random.Random(int(sys.argv[1])); N = 300; out = {}
T, Fa = boolean(True), boolean(False)
def q(fn, rows): return vals(ask(fn, rows))
def eq(a, b): return isinstance(a, float) and isinstance(b, float) and bits(a) == bits(b)
def rep(name, pairs, extra=""):
    n = sum(eq(a, b) for a, b in pairs); d = [abs(ulps(a, b)) for a, b in pairs if isinstance(a, float) and isinstance(b, float) and a * b > 0 and not eq(a, b)]
    out[name] = [n, len(pairs)]; print(f"{name}: {n} / {len(pairs)}" + (f"  (largest difference {max(d)} ulp)" if d else "") + extra)
D = lambda a, b: float(st(F(a) / F(b))); A_ = lambda a, b: float(st(F(a) + F(b))); M_ = lambda a, b: float(st(F(a) * F(b)))
# --- bonds: MDURATION against Excel's DURATION
bonds = []
for _ in range(N):
    s = float(rng.randint(40000, 44000)); m = s + rng.randint(200, 9000); f = float(rng.choice([1, 2, 4])); b = float(rng.randint(0, 4))
    bonds.append([s, m, round(rng.uniform(0.005, 0.12), 4), round(rng.uniform(0.005, 0.15), 4), f, b])
du, md = q("DURATION", bonds), q("MDURATION", bonds)
rep("MDURATION = DURATION / (1 + yld/freq), each op stored", [(md[i], D(du[i], A_(1.0, D(r[3], r[4]))) if isinstance(du[i], float) else None) for i, r in enumerate(bonds)])
rep("MDURATION = DURATION / ((freq + yld)/freq)", [(md[i], D(du[i], D(A_(r[4], r[3]), r[4])) if isinstance(du[i], float) else None) for i, r in enumerate(bonds)])
rep("MDURATION = DURATION * freq / (freq + yld)", [(md[i], D(M_(du[i], r[4]), A_(r[4], r[3])) if isinstance(du[i], float) else None) for i, r in enumerate(bonds)])
# --- beta family
pts = [(rng.uniform(0.01, 0.99), rng.choice([0.5, 1.5, 2.0, 3.7, 8.0, 20.0, 55.5]), rng.choice([0.5, 1.5, 2.0, 4.2, 9.0, 30.0])) for _ in range(N)]
bd = q("BETA.DIST", [[x, a, b, T] for x, a, b in pts]); old = q("BETADIST", [[x, a, b] for x, a, b in pts])
rep("BETADIST(x,a,b) = BETA.DIST(x,a,b,TRUE)", list(zip(old, bd)))
sc = [(x, a, b, rng.choice([1.0, 3.0, -2.5]), rng.choice([7.0, 10.0, 12.5])) for x, a, b in pts]
o2 = q("BETADIST", [[lo + x * (hi - lo), a, b, lo, hi] for x, a, b, lo, hi in sc]); n2 = q("BETA.DIST", [[lo + x * (hi - lo), a, b, T, lo, hi] for x, a, b, lo, hi in sc])
rep("BETADIST(x,a,b,A,B) = BETA.DIST(x,a,b,TRUE,A,B)", list(zip(o2, n2)))
n3 = q("BETA.DIST", [[D(float(st(F(lo + x * (hi - lo)) - F(lo))), float(st(F(hi) - F(lo)))), a, b, T] for x, a, b, lo, hi in sc])
rep("BETA.DIST(x,a,b,TRUE,A,B) = BETA.DIST((x-A)/(B-A),a,b,TRUE)", list(zip(n2, n3)))
# --- F against beta
fp = [(rng.choice([0.05, 0.3, 0.9, 1.0, 1.7, 4.0, 12.0]) * rng.uniform(0.9, 1.1), float(rng.choice([1, 2, 3, 5, 8, 13, 40])), float(rng.choice([1, 2, 4, 7, 12, 30, 100]))) for _ in range(N)]
fd = q("F.DIST", [[x, d1, d2, T] for x, d1, d2 in fp]); frt = q("F.DIST.RT", [[x, d1, d2] for x, d1, d2 in fp]); fold = q("FDIST", [[x, d1, d2] for x, d1, d2 in fp])
rep("FDIST = F.DIST.RT", list(zip(fold, frt)))
def half(d): return d / 2
z1 = [D(M_(d1, x), A_(M_(d1, x), d2)) for x, d1, d2 in fp]; z2 = [D(d2, A_(d2, M_(d1, x))) for x, d1, d2 in fp]
b1 = q("BETA.DIST", [[z, half(d1), half(d2), T] for z, (x, d1, d2) in zip(z1, fp)]); b2 = q("BETA.DIST", [[z, half(d2), half(d1), T] for z, (x, d1, d2) in zip(z2, fp)])
rep("F.DIST(x,d1,d2,TRUE) = BETA.DIST(d1x/(d1x+d2), d1/2, d2/2)", list(zip(fd, b1)))
rep("F.DIST.RT(x,d1,d2) = BETA.DIST(d2/(d2+d1x), d2/2, d1/2)", list(zip(frt, b2)))
rep("F.DIST(x,d1,d2,TRUE) = 1 - F.DIST.RT (stored)", [(a, float(st(1 - F(b))) if isinstance(b, float) else None) for a, b in zip(fd, frt)])
rep("F.DIST.RT = 1 - F.DIST(..,TRUE) (stored)", [(b, float(st(1 - F(a))) if isinstance(a, float) else None) for a, b in zip(fd, frt)])
# --- gamma family
gp = [(rng.choice([0.2, 1.0, 3.3, 9.0, 25.0, 80.0]) * rng.uniform(0.8, 1.2), rng.choice([0.5, 1.0, 2.0, 3.0, 4.5, 10.0, 37.0])) for _ in range(N)]
g1 = q("GAMMA.DIST", [[x, a, 1.0, T] for x, a in gp]); g2 = q("GAMMA.DIST", [[2 * x, a, 2.0, T] for x, a in gp]); ch = q("CHISQ.DIST", [[2 * x, 2 * a, T] for x, a in gp])
gold = q("GAMMADIST", [[x, a, 1.0, T] for x, a in gp])
rep("GAMMADIST = GAMMA.DIST (cdf)", list(zip(gold, g1)))
rep("GAMMA.DIST(2x,a,2,TRUE) = GAMMA.DIST(x,a,1,TRUE)  (x/beta exact)", list(zip(g2, g1)))
rep("CHISQ.DIST(2x,2a,TRUE) = GAMMA.DIST(x,a,1,TRUE)", list(zip(ch, g1)))
g3 = q("GAMMA.DIST", [[3 * x, a, 3.0, T] for x, a in gp]); g3b = q("GAMMA.DIST", [[D(float(3 * x), 3.0), a, 1.0, T] for x, a in gp])
rep("GAMMA.DIST(y,a,3,TRUE) = GAMMA.DIST(stored y/3,a,1,TRUE)", list(zip(g3, g3b)))
pd1 = q("GAMMA.DIST", [[x, a, 1.0, Fa] for x, a in gp]); pd3 = q("GAMMA.DIST", [[3 * x, a, 3.0, Fa] for x, a in gp]); pd3b = q("GAMMA.DIST", [[D(float(3 * x), 3.0), a, 1.0, Fa] for x, a in gp])
rep("pdf: GAMMA.DIST(y,a,3,FALSE) = GAMMA.DIST(stored y/3,a,1,FALSE) / 3 (stored)", [(p, D(b, 3.0) if isinstance(b, float) else None) for p, b in zip(pd3, pd3b)])
ip = [(rng.uniform(0.02, 0.98), a) for _, a in gp]
i1 = q("GAMMA.INV", [[p, a, 1.0] for p, a in ip]); i3 = q("GAMMA.INV", [[p, a, 3.0] for p, a in ip]); i2 = q("GAMMA.INV", [[p, a, 2.0] for p, a in ip]); ci = q("CHISQ.INV", [[p, 2 * a] for p, a in ip])
rep("GAMMA.INV(p,a,3) = 3 * GAMMA.INV(p,a,1) (stored)", [(b, M_(3.0, a) if isinstance(a, float) else None) for a, b in zip(i1, i3)])
rep("GAMMA.INV(p,a,2) = 2 * GAMMA.INV(p,a,1) (exact scale)", [(b, 2 * a if isinstance(a, float) else None) for a, b in zip(i1, i2)])
rep("CHISQ.INV(p,2a) = GAMMA.INV(p,a,2)", list(zip(ci, i2)))
back = q("GAMMA.DIST", [[x, a, 1.0, T] for x, (p, a) in zip(i1, ip)])
rep("GAMMA.DIST(GAMMA.INV(p,a,1),a,1,TRUE) = p", [(b, p) for b, (p, a) in zip(back, ip)])
# --- hypergeometric
hp = []
for _ in range(N):
    Np = rng.randint(20, 400); K = rng.randint(1, Np - 1); n = rng.randint(1, Np - 1); lo, hi = max(0, n + K - Np), min(n, K)
    hp.append((float(rng.randint(lo, hi)), float(n), float(K), float(Np)))
h1 = q("HYPGEOM.DIST", [[k, n, K, Np, Fa] for k, n, K, Np in hp]); h2 = q("HYPGEOM.DIST", [[k, K, n, Np, Fa] for k, n, K, Np in hp])
h3 = q("HYPGEOM.DIST", [[n - k, n, Np - K, Np, Fa] for k, n, K, Np in hp]); h4 = q("HYPGEOMDIST", [[k, n, K, Np] for k, n, K, Np in hp])
rep("HYPGEOMDIST = HYPGEOM.DIST(..,FALSE)", list(zip(h4, h1)))
rep("pmf(k;n,K,N) = pmf(k;K,n,N)  (swap sample and successes)", list(zip(h1, h2)))
rep("pmf(k;n,K,N) = pmf(n-k;n,N-K,N)  (swap successes and failures)", list(zip(h1, h3)))
c1 = q("HYPGEOM.DIST", [[k, n, K, Np, T] for k, n, K, Np in hp]); c0 = q("HYPGEOM.DIST", [[k - 1, n, K, Np, T] for k, n, K, Np in hp])
rep("cdf(k) = cdf(k-1) + pmf(k) (stored)", [(c, A_(a, p) if isinstance(a, float) and isinstance(p, float) else None) for c, a, p in zip(c1, c0, h1)])
save(f"r1_relative_oracles_{sys.argv[1]}.json", out)
