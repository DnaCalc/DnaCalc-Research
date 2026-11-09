"""E1c: accuracy of GAMMALN's core on [0.7,3): error in ulps against a 40-digit lgamma, per sub-band."""
import random
from common import *
try:
    import mpmath
    mpmath.mp.dps = 40
    true = lambda x: float(mpmath.loggamma(mpmath.mpf(x)))
    src = "mpmath 40 digits"
except ImportError:
    import math
    true = math.lgamma
    src = "math.lgamma (double)"
random.seed(13)
xs = sorted(random.uniform(0.7, 3.0) for _ in range(600))
gln = vals(ask("GAMMALN", [[x] for x in xs]))
rows = []
bands = {}
for x, g in zip(xs, gln):
    t = true(x)
    if g == 0 or t == 0:
        continue
    rel = (g - t) / abs(t)
    u = ulps(t, g)
    key = "0.7-1.0" if x < 1 else "1.0-1.5" if x < 1.5 else "1.5-2.0" if x < 2 else "2.0-2.5" if x < 2.5 else "2.5-3.0"
    d = bands.setdefault(key, [])
    d.append((x, rel, u))
    rows.append({"x": h(x), "gammaln": h(g), "true": h(t), "ulps": u, "rel": rel})
print("truth:", src)
for k in sorted(bands):
    d = bands[k]
    rels = [abs(r) for _, r, _ in d]
    us = [abs(u) for _, _, u in d]
    print(f"{k}: n={len(d)} max|rel|={max(rels):.3e} median|rel|={sorted(rels)[len(rels)//2]:.3e} "
          f"max|ulp|={max(us)} exact={sum(u==0 for u in us)} within1={sum(u<=1 for u in us)}")
save("e1c_gammaln_core_error.json", {"truth": src, "rows": rows})
