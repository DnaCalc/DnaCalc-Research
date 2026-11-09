"""E2: (a) is NORM.S.INV(1-p) == -NORM.S.INV(p) bit for bit, by region? (b) does NORM.INV(p,0,s) leak an extended
intermediate, i.e. differ from RN(s*NORM.S.INV(p))?"""
import random
from collections import Counter
from common import *

random.seed(21)
out = {"mirror": {}, "leak": {}}
regions = [("upper_tail_p>15/16", 0.9375, 0.995), ("piece_C_(0.925,15/16]", 0.925, 0.9375),
           ("upper_B_(0.85<t<=0.925)", 0.85, 0.925), ("upper_A_(0.5,0.85)", 0.5, 0.85)]
for name, lo, hi in regions:
    ph = [random.uniform(lo, hi) for _ in range(300)]
    pl = [1 - p for p in ph]                      # exact for p >= 0.5
    a = vals(ask("NORM.S.INV", [[p] for p in pl]))
    b = vals(ask("NORM.S.INV", [[p] for p in ph]))
    diffs = Counter(ulps(-y, x) for x, y in zip(a, b))
    print(name, "exact mirror", diffs.get(0, 0), "/300; ulp(NORM.S.INV(p) + NORM.S.INV(1-p)) histogram:", dict(sorted(diffs.items())))
    out["mirror"][name] = {"p_hi": [h(p) for p in ph], "inv_lo": [h(v) for v in a], "inv_hi": [h(v) for v in b],
                           "hist": dict(diffs)}
# dense windows around 1/16 and 15/16: is there a split at 1/16 on the lower side too?
for centre in (0.0625, 0.9375):
    ps = [nxt(centre, i) for i in range(-48, 49)]
    v = vals(ask("NORM.S.INV", [[p] for p in ps]))
    m = vals(ask("NORM.S.INV", [[1 - p] for p in ps]))
    eq = [x == -y for x, y in zip(v, m)]
    print(f"window at {centre}: mirror-exact pattern:", "".join("1" if e else "0" for e in eq))
    out["mirror"][f"window_{centre}"] = {"p": [h(p) for p in ps], "inv": [h(x) for x in v], "inv_mirror": [h(x) for x in m]}

# (b) NORM.INV leak
ps = [random.uniform(0.08, 0.92) for _ in range(400)]
ss = [random.choice([3.0, 5.0, 7.0, 1.0, random.uniform(0.3, 9)]) for _ in ps]
z = vals(ask("NORM.S.INV", [[p] for p in ps]))
ni = vals(ask("NORM.INV", [[p, 0.0, s] for p, s in zip(ps, ss)]))
pred = [s * zz for s, zz in zip(ss, z)]
d = Counter(ulps(pr, n) for pr, n in zip(pred, ni))
print("NORM.INV(p,0,s) vs RN(s*NORM.S.INV(p)):", dict(sorted(d.items())))
ni1 = vals(ask("NORM.INV", [[p, 0.0, 1.0] for p in ps]))
print("NORM.INV(p,0,1) == NORM.S.INV(p):", sum(a == b for a, b in zip(ni1, z)), "/", len(ps))
ms = [random.uniform(-3, 3) for _ in ps]
nim = vals(ask("NORM.INV", [[p, m, s] for p, m, s in zip(ps, ms, ss)]))
pred2 = [m + s * zz for m, s, zz in zip(ms, ss, z)]
d2 = Counter(ulps(pr, n) for pr, n in zip(pred2, nim))
print("NORM.INV(p,m,s) vs RN(m + RN(s*z)):", dict(sorted(d2.items())))
out["leak"] = {"p": [h(p) for p in ps], "s": [h(s) for s in ss], "m": [h(m) for m in ms], "norm_s_inv": [h(v) for v in z],
               "norm_inv_0_s": [h(v) for v in ni], "norm_inv_0_1": [h(v) for v in ni1], "norm_inv_m_s": [h(v) for v in nim],
               "hist_0_s": dict(d), "hist_m_s": dict(d2)}
out["build"] = meta(ask("NORM.S.INV", [[0.5]]))
save("e2_norm_inverse_symmetry_leak.json", out)
