"""N10: the same designed separation for NPV's division and addition.
(a) NPV(b - 1, v) = v / b with RN53(v/b) != RN53(RN64(v/b));  (b) NPV(0, v1, v2) = v1 + v2 with RN53 != RN53(RN64)."""
from common import *
rng = random.Random(5210); rows = []; plain = []; x87 = []; kind = []
while kind.count("div") < 60:
    b = 1 + rng.getrandbits(52) / 2 ** 52; v = rand_double(rng, -3, 12, sign=True)
    a, c = rn(F(v) / F(b)), st(F(v) / F(b))
    if a != c:
        rows.append([b - 1.0, v]); plain.append(float(a)); x87.append(float(c)); kind.append("div")
while kind.count("add") < 60:
    v1 = rand_double(rng, 0, 4, sign=True); v2 = rand_double(rng, -20, -12, sign=True)
    a, c = rn(F(v1) + F(v2)), st(F(v1) + F(v2))
    if a != c:
        rows.append([0.0, v1, v2]); plain.append(float(a)); x87.append(float(c)); kind.append("add")
save("n10_npv_div_add_frozen.json", {"rows": [[h(x) for x in r] for r in rows], "plain": [h(x) for x in plain], "x87": [h(x) for x in x87]})
got = []
for k in ("div", "add"):
    idx = [i for i, t in enumerate(kind) if t == k]
    g = vals(ask("NPV", [rows[i] for i in idx])); got += g
    print(k, "plain binary64:", sum(bits(g[j]) == bits(plain[i]) for j, i in enumerate(idx)),
          "x87 stored:", sum(bits(g[j]) == bits(x87[i]) for j, i in enumerate(idx)), "of", len(idx))
save("n10_npv_div_add.json", {"excel": [h(g) for g in got]})
