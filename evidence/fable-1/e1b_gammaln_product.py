"""E1b: above 3, is GAMMALN(x) = RN(GAMMALN(x-k) + LN(stored product (x-1)...(x-k))) down to [2,3)? Same for GAMMA."""
import random
from common import *

random.seed(12)
out = {"build": None, "bands": {}}
for band in range(3, 8):
    xs = sorted(random.uniform(band, band + 1) for _ in range(120))
    k = band - 2                      # steps down to [2,3)
    gln = vals(ask("GAMMALN", [[x] for x in xs]))
    core = vals(ask("GAMMALN", [[x - k] for x in xs]))
    prods_stored = []
    for x in xs:
        p = 1.0
        for j in range(1, k + 1):
            p = p * (x - j)           # binary64 stored product, low factor last? no: (x-1) first
        prods_stored.append(p)
    prods_rev = []
    for x in xs:
        p = 1.0
        for j in range(k, 0, -1):
            p = p * (x - j)
        prods_rev.append(p)
    ln_prod = vals(ask("LN", [[p] for p in prods_stored]))
    ln_prod_rev = vals(ask("LN", [[p] for p in prods_rev]))
    ln_each = [vals(ask("LN", [[x - j] for x in xs])) for j in range(1, k + 1)]
    a = sum(gln[i] == core[i] + ln_prod[i] for i in range(len(xs)))
    a2 = sum(gln[i] == core[i] + ln_prod_rev[i] for i in range(len(xs)))
    # (b) sum of logs, stored, added high to low then to core
    b = 0
    c = 0
    for i in range(len(xs)):
        s = 0.0
        for j in range(k):
            s = s + ln_each[j][i]
        b += gln[i] == core[i] + s
        t = core[i]
        for j in range(k):
            t = t + ln_each[j][i]
        c += gln[i] == t
    gam = vals(ask("GAMMA", [[x] for x in xs]))
    gcore = vals(ask("GAMMA", [[x - k] for x in xs]))
    g1 = sum(gam[i] == gcore[i] * prods_stored[i] for i in range(len(xs)))
    g2 = 0
    for i in range(len(xs)):
        t = gcore[i]
        for j in range(k, 0, -1):
            t = t * (xs[i] - j)
        g2 += gam[i] == t
    g3 = 0
    for i in range(len(xs)):
        t = gcore[i]
        for j in range(1, k + 1):
            t = t * (xs[i] - j)
        g3 += gam[i] == t
    within1 = sum(abs(ulps(gln[i], core[i] + ln_prod[i])) <= 1 for i in range(len(xs)))
    print(f"[{band},{band+1}) n={len(xs)} GAMMALN: core+LN(prod asc) {a}, core+LN(prod desc) {a2}, core+sum(ln) {b}, "
          f"seq adds {c}, within1ulp {within1} | GAMMA: core*prod {g1}, seq desc {g2}, seq asc {g3}")
    out["bands"][band] = {"x": [h(x) for x in xs], "gammaln": [h(v) for v in gln], "gammaln_core": [h(v) for v in core],
                          "gamma": [h(v) for v in gam], "gamma_core": [h(v) for v in gcore],
                          "ln_prod_asc": [h(v) for v in ln_prod], "counts": [a, a2, b, c, g1, g2, g3]}
out["build"] = meta(ask("GAMMALN", [[2.0]]))
save("e1b_gammaln_product.json", out)
