"""E3: at tiny x the regularized beta is pdf*x/a to the last bit. Does Excel's BETA.DIST cdf equal one rounding of its
own pdf times x over a (a shared prefactor routine), and the same for GAMMA.DIST? Plus the pdf against public leaves."""
import random
from collections import Counter
from fractions import Fraction
from common import *

random.seed(31)
pairs = [(2, 3), (4, 7), (8, 9), (1.5, 2.5), (3.3, 4.4), (10, 11), (2.5, 30), (0.5, 0.5), (1.2, 7.7), (20, 5)]
rows = []
for a, b in pairs:
    xs = [2.0 ** -e * (1 + random.random()) for e in range(60, 76)]
    cdf = vals(ask("BETA.DIST", [[x, a, b, boolean(True)] for x in xs]))
    pdf = vals(ask("BETA.DIST", [[x, a, b, boolean(False)] for x in xs]))
    gc = vals(ask("GAMMA.DIST", [[x, a, 1.0, boolean(True)] for x in xs]))
    gp = vals(ask("GAMMA.DIST", [[x, a, 1.0, boolean(False)] for x in xs]))
    for x, c, p, c2, p2 in zip(xs, cdf, pdf, gc, gp):
        m = {"m1_RN(RN(pdf*x)/a)": (p * x) / a, "m2_RN(pdf*RN(x/a))": p * (x / a), "m3_RN(RN(pdf/a)*x)": (p / a) * x,
             "m4_RN(exact pdf*x/a)": float(Fraction(p) * Fraction(x) / Fraction(a))}
        g = {"g1_RN(RN(pdf*x)/a)": (p2 * x) / a, "g2_RN(pdf*RN(x/a))": p2 * (x / a), "g3_RN(RN(pdf/a)*x)": (p2 / a) * x,
             "g4_RN(exact)": float(Fraction(p2) * Fraction(x) / Fraction(a))}
        rows.append({"a": a, "b": b, "x": h(x), "beta_cdf": h(c), "beta_pdf": h(p), "gamma_cdf": h(c2), "gamma_pdf": h(p2),
                     "beta_ulps": {k: ulps(v, c) for k, v in m.items()}, "gamma_ulps": {k: ulps(v, c2) for k, v in g.items()}})
for key in ("beta_ulps", "gamma_ulps"):
    names = rows[0][key].keys()
    for n in names:
        hist = Counter(r[key][n] for r in rows)
        print(key, n, "exact", hist.get(0, 0), "/", len(rows), "hist", dict(sorted((k, v) for k, v in hist.items() if abs(k) <= 3)))
# pdf vs public leaves (tiny x: the (b-1)ln(1-x) term is below the exponent's half ulp)
xs = [from_bits(r["x"]) for r in rows]
A = [r["a"] for r in rows]
B = [r["b"] for r in rows]
lnx = vals(ask("LN", [[x] for x in xs]))
ga = vals(ask("GAMMALN", [[a] for a in A]))
gb = vals(ask("GAMMALN", [[b] for b in B]))
gab = vals(ask("GAMMALN", [[a + b] for a, b in zip(A, B)]))
forms = {
    "f1 exp(RN((a-1)lnx) - RN(RN(ga+gb)-gab))": [((a - 1) * l) - ((x_ + y_) - z_) for a, l, x_, y_, z_ in zip(A, lnx, ga, gb, gab)],
    "f2 exp(RN(RN((a-1)lnx) + gab) - ga - gb) l-to-r": [(((a - 1) * l) + z_) - x_ - y_ for a, l, x_, y_, z_ in zip(A, lnx, ga, gb, gab)],
    "f3 exp(RN((a-1)lnx) - ga - gb + gab) l-to-r": [((a - 1) * l) - x_ - y_ + z_ for a, l, x_, y_, z_ in zip(A, lnx, ga, gb, gab)],
}
pdf = [from_bits(r["beta_pdf"]) for r in rows]
for name, expo in forms.items():
    e = vals(ask("EXP", [[v] for v in expo]))
    hist = Counter(ulps(v, p) for v, p in zip(e, pdf))
    print("pdf vs", name, "exact", hist.get(0, 0), "/", len(rows), dict(sorted((k, c) for k, c in hist.items() if abs(k) <= 4)))
save("e3_beta_pdf_cdf_sharing.json", {"rows": rows, "build": meta(ask("BETA.DIST", [[0.5, 2.0, 3.0, boolean(True)]]))})
