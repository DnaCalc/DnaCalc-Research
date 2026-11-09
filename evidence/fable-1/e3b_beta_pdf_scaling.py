"""E3b: for integer a and tiny x, pdf(x) = x^(a-1) K(a,b) (1+O(x)). If the pdf used POWER(x, a-1) the ratio
pdf(2^k x)/pdf(x) would be exactly 2^(k(a-1)) for dyadic x; with EXP((a-1)*LN(x)) it is not. Same test for the cdf
(which R44 already found to be an exponent calculation)."""
import random
from collections import Counter
from common import *

random.seed(33)
out = []
hist_pdf, hist_cdf = Counter(), Counter()
for a, b in [(2, 3), (3, 5), (5, 2), (8, 9), (12, 12)]:
    xs = [2.0 ** -e for e in range(60, 76)]           # exact powers of two
    pdf = vals(ask("BETA.DIST", [[x, float(a), float(b), boolean(False)] for x in xs]))
    cdf = vals(ask("BETA.DIST", [[x, float(a), float(b), boolean(True)] for x in xs]))
    for i in range(len(xs) - 1):
        # ratio of consecutive answers against the exact 2^(a-1) (pdf) and 2^a (cdf)
        rp = pdf[i] / pdf[i + 1]
        rc = cdf[i] / cdf[i + 1]
        hist_pdf[ulps(2.0 ** (a - 1), rp)] += 1
        hist_cdf[ulps(2.0 ** a, rc)] += 1
    out.append({"a": a, "b": b, "x": [h(x) for x in xs], "pdf": [h(v) for v in pdf], "cdf": [h(v) for v in cdf]})
print("pdf(x)/pdf(x/2) vs exact 2^(a-1), ulps:", dict(sorted(hist_pdf.items())))
print("cdf(x)/cdf(x/2) vs exact 2^a, ulps:", dict(sorted(hist_cdf.items())))
save("e3b_beta_pdf_scaling.json", {"rows": out, "hist_pdf": dict(hist_pdf), "hist_cdf": dict(hist_cdf),
                                   "build": meta(ask("BETA.DIST", [[0.5, 2.0, 3.0, boolean(False)]]))})
