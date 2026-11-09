"""E1e (offline): absolute error of GAMMALN's core against math.lgamma, binned in x; does it vanish at the roots?"""
import json, math
from common import from_bits
rows = json.load(open("e1c_gammaln_core_error.json"))["rows"]
bins = {}
for r in rows:
    x, g, t = from_bits(r["x"]), from_bits(r["gammaln"]), from_bits(r["true"])
    k = round((x - 0.7) / 0.1)
    bins.setdefault(k, []).append((abs(g - t), g - t, x))
print("x-bin  n  max|abs err|  mean abs err  mean signed  (1e-17 units)")
for k in sorted(bins):
    d = bins[k]
    print(f"{0.7 + 0.1*k:4.1f} {len(d):3d} {max(a for a,_,_ in d)/1e-17:8.1f} {sum(a for a,_,_ in d)/len(d)/1e-17:8.1f} "
          f"{sum(s for _,s,_ in d)/len(d)/1e-17:8.1f}")
