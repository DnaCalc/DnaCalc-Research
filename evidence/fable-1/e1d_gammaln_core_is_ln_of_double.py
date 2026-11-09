"""E1d: in the core [0.7,3), is GAMMALN(x) = RN(ln(G)) for a binary64 G? Near the roots the output's ulp is far
finer than G's, so exp(GAMMALN) must land on the binary64 grid if G is a stored double. Also: is G the public GAMMA?"""
import json
from decimal import Decimal, getcontext
from fractions import Fraction
import struct
from common import *

getcontext().prec = 45
with open("e1c_gammaln_core_error.json") as stream:
    rows = json.load(stream)["rows"]
xs = [from_bits(r["x"]) for r in rows]
gln = [from_bits(r["gammaln"]) for r in rows]
gam = vals(ask("GAMMA", [[x] for x in xs]))


def ulp(v):
    return abs(nxt(v, 1) - v)


out = []
stats = {"n": 0, "on_grid": 0, "chance_sum": 0.0, "G_equals_public_GAMMA": 0, "G_within1_public": 0, "tested": 0}
for x, g, pub in zip(xs, gln, gam):
    if g == 0:
        continue
    G = Decimal(g).exp()
    Gf = float(G)
    cands = [nxt(Gf, i) for i in (-2, -1, 0, 1, 2)]
    # distance from G to the nearest double, in units of ulp(G)
    dist = min(abs(G - Decimal(c)) for c in cands)
    nearest = min(cands, key=lambda c: abs(G - Decimal(c)))
    du = float(dist / Decimal(ulp(Gf)))
    # preimage half-width of RN(ln(.)) = g, in G units, relative to ulp(G)
    half = Decimal(ulp(g)) / 2 * G
    half_u = float(half / Decimal(ulp(Gf)))
    decisive = half_u < 0.25         # output pins G to a quarter ulp or better
    on_grid = dist <= half + Decimal("1e-19")
    r = {"x": h(x), "gammaln": h(g), "G_nearest": h(nearest), "dist_ulpG": du, "preimage_halfwidth_ulpG": half_u,
         "decisive": decisive, "on_grid": on_grid, "public_GAMMA": h(pub), "nearest_eq_public": nearest == pub,
         "nearest_vs_public_ulps": ulps(pub, nearest)}
    out.append(r)
    if decisive:
        stats["tested"] += 1
        stats["on_grid"] += on_grid
        stats["chance_sum"] += 2 * half_u          # probability a random real lands within +-half of a grid point
        stats["G_equals_public_GAMMA"] += nearest == pub
        stats["G_within1_public"] += abs(ulps(pub, nearest)) <= 1
    stats["n"] += 1
print(stats)
print("expected on-grid by chance:", round(stats["chance_sum"], 1), "of", stats["tested"])
# per band
for lo, hi in [(0.7, 1), (1, 1.5), (1.5, 2), (2, 2.5), (2.5, 3)]:
    sel = [r for r in out if lo <= from_bits(r["x"]) < hi and r["decisive"]]
    if sel:
        print(f"[{lo},{hi}): decisive {len(sel)}, on grid {sum(r['on_grid'] for r in sel)}, "
              f"chance {sum(2*r['preimage_halfwidth_ulpG'] for r in sel):.1f}, G==GAMMA {sum(r['nearest_eq_public'] for r in sel)}")
save("e1d_gammaln_core_is_ln_of_double.json", {"stats": stats, "rows": out, "build": meta(ask("GAMMA", [[2.0]]))})
