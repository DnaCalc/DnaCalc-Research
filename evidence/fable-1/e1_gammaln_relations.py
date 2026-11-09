"""E1: is GAMMALN below 8 LN(GAMMA)? Where do one-step stored recurrences hold? Exact join near 0.7."""
import random
from common import *

random.seed(11)
xs = sorted([random.uniform(0.7, 8.0) for _ in range(400)] + [i / 8 for i in range(6, 64)])
gln = vals(ask("GAMMALN", [[x] for x in xs]))
gam = vals(ask("GAMMA", [[x] for x in xs]))
ln_gam = vals(ask("LN", [[g] for g in gam]))
exp_gln = vals(ask("EXP", [[v] for v in gln]))
gln_up = vals(ask("GAMMALN", [[x + 1] for x in xs]))      # GAMMALN(x+1)
ln_x = vals(ask("LN", [[x] for x in xs]))
gln_dn = vals(ask("GAMMALN", [[x - 1] for x in xs]))      # GAMMALN(x-1), only meaningful for x>1
ln_xm1 = vals(ask("LN", [[x - 1] for x in xs if x > 1]))
out = {"build": meta(ask("GAMMALN", [[2.0]])), "rows": []}
by_band = {}
for i, x in enumerate(xs):
    band = int(x) if x >= 1 else 0
    rec_up = (gln[i] == gln_up[i] - ln_x[i])
    rec_dn = None
    if x > 1:
        j = sum(1 for y in xs[:i] if y > 1)
        rec_dn = (gln[i] == gln_dn[i] + ln_xm1[j])
    d = by_band.setdefault(band, {"n": 0, "ln_gamma": 0, "exp_gln": 0, "rec_up": 0, "rec_dn": 0, "ln_gamma_within2": 0})
    d["n"] += 1
    d["ln_gamma"] += gln[i] == ln_gam[i]
    d["ln_gamma_within2"] += abs(ulps(gln[i], ln_gam[i])) <= 2 if gln[i] != 0 and ln_gam[i] != 0 else False
    d["exp_gln"] += gam[i] == exp_gln[i]
    d["rec_up"] += rec_up
    d["rec_dn"] += bool(rec_dn)
    out["rows"].append({"x": h(x), "gammaln": h(gln[i]), "gamma": h(gam[i]), "ln_gamma": h(ln_gam[i]),
                        "exp_gammaln": h(exp_gln[i]), "gammaln_x_plus_1": h(gln_up[i]), "ln_x": h(ln_x[i]),
                        "rec_up": rec_up, "rec_dn": rec_dn})
print("band: n, GAMMALN==LN(GAMMA), within2ulp, GAMMA==EXP(GAMMALN), stored recurrence up, down")
for b in sorted(by_band):
    d = by_band[b]
    print(b, d["n"], d["ln_gamma"], d["ln_gamma_within2"], d["exp_gln"], d["rec_up"], d["rec_dn"])
out["by_band"] = by_band

# Exact join of the upward recurrence near 0.7: bisect on x in bits between 0.60 and 0.80.
def rec_up_holds(x):
    a = vals(ask("GAMMALN", [[x], [x + 1]]))
    l = vals(ask("LN", [[x]]))[0]
    return a[0] == a[1] - l

grid = [0.60 + 0.01 * k for k in range(21)]
flags = [rec_up_holds(x) for x in grid]
print("recurrence-up on grid 0.60..0.80:", "".join("1" if f else "0" for f in flags))
lo, hi = None, None
for k in range(len(grid) - 1):
    if flags[k] and not flags[k + 1]:
        lo, hi = grid[k], grid[k + 1]
        break
join = None
if lo is not None:
    # bisect in the real line (floats), 40 steps
    for _ in range(45):
        mid = (lo + hi) / 2
        if rec_up_holds(mid):
            lo = mid
        else:
            hi = mid
    join = (lo, hi)
    print("join: last holds", repr(lo), h(lo), "first fails", repr(hi), h(hi))
    print("candidates: 87/128 =", 87 / 128, " 0.7 =", 0.7, " ln2 =", 0.6931471805599453)
out["join_up"] = [h(lo), h(hi)] if join else None
out["grid_flags"] = flags
save("e1_gammaln_relations.json", out)
