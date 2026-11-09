"""E4: YIELD's small-cash limit. Two-period zero-coupon and coupon bonds, price and redemption scaled jointly by powers
of two. If the small-cash answers are 'nice' numbers (few significant bits) the solver returns a seed or a bracket
midpoint; if they are the correctly rounded roots, something else stops it."""
import math
from common import *

S, M = 43831.0, 44562.0          # 2020-01-01, 2022-01-01: exactly two annual periods, basis 0
designs = [(0.0, 25.0, 100.0), (0.0, 30.0, 100.0), (0.0, 37.0, 100.0), (0.0, 50.0, 100.0), (0.0, 80.0, 100.0),
           (0.05, 95.0, 100.0), (0.08, 103.7, 100.0)]
scales = [-500, -400, -300, -100, 0, 100, 300]


def root(rate, pr, red):
    c = 100.0 * rate
    f = lambda y: c / (1 + y) + (c + red) / (1 + y) ** 2 - pr
    lo, hi = -0.99, 50.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def sig_bits(v):
    if v == 0 or not math.isfinite(v):
        return None
    m = struct.unpack(">Q", struct.pack(">d", v))[0] & ((1 << 52) - 1)
    return 53 - (((m | (1 << 52)) & -((m | (1 << 52)))).bit_length() - 1)


import struct
rows = [[S, M, r, pr * 2.0 ** e, red * 2.0 ** e, 1.0, 0.0] for (r, pr, red) in designs for e in scales]
ans = vals(ask("YIELD", rows))
out = []
i = 0
for (r, pr, red) in designs:
    tr = root(r, pr, red)
    line = []
    for e in scales:
        v = ans[i]
        i += 1
        if isinstance(v, float):
            line.append(f"2^{e}: {v!r} ({h(v)}, {sig_bits(v)} bits, {ulps(tr, v):+d} ulp)")
        else:
            line.append(f"2^{e}: {v}")
        out.append({"rate": r, "pr": pr, "red": red, "scale_exp": e, "answer": h(v) if isinstance(v, float) else v,
                    "root_py": h(tr), "ulps_from_root": ulps(tr, v) if isinstance(v, float) else None})
    print(f"rate {r} pr {pr} red {red}: root {tr!r}")
    for l in line:
        print("   ", l)
save("e4_yield_tiny_scale.json", {"rows": out, "build": meta(ask("YIELD", [rows[0]]))})
