"""E7: the fifteen-digit reading's half decision at large scales. Find doubles x whose exact x*10^n (n = 14 - e10)
has fractional part within a few 1e-4 of one half, ask ROUNDDOWN(x, n) and read which way the 15th digit went,
against the exact signed distance from the half in units of 2^-64 of the scaled value. Three scale classes:
one-chunk (|n| <= 27), two-chunk tiny (n ~ 314) and two-chunk huge (n ~ -286)."""
from fractions import Fraction as F
from common import *
import math, random

random.seed(71)


def design(x0, want=40, band=3e-4):
    """Scan consecutive doubles from x0 for exact scaled fractional parts near one half."""
    e10 = math.floor(math.log10(abs(x0)))
    n = 14 - e10
    scale = F(10) ** n if n >= 0 else F(1, 10 ** (-n))
    found = []
    i = 0
    while len(found) < want and i < 400000:
        x = nxt(x0, i)
        s = F(x) * scale
        frac = s - (s.numerator // s.denominator)
        d = frac - F(1, 2)
        if abs(d) < band:
            unit = s * F(2) ** -64                        # one unit of the 64-bit scaled product
            found.append((x, n, float(d / unit), float(d), s.numerator // s.denominator))
        i += 1
    return found


cases = {"one_chunk_1e-20": 1.2345678901234e-20, "two_chunk_tiny_1e-300": 1.2345678901234e-300,
         "two_chunk_huge_1e300": 1.2345678901234e300, "one_chunk_1e20": 1.2345678901234e20}
out = {}
for name, x0 in cases.items():
    des = design(x0)
    rows = [[x, float(n)] for x, n, *_ in des]
    ans = vals(ask("ROUNDDOWN", rows))
    recs = []
    for (x, n, dunits, dreal, ip), a in zip(des, ans):
        if not isinstance(a, float):
            recs.append((dunits, None)); continue
        scale = F(10) ** n if n >= 0 else F(1, 10 ** (-n))
        read = F(a) * scale
        read_int = round(read)
        up = read_int == ip + 1        # the 15-digit reading rounded the half up
        recs.append((dunits, up, h(x), h(a), n))
    recs.sort(key=lambda r: r[0])
    print(name, "n =", des[0][1], "| signed distance from half (units of 2^-64 of scaled x) -> up?")
    print("   " + " ".join(f"{d:+.2f}:{'U' if u else 'D'}" for d, u, *rest in recs))
    out[name] = [{"x": r[2], "n": r[4], "dist_units": r[0], "rounddown": r[3], "up": r[1]} for r in recs if r[1] is not None]
save("e7_trunc_half_decision.json", {"cases": out, "build": meta(ask("ROUNDDOWN", [[1.5, 0.0]]))})
