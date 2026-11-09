"""T4: the half decision when the nibble-chain product (K2) lands exactly on a half, one 2^-64 unit below it, or one unit
above it. Designed x for n = 314 (chunks 10, 48, 256) and n = -286 (chunks 14, 16, 256 as 10^-c factors), 16 points per class,
asked of ROUNDDOWN(x, n). The ld1 bank reading said: exact half -> Excel down (24 rows), one unit below -> Excel up (26 rows)."""
from common import *
import math, random, struct
from t1_trunc_constructions import r64, chunks, chain
random.seed(11)


def k2(x, n):
    return chain(x, n, "low", "mul")


def unit_of(v):
    e = v.numerator.bit_length() - v.denominator.bit_length()
    e = e - (1 if v < F(2) ** e else 0)
    return F(2) ** (e - 63)


def design(n, want):
    """x whose K2 value is N + 1/2 + k units (k in {-1, 0, +1}) for 15-digit integers N."""
    found = {-1: [], 0: [], 1: []}
    scale = F(10) ** n if n >= 0 else F(1, 10 ** (-n))
    tries = 0
    while any(len(v) < want for v in found.values()) and tries < 4000:
        tries += 1
        N = random.randint(10 ** 14, 10 ** 15 - 1)
        target = (F(N) + F(1, 2)) / scale
        x0 = float(target)
        for i in range(-40, 41):
            x = nxt(x0, i)
            v = k2(x, n)
            u = unit_of(v)
            k = (v - F(N) - F(1, 2)) / u
            if k.denominator == 1 and -1 <= k <= 1 and len(found[int(k)]) < want:
                found[int(k)].append((x, N))
    return found


results = {}
for n in (314, -286, 154, -126):
    des = design(n, 16)
    pts = [(k, x, N) for k in (-1, 0, 1) for x, N in des[k]]
    ans = vals(ask("ROUNDDOWN", [[x, float(n)] for k, x, N in pts]))
    score = {k: [0, 0] for k in (-1, 0, 1)}
    recs = []
    sc = F(10) ** n if n >= 0 else F(1, 10 ** (-n))
    for (k, x, N), a in zip(pts, ans):
        if not isinstance(a, float):
            continue
        up = round(F(a) * sc) == N + 1
        score[k][0] += 1; score[k][1] += up
        recs.append({"x": h(x), "n": n, "k2_offset_units": k, "rounddown": h(a), "up": up})
    print(f"n={n:5}: K2 one unit below half: up {score[-1][1]}/{score[-1][0]};  exactly half: up {score[0][1]}/{score[0][0]};"
          f"  one unit above: up {score[1][1]}/{score[1][0]}")
    results[str(n)] = recs
save("t4_trunc_ties.json", {"results": results, "build": meta(ask("ROUNDDOWN", [[1.5, 0.0]]))})
