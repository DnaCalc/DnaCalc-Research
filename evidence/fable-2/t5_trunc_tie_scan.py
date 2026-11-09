"""T5: more tie-class points for the nibble chain K2: consecutive doubles whose exact x 10^n is within 3 units of a half are
filtered by integer arithmetic, then K2's offset from the half in 2^-64 units is read; classes -1, 0, +1 are asked of ROUNDDOWN."""
from common import *
import math, struct, sys
from t1_trunc_constructions import chunks, chain

def k2(x, n):
    return chain(x, n, "low", "mul")

def unit_of(v):
    e = v.numerator.bit_length() - v.denominator.bit_length()
    e = e - (1 if v < F(2) ** e else 0)
    return F(2) ** (e - 63)

def scan(x0, n, want, limit=3000000):
    found = {-1: [], 0: [], 1: []}
    bits0 = struct.unpack(">q", struct.pack(">d", x0))[0]
    p10 = 10 ** abs(n)
    for i in range(limit):
        x = struct.unpack(">d", struct.pack(">q", bits0 + i))[0]
        mant, e = math.frexp(x); m = int(mant * 2 ** 53); e -= 53
        num, den = (m * p10, 1) if n >= 0 else (m, p10)
        if e >= 0: num <<= e
        else: den <<= -e
        ip = num // den
        # distance of the fraction from a half, in units of V 2^-64, with V ~ 2^(bitlen-1): cheap bound first
        rem2 = 2 * (num - ip * den) - den          # 2 den (frac - 1/2)
        unit2 = (2 * num) >> 64                    # 2 den V 2^-64 approx (den V = num)
        if abs(rem2) > 3 * unit2:
            continue
        v = k2(x, n)
        k = (v - F(ip) - F(1, 2)) / unit_of(v)
        if k.denominator == 1 and -1 <= k <= 1 and len(found[int(k)]) < want:
            found[int(k)].append((x, ip))
        if all(len(t) >= want for t in found.values()):
            break
    return found

results = {}
for x0, n in [(1.2345678901234e-300, 314), (1.2345678901234e300, -286), (4.4e-140, 154), (4.4e140, -126), (2.2e-20, 34), (2.2e20, -6)]:
    des = scan(x0, n, 20)
    pts = [(k, x, ip) for k in (-1, 0, 1) for x, ip in des[k]]
    if not pts:
        print(f"n={n}: none"); continue
    ans = vals(ask("ROUNDDOWN", [[x, float(n)] for k, x, ip in pts]))
    sc = F(10) ** n if n >= 0 else F(1, 10 ** (-n))
    score = {k: [0, 0] for k in (-1, 0, 1)}; recs = []
    for (k, x, ip), a in zip(pts, ans):
        if not isinstance(a, float):
            continue
        up = round(F(a) * sc) == ip + 1
        score[k][0] += 1; score[k][1] += up
        recs.append({"x": h(x), "n": n, "k2_offset_units": k, "rounddown": h(a), "up": up})
    print(f"n={n:5}: K2 one unit below half: up {score[-1][1]}/{score[-1][0]};  exactly half: up {score[0][1]}/{score[0][0]};"
          f"  one unit above: up {score[1][1]}/{score[1][0]}", flush=True)
    results[str(n)] = recs
save("t5_trunc_tie_scan.json", {"results": results, "build": meta(ask("ROUNDDOWN", [[1.5, 0.0]]))})
