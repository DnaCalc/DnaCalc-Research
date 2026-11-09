"""T2: designed x where the candidate constructions of x * 10^n disagree on the half decision, asked of ROUNDDOWN(x, n).
Candidates (64-bit significands): K2 nibble chunks low first, RN64 each product, negative n multiplies by RN64(10^-c);
K2z the same with every product truncated (RZ); K2h half-up ties; K1 one RN64 product/quotient with the nearest 64-bit 10^|n|;
K7 chunk constants multiplied first (RN64), then one product; K4 chunks high first; K3 low first, divide for n < 0;
K8 low first, RN64, but the FIRST product (x by the low chunk) exact then the chain (x87 keeps 64 bits, so same as K2 unless
the low chunk is absent)."""
from common import *
import math, random, sys
from t1_trunc_constructions import rz, r64, chunks, K1, K7, chain

def chain_round(x, n, rounder):
    v = F(x)
    for c in chunks(n):
        f = r64(F(10) ** c) if n >= 0 else r64(F(1, 10 ** c))
        v = rounder(v * f)
    return v

def rhu(q, bits=64):
    """Round to nearest, ties up (away from zero for positives)."""
    e = q.numerator.bit_length() - q.denominator.bit_length()
    k = e - bits + 1
    scaled = q / F(2) ** k if k >= 0 else q * F(2) ** (-k)
    while scaled >= 2 ** bits:
        k += 1; scaled /= 2
    while scaled < 2 ** (bits - 1):
        k -= 1; scaled *= 2
    n = scaled.numerator // scaled.denominator
    if scaled - n >= F(1, 2):
        n += 1
    return F(n) * (F(2) ** k if k >= 0 else F(1, 2 ** (-k)))

cands = {"K2": lambda x, n: chain(x, n, "low", "mul"), "K2z": lambda x, n: chain_round(x, n, rz), "K2h": lambda x, n: chain_round(x, n, rhu),
         "K1": K1, "K7": K7, "K4": lambda x, n: chain(x, n, "high", "mul"), "K3": lambda x, n: chain(x, n, "low", "div")}

def decision(v):
    return (v - (v.numerator // v.denominator)) >= F(1, 2)

def design(x0, n, want, scan=400000, band_units=6):
    """Consecutive doubles from x0 whose exact x*10^n lies within band_units of 2^-64 V of a half, keeping those where
    the candidates disagree."""
    out = []
    bits0 = struct.unpack(">q", struct.pack(">d", x0))[0]
    for i in range(scan):
        x = struct.unpack(">d", struct.pack(">q", bits0 + i))[0]
        mant, e = math.frexp(x)
        m = int(mant * 2 ** 53); e -= 53          # x = m * 2^e
        # exact V = m * 10^n * 2^e as a fraction with power-of-two denominator when n >= 0
        if n >= 0:
            num = m * 10 ** n; den = 1
        else:
            num = m; den = 10 ** (-n)
        if e >= 0:
            num <<= e
        else:
            den <<= -e
        # fractional distance from the half in units of V * 2^-64
        ip = num // den; frac = F(num - ip * den, den)
        V = F(num, den)
        unit = V / 2 ** 64
        dist = (frac - F(1, 2)) / unit
        if abs(dist) > band_units:
            continue
        vals_ = {k: f(x, n) for k, f in cands.items()}
        decs = {k: decision(v) for k, v in vals_.items()}
        if len(set(decs.values())) > 1:
            out.append((x, float(dist), decs))
            if len(out) >= want:
                break
    return out

import struct
plan = [(1.2345678901234e-300, 314), (1.2345678901234e300, -286), (1.2345678901234e-20, 34), (1.2345678901234e20, -6),
        (7.3e-140, 154), (7.3e140, -126)]
results = {}
for x0, n in plan:
    des = design(x0, n, 24)
    if not des:
        print(f"n={n}: no disagreement found"); continue
    ans = vals(ask("ROUNDDOWN", [[x, float(n)] for x, _, _ in des]))
    score = {k: 0 for k in cands}
    recs = []
    for (x, dist, decs), a in zip(des, ans):
        if not isinstance(a, float):
            continue
        sc = F(10) ** n if n >= 0 else F(1, 10 ** (-n))
        read = F(a) * sc
        V = F(x) * sc
        up = round(read) == (V.numerator // V.denominator) + 1
        for k in cands:
            score[k] += decs[k] == up
        recs.append({"x": h(x), "n": n, "dist_units": dist, "rounddown": h(a), "up": up, "predictions": decs})
    print(f"n={n:5} points {len(recs):2}: " + "  ".join(f"{k} {v}" for k, v in score.items()))
    results[str(n)] = recs
save("t2_trunc_designed.json", {"results": results, "build": meta(ask("ROUNDDOWN", [[1.5, 0.0]]))})
