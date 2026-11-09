"""T1 (offline on fable-1 e7 answers): which construction of x * 10^n (n = 14 - e10, |n| > 27) predicts the measured
half decisions? Each candidate gives a 64-bit value V; the directed reading goes up iff frac(V) >= 1/2 (ties up).
Candidates: K1 nearest 64-bit constant for 10^|n|, one RN64 product (n > 0) or quotient (n < 0) [current calc-func];
K2 nibble chunks low first, each stored 10^c (RN64), product RN64 per step, n < 0 multiplies by RN64(10^-c) [S4];
K3 as K2 but n < 0 divides by RN64(10^c); K4/K5 as K2/K3 with high chunks first; K6 K1 with RZ (truncating) rounding;
K7 chunks low first with the constant built first (product of stored chunk powers, RN64 each), then one product."""
from common import *
import json, math


def rz(q, bits=64):
    """Round a positive Fraction toward zero to `bits` significand bits."""
    e = q.numerator.bit_length() - q.denominator.bit_length()
    k = e - bits + 1
    scaled = q / F(2) ** k if k >= 0 else q * F(2) ** (-k)
    while scaled >= 2 ** bits:
        k += 1; scaled /= 2
    while scaled < 2 ** (bits - 1):
        k -= 1; scaled *= 2
    n = scaled.numerator // scaled.denominator
    return F(n) * (F(2) ** k if k >= 0 else F(1, 2 ** (-k)))


def r64(q):
    return rn(q, 64)


def chunks(n):
    c = []; mag = abs(n); shift = 0
    while mag:
        part = (mag & 15) << shift
        if part:
            c.append(part)
        mag >>= 4; shift += 4
    return c


def K1(x, n):
    if n >= 0:
        return r64(F(x) * r64(F(10) ** n))
    return r64(F(x) / r64(F(10) ** (-n)))


def K6(x, n):
    if n >= 0:
        return rz(F(x) * r64(F(10) ** n))
    return rz(F(x) / r64(F(10) ** (-n)))


def chain(x, n, order, negative):
    v = F(x)
    cs = chunks(n)
    if order == "high":
        cs = cs[::-1]
    for c in cs:
        if n >= 0:
            v = r64(v * r64(F(10) ** c))
        elif negative == "mul":
            v = r64(v * r64(F(1, 10 ** c)))
        else:
            v = r64(v / r64(F(10) ** c))
    return v


def K7(x, n):
    k = F(1)
    for c in chunks(n):
        k = r64(k * r64(F(10) ** c))
    return r64(F(x) * k) if n >= 0 else r64(F(x) / k)


cands = {"K1 nearest, one op": K1, "K2 low, mul 10^-c": lambda x, n: chain(x, n, "low", "mul"),
         "K3 low, div 10^c": lambda x, n: chain(x, n, "low", "div"), "K4 high, mul 10^-c": lambda x, n: chain(x, n, "high", "mul"),
         "K5 high, div 10^c": lambda x, n: chain(x, n, "high", "div"), "K6 nearest, RZ": K6, "K7 constant first": K7}
data = json.load(open("../fable-1/e7_trunc_half_decision.json"))["cases"]
for name, rows in data.items():
    print(name, "n =", rows[0]["n"], "chunks", chunks(int(rows[0]["n"])))
    for cname, f in cands.items():
        agree = 0; ties = 0
        for r in rows:
            x = fb(r["x"]); n = int(r["n"])
            v = f(x, n)
            frac = v - (v.numerator // v.denominator)
            up = frac >= F(1, 2)
            ties += frac == F(1, 2)
            agree += (up == r["up"])
        print(f"   {cname:22} agree {agree:2}/{len(rows)}  exact-half products {ties}")
