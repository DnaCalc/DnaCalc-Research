"""C1 (offline, on evidence/ld1 pv scans): the last period's principal P at r = 0.75, n = 1000 (q = -pv exactly) as a
two-stage product st(st(pv op1 A) op2 B) with B a 53- or 64-bit constant near a loan constant and A a free real,
by exact interval intersection over A. st = RN53(RN64(.)). Also the x87 single product st(pv * A)."""
from fractions import Fraction as F
import struct, sys, itertools
sys.path.insert(0, ".")
from common import rn, st, fb, h

def pre53(v):
    """Reals whose st (RN64 then RN53) is the double v: use the RN53 preimage widened by nothing (RN64 first only
    moves a measure-zero set); return the closed RN53 interval."""
    v = F(v)
    lo = struct.unpack(">d", struct.pack(">q", struct.unpack(">q", struct.pack(">d", float(v)))[0] - 1))[0]
    hi = struct.unpack(">d", struct.pack(">q", struct.unpack(">q", struct.pack(">d", float(v)))[0] + 1))[0]
    return (v + F(lo)) / 2, (v + F(hi)) / 2

def neighbours(x, bits, k=3):
    """The 2k+1 values of the `bits`-bit grid around the real x."""
    c = rn(F(x), bits)
    e = c.numerator.bit_length() - c.denominator.bit_length()
    ulp = F(2) ** (e - bits + 1) if e - bits + 1 >= 0 else F(1, 2 ** (bits - 1 - e))
    return [c + i * ulp for i in range(-k, k + 1)]

def intersect(ivs, new):
    out = []
    for a in ivs:
        for b in new:
            lo, hi = max(a[0], b[0]), min(a[1], b[1])
            if lo <= hi:
                out.append((lo, hi))
    return out

def fit_first_stage(samples, inner, outer):
    """P = st(outer(T)), T = st(inner(pv, A)). For each sample, the T doubles with st(outer(T)) == P, then the A interval
    from T's preimage; intersect over samples. inner: (pv, A) -> value (A free real), so A = inner_inverse(pv, T)."""
    current = None
    for pv, P in samples:
        # scan T doubles around the centre
        centre = P / outer(1.0) if outer(1.0) != 0 else P
        base = struct.unpack(">q", struct.pack(">d", float(centre)))[0]
        ivs = []
        for i in range(-600, 601):
            T = struct.unpack(">d", struct.pack(">q", base + i))[0]
            if float(st(F(outer(T)))) == P if isinstance(outer(T), F) else outer(T) == P:
                tlo, thi = pre53(T)
                a, b = inner(pv, tlo), inner(pv, thi)
                ivs.append((min(a, b), max(a, b)))
        if not ivs:
            return "no T", None
        ivs.sort()
        current = ivs if current is None else intersect(current, ivs)
        if not current:
            return "empty", None
    return f"{len(current)} interval(s)", current

rows = [l.split() for l in open("../ld1/pvscan-big1.txt") if l.strip()]
samples = [(abs(fb(a)), abs(fb(c))) for a, b, c in rows]
pmt = [abs(fb(b)) for a, b, c in rows]
print(len(samples), "samples; P/pv =", float(samples[0][1] / samples[0][0]), " 3/7 =", 3 / 7)
r = 0.75
consts = {"r": r, "1/(1+r)": 1 / 1.75, "r/(1+r)": 0.75 / 1.75, "1+r": 1.75, "1/r": 1 / 0.75, "(1+r)/r": 1.75 / 0.75}
# single product st(pv*A) for a free real A
status, ivs = fit_first_stage(samples, lambda pv, T: T / F(pv), lambda T: T)
print("single st(pv*A):", status)
# two-stage shapes
for name, cval in consts.items():
    for bits in (53, 64):
        for B in neighbours(cval, bits):
            Bf = B
            shapes = {
                "st(st(pv*A)*B)": (lambda pv, T: T / F(pv), lambda T: st(F(T) * Bf)),
                "st(st(pv*A)/B)": (lambda pv, T: T / F(pv), lambda T: st(F(T) / Bf)),
                "st(st(pv/A)*B)": (lambda pv, T: F(pv) / T, lambda T: st(F(T) * Bf)),
                "st(st(pv/A)/B)": (lambda pv, T: F(pv) / T, lambda T: st(F(T) / Bf)),
            }
            for sname, (inner, outer) in shapes.items():
                status, ivs = fit_first_stage(samples[:400], inner, outer)
                if ivs:
                    status2, ivs2 = fit_first_stage(samples, inner, outer)
                    lo, hi = ivs2[0] if ivs2 else (None, None)
                    print(f"  {sname} B={name}{'' if B == rn(F(cval), bits) else ' +-'} ({bits}b, {float(B)!r}): {status2}",
                          f"A in [{float(lo)!r}, {float(hi)!r}]" if ivs2 else "")
print("done")
