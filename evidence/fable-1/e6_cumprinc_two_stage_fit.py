"""E6 (offline, on evidence/ld1 pv scans): the last period's principal at r=0.75 is pv*r/(1+r)*(1+r)^n/((1+r)^n-1).
ld1 showed it is not one rounding of pv*A for any real A. Test two-stage shapes P = g(RN(pv*A)) with g a fixed
second rounding, by exact interval intersection over A. A non-empty A interval confirms the shape."""
from fractions import Fraction as F
import struct, sys


def fb(s):
    return struct.unpack(">d", bytes.fromhex(s))[0]


def ulp(v):
    m, e = __import__("math").frexp(abs(v))
    return 2.0 ** (e - 53)


def rn(x, bits=53):
    """Round a Fraction to a binary float with `bits` significand bits, ties to even."""
    if x == 0:
        return F(0)
    s = -1 if x < 0 else 1
    x = abs(x)
    e = x.numerator.bit_length() - x.denominator.bit_length()
    # scale so that 2^(bits-1) <= x*2^-k < 2^bits
    k = e - bits + 1
    q = x / F(2) ** k if k >= 0 else x * F(2) ** (-k)
    while q >= 2 ** bits:
        k += 1; q /= 2
    while q < 2 ** (bits - 1):
        k -= 1; q *= 2
    n = q.numerator // q.denominator
    r = q - n
    if r > F(1, 2) or (r == F(1, 2) and n % 2 == 1):
        n += 1
    return s * F(n) * (F(2) ** k if k >= 0 else F(1, 2 ** (-k)))


def preimage(v):
    """Closed rational interval of reals that round (RN53) to the double v."""
    v = F(v)
    lo = (v + F(struct.unpack(">d", struct.pack(">q", struct.unpack(">q", struct.pack(">d", float(v)))[0] - 1))[0])) / 2
    hi = (v + F(struct.unpack(">d", struct.pack(">q", struct.unpack(">q", struct.pack(">d", float(v)))[0] + 1))[0])) / 2
    return lo, hi


def doubles_in(lo, hi):
    out = []
    v = float(lo)
    for i in range(-2, 3):
        c = struct.unpack(">d", struct.pack(">q", struct.unpack(">q", struct.pack(">d", v))[0] + i))[0]
        if lo <= F(c) <= hi:
            out.append(c)
    return out


def fit(samples, g, label):
    """samples: (pv, P). g: float T1 -> float. Intersect the union of A intervals per sample."""
    current = None
    for pv, P in samples:
        lo, hi = preimage(P)
        # candidate T1 doubles whose g lands exactly on P: scan a neighbourhood of T1 ~ pv*A0
        cands = []
        T = P * 1.75 if g.__name__ != "ident" else P    # centre: T1 ~ PMT = P*(1+r); scan +-400 ulps
        base = struct.unpack(">q", struct.pack(">d", T))[0]
        for i in range(-400, 401):
            t = struct.unpack(">d", struct.pack(">q", base + i))[0]
            if g(t) == P:
                cands.append(t)
        if not cands:
            return label, "no T1 candidate", None
        ivs = []
        for t in cands:
            tlo, thi = preimage(t)
            ivs.append((tlo / F(pv), thi / F(pv)))
        ivs.sort()
        if current is None:
            current = ivs
        else:
            new = []
            for a in current:
                for b in ivs:
                    lo2, hi2 = max(a[0], b[0]), min(a[1], b[1])
                    if lo2 <= hi2:
                        new.append((lo2, hi2))
            current = new
            if not current:
                return label, "empty", None
    return label, f"{len(current)} interval(s)", current


def ident(t):
    return t


def div175(t):
    return t / 1.75


def mul_recip(t):
    return t * (1 / 1.75)


def mul_recip64(t):
    return float(rn(F(t) * rn(F(1) / F(7, 4), 64)))


def minus_interest(t):
    return t - t * (0.75 / 1.75)


def minus_interest2(t):
    return t - (t * 0.75) / 1.75


for fname in ("../ld1/pvscan-big1.txt", "../ld1/pvscan-40.txt"):
    rows = [l.split() for l in open(fname) if l.strip()]
    samples = [(abs(fb(a)), abs(fb(c))) for a, b, c in rows]   # sign-symmetric: work with magnitudes
    pmt_pub = [abs(fb(b)) for a, b, c in rows]
    print(fname, len(samples), "samples")
    for g in (ident, div175, mul_recip, mul_recip64, minus_interest, minus_interest2):
        label, status, ivs = fit(samples, g, g.__name__)
        line = f"  P = {label}(RN(pv*A)): {status}"
        if ivs:
            lo, hi = ivs[0]
            line += f"  A in [{float(lo)!r}, {float(hi)!r}] width {float(hi-lo):.3e}"
            # does RN(pv*A) equal the public PMT?
            A = (lo + hi) / 2
            eq = sum(float(F(pv) * A) == pm for (pv, _), pm in zip(samples, pmt_pub))
            line += f"; RN(pv*A)==public PMT on {eq}/{len(samples)}"
        print(line)
