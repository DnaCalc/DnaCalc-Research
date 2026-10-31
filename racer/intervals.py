"""Exact finite IEEE rounding preimages; research, never Excel observations.

Implementer; tier 1 IEEE 754 nearest/even and rational algebra. No restricted
source read. Endpoints retain their open/closed flags through every inverse.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
import math
import random
import struct
from functools import lru_cache


@lru_cache(maxsize=None)
def pow2(k):
    return Q(1 << k) if k >= 0 else Q(1, 1 << -k)


def log2floor(x):
    assert x > 0
    k = x.numerator.bit_length() - x.denominator.bit_length()
    return k - (x < pow2(k))


@dataclass(frozen=True)
class Format:
    p: int
    emin: int
    emax: int

    @property
    @lru_cache(maxsize=None)
    def maximum(self):
        return (2 - pow2(1 - self.p)) * pow2(self.emax)

    def floor(self, x):
        if x < 0:
            y = self.ceil(-x)
            return -y if y is not None else None
        if x > self.maximum:
            return self.maximum
        if x == 0:
            return Q(0)
        step = pow2(max(self.emin, log2floor(x)) - self.p + 1)
        return (x // step) * step

    def ceil(self, x):
        if x < 0:
            y = self.floor(-x)
            # floor above +MAX is MAX, but negative ceil stays -MAX.
            return -y
        if x > self.maximum:
            return None
        y = self.floor(x)
        if y == x:
            return y
        return self.next(y)

    def next(self, x):
        assert -self.maximum <= x <= self.maximum
        if x == self.maximum:
            return None
        if x == 0:
            return pow2(self.emin - self.p + 1)
        if x < 0:
            return -self.prev(-x)
        return x + pow2(max(self.emin, log2floor(x)) - self.p + 1)

    def prev(self, x):
        if x == -self.maximum:
            return None
        if x <= 0:
            return -self.next(-x)
        e = log2floor(x)
        step = pow2(max(self.emin, e - (x == pow2(e))) - self.p + 1)
        return x - step

    def even(self, x):
        if x == 0:
            return True
        step = pow2(max(self.emin, log2floor(abs(x))) - self.p + 1)
        return int(x / step) % 2 == 0

    def rn(self, x):
        if x < 0:
            y = self.rn(-x)
            return -y if y is not None else None
        lo, hi = self.floor(x), self.ceil(x)
        if hi is None:
            # The ideal successor to MAX for overflow is 2^(emax+1).
            return lo if x < (lo + pow2(self.emax + 1)) / 2 else None
        dl, dh = x - lo, hi - x
        return lo if dl < dh or (dl == dh and self.even(lo)) else hi

    @lru_cache(maxsize=131072)
    def preimage(self, i):
        if i is None:
            return None
        low, high = self.ceil(i.lo), self.floor(i.hi)
        if low is None or high is None:
            return None
        if low == i.lo and not i.lc:
            low = self.next(low)
        if high == i.hi and not i.hc:
            high = self.prev(high)
        if low is None or high is None or low > high:
            return None
        prv, nxt = self.prev(low), self.next(high)
        if prv is None:
            prv = -pow2(self.emax + 1)
        if nxt is None:
            nxt = pow2(self.emax + 1)
        return Interval((prv + low) / 2, (high + nxt) / 2,
                        self.even(low), self.even(high))


B64 = Format(53, -1022, 1023)
E64 = Format(64, -16382, 16383)


@dataclass(frozen=True)
class Interval:
    lo: Q
    hi: Q
    lc: bool = True
    hc: bool = True

    def intersect(self, other):
        if other is None:
            return None
        low, high = max(self.lo, other.lo), min(self.hi, other.hi)
        lc = (self.lc if low == self.lo else True) and (other.lc if low == other.lo else True)
        hc = (self.hc if high == self.hi else True) and (other.hc if high == other.hi else True)
        return Interval(low, high, lc, hc) if low < high or (low == high and lc and hc) else None

    def affine(self, a, b=Q(0)):
        assert a != 0
        if a > 0:
            return Interval(a * self.lo + b, a * self.hi + b, self.lc, self.hc)
        return Interval(a * self.hi + b, a * self.lo + b, self.hc, self.lc)

    def reciprocal(self, c):
        assert self.lo * self.hi > 0 and c != 0, 'nonzero sign interval required'
        a, b = c / self.lo, c / self.hi
        return Interval(a, b, self.lc, self.hc) if a < b else Interval(b, a, self.hc, self.lc)

    def contains(self, x):
        return (x > self.lo or (x == self.lo and self.lc)) and (x < self.hi or (x == self.hi and self.hc))

    def json(self):
        def pair(x):
            return [str(x.numerator), str(x.denominator)]
        return dict(lower=pair(self.lo), upper=pair(self.hi), lower_closed=self.lc, upper_closed=self.hc)


def number(bits):
    return struct.unpack('>d', int(bits.removeprefix('0x'), 16).to_bytes(8, 'big'))[0]


def stored(x):
    y = E64.rn(x)
    return B64.rn(y) if y is not None else None


def rounding_inverse(i, model):
    if model == 0:
        return B64.preimage(i)
    if model == 1:
        return E64.preimage(B64.preimage(i))
    assert model == 2
    return E64.preimage(i)


def inverse_op(i, model, op, c, variable_left):
    """Invert a single arithmetic node, including strict input stores."""
    i = rounding_inverse(i, model)
    if i is None:
        return None
    if model == 0:
        c = B64.rn(c)
    assert c is not None and c != 0
    if op == 'mul':
        i = i.affine(1 / c)
    elif variable_left:
        i = i.affine(c)
    else:
        i = i.reciprocal(c)
    return B64.preimage(i) if model == 0 else i


def self_test():
    rng = random.Random(0x524e375f494e54)
    n = 0
    for bits in [0, 1, 2, 0xfffffffffffff, 0x10000000000000, 0x3ff0000000000000,
                 0x7fefffffffffffff] + [rng.randrange(0x7ff0000000000000) for _ in range(2000)]:
        f = number(f'{bits:016x}')
        x = Q(f)
        assert B64.rn(x) == x
        if bits < 0x7fefffffffffffff:
            nxt = Q(math.nextafter(f, math.inf))
            assert B64.next(x) == nxt
            mid = (x + nxt) / 2
            assert B64.rn(mid) == (x if bits % 2 == 0 else nxt)
        if bits:
            prv = Q(math.nextafter(f, -math.inf))
            assert B64.prev(x) == prv
        interval = B64.preimage(Interval(x, x))
        assert interval.contains(x)
        assert interval.lc == (bits % 2 == 0)
        assert interval.hc == (bits % 2 == 0)
        assert B64.preimage(Interval(-x, -x)).contains(-x)
        n += 1
    # Exact absence between adjacent representables, with each endpoint open.
    assert B64.preimage(Interval(Q(1), Q(1) + pow2(-52), False, False)) is None
    print({'rounding_preimage_controls': n, 'arithmetic': 'exact Fraction, ties and endpoint flags'})


if __name__ == '__main__':
    self_test()
