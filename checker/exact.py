"""Independent finite nearest/even certificate arithmetic.

Implementer: IEEE 754 rounding and own integer/rational derivation only. No
imports from the racer, its interval library, or DnaCalc's numeric substrate.
Finite values use a signed integer coordinate on the representable lattice.
Zero signs are intentionally absent: callers must separately justify any zero
publication policy. NaNs/infinities and unsupported inversions raise errors.
"""
from dataclasses import dataclass
from fractions import Fraction as Rat
from functools import lru_cache


@lru_cache(maxsize=4096)
def two(exponent):
    return Rat(1 << exponent) if exponent >= 0 else Rat(1, 1 << -exponent)


def exponent_of(x):
    """floor(log2(x)), without a floating-point logarithm or conversion."""
    if x <= 0:
        raise ValueError('positive rational required')
    e = x.numerator.bit_length() - x.denominator.bit_length()
    return e if x >= two(e) else e - 1


@dataclass(frozen=True)
class Span:
    low: Rat
    high: Rat
    left: bool = True
    right: bool = True

    def __post_init__(self):
        if self.low > self.high or (self.low == self.high and not (self.left and self.right)):
            raise ValueError('empty span: represent it by None')

    def contains(self, x):
        return (x > self.low or (x == self.low and self.left)) and (
            x < self.high or (x == self.high and self.right))

    def meet(self, other):
        if other is None:
            return None
        lo, hi = max(self.low, other.low), min(self.high, other.high)
        lc = (self.left if lo == self.low else True) and (other.left if lo == other.low else True)
        hc = (self.right if hi == self.high else True) and (other.right if hi == other.high else True)
        return Span(lo, hi, lc, hc) if lo < hi or (lo == hi and lc and hc) else None

    def linear(self, scale, offset=Rat(0)):
        if scale == 0:
            raise ValueError('inverse needs a nonzero scale')
        if scale > 0:
            return Span(scale*self.low+offset, scale*self.high+offset, self.left, self.right)
        return Span(scale*self.high+offset, scale*self.low+offset, self.right, self.left)

    def reciprocal(self, numerator):
        if numerator == 0 or self.low*self.high <= 0:
            raise ValueError('reciprocal requires a nonzero same-sign span')
        lo, hi = numerator/self.low, numerator/self.high
        if lo < hi:
            return Span(lo, hi, self.left, self.right)
        return Span(hi, lo, self.right, self.left)

    def json(self):
        return dict(lower=[str(self.low.numerator), str(self.low.denominator)],
                    upper=[str(self.high.numerator), str(self.high.denominator)],
                    lower_closed=self.left, upper_closed=self.right)


@dataclass(frozen=True)
class Lattice:
    precision: int
    emin: int
    emax: int

    @property
    def half(self):
        return 1 << (self.precision-1)

    @property
    def last(self):
        return (self.emax-self.emin+2)*self.half-1

    def value(self, coordinate):
        """Decode an ordered coordinate; ±(last+1) are ideal overflow neighbors."""
        if abs(coordinate) > self.last+1:
            raise ValueError('coordinate out of range')
        if coordinate < 0:
            return -self.value(-coordinate)
        if coordinate < self.half:
            return coordinate*two(self.emin-self.precision+1)
        binade, fraction = divmod(coordinate, self.half)
        return (self.half+fraction)*two(self.emin+binade-self.precision)

    def floor_coordinate(self, x):
        """Unclipped at the negative end; caller checks finite bounds."""
        if x < 0:
            return -self.ceiling_coordinate(-x)
        if x >= self.value(self.last+1):
            return self.last+1
        if x == 0:
            return 0
        e = max(self.emin, exponent_of(x))
        m = x // two(e-self.precision+1)
        return (e-self.emin)*self.half+int(m)

    def ceiling_coordinate(self, x):
        if x < 0:
            return -self.floor_coordinate(-x)
        k = self.floor_coordinate(x)
        return k if self.value(k) == x else k+1

    def coordinates(self, span):
        """First/last finite lattice points in a span, respecting open endpoints."""
        if span is None:
            return None
        if span.low > self.value(self.last) or span.high < self.value(-self.last):
            return None
        a = max(-self.last, self.ceiling_coordinate(span.low))
        b = min(self.last, self.floor_coordinate(span.high))
        a += self.value(a) == span.low and not span.left
        b -= self.value(b) == span.high and not span.right
        return (a, b) if a <= b else None

    def pullback(self, span):
        """Union of the rounding cells of all lattice points in span."""
        bounds = self.coordinates(span)
        if bounds is None:
            return None
        a, b = bounds
        return Span((self.value(a-1)+self.value(a))/2,
                    (self.value(b)+self.value(b+1))/2, a % 2 == 0, b % 2 == 0)

    def round(self, x):
        """Finite nearest/even result, or None for an overflow to infinity."""
        if x < 0:
            r = self.round(-x)
            return -r if r is not None else None
        if x >= self.value(self.last+1):
            return None
        k = self.floor_coordinate(x)
        lower, upper = self.value(k), self.value(k+1)
        if x-lower < upper-x or (x-lower == upper-x and k % 2 == 0):
            return lower
        return upper if k < self.last else None


BINARY64 = Lattice(53, -1022, 1023)
EXTENDED64 = Lattice(64, -16382, 16383)
EXTENDED53 = Lattice(53, -16382, 16383)


def binary64(bits):
    """Decode bits directly as a rational; no decimal/f64 serialization."""
    b = int(bits, 16) if isinstance(bits, str) else bits
    if not 0 <= b < 1 << 64:
        raise ValueError('not a u64')
    e, f = (b >> 52) & 2047, b & ((1 << 52)-1)
    if e == 2047:
        raise ValueError('nonfinite binary64')
    m = f if e == 0 else (1 << 52)+f
    return (-1 if b >> 63 else 1)*m*two((-1022 if e == 0 else e-1023)-52)


def round_model(x, model):
    if model == 'strict':
        return BINARY64.round(x)
    fmt = EXTENDED53 if model.endswith('pc53') else EXTENDED64
    if model not in ('stored_pc64', 'resident_pc64', 'stored_pc53', 'resident_pc53'):
        raise ValueError('unknown arithmetic model')
    y = fmt.round(x)
    return BINARY64.round(y) if y is not None and model.startswith('stored') else y


def pull_rounding(span, model):
    if span is None:
        return None
    if model == 'strict':
        return BINARY64.pullback(span)
    if model not in ('stored_pc64', 'resident_pc64', 'stored_pc53', 'resident_pc53'):
        raise ValueError('unknown arithmetic model')
    fmt = EXTENDED53 if model.endswith('pc53') else EXTENDED64
    return fmt.pullback(BINARY64.pullback(span) if model.startswith('stored') else span)


def operation(x, y, op, model):
    if model == 'strict':
        x, y = BINARY64.round(x), BINARY64.round(y)
    if x is None or y is None:
        raise ValueError('nonfinite operand')
    if op == 'add':
        z = x+y
    elif op == 'mul':
        z = x*y
    elif op == 'div':
        if y == 0:
            raise ValueError('zero divisor')
        z = x/y
    else:
        raise ValueError('unknown operation')
    r = round_model(z, model)
    if r is None:
        raise ValueError('nonfinite operation')
    return r


def inverse(span, constant, op, variable_left, model):
    """Pull one arithmetic node back, including strict input stores."""
    span = pull_rounding(span, model)
    if span is None:
        return None
    if model == 'strict':
        constant = BINARY64.round(constant)
    if constant is None:
        raise ValueError('nonfinite constant')
    if op == 'add':
        span = span.linear(Rat(1), -constant)
    elif op == 'mul' and constant:
        span = span.linear(1/constant)
    elif op == 'div' and constant:
        span = span.linear(constant) if variable_left else span.reciprocal(constant)
    else:
        raise ValueError('unsupported or zero-constant inversion')
    return BINARY64.pullback(span) if model == 'strict' else span
