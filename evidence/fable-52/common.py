"""Helper for the fable-52 probes (Claude Fable 5.2, 2026-11-09): typed questions to the Excel question service
(interactive priority, read-only), exact rational rounding to 64 and 53 significand bits, bit helpers.
The rounding helpers follow evidence/fable-2/common.py so both lanes' scripts read the same way."""
import json, os, random, struct, sys
from fractions import Fraction as F
sys.path.insert(0, r"C:\Work\DnaCalc\DnaCalc\tools\excel-service")
from client import Client, number, boolean, from_bits  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
EX = Client()


def h(x):
    return struct.pack(">d", float(x)).hex()


def bits(x):
    return struct.unpack(">q", struct.pack(">d", float(x)))[0]


def ulps(a, b):
    """b - a in units of the last place (finite doubles of one sign)."""
    return bits(b) - bits(a)


def nxt(x, i=1):
    return struct.unpack(">d", struct.pack(">q", bits(x) + i))[0]


def ask(fn, rows, **kw):
    rows2 = [[number(v) if isinstance(v, (int, float)) else v for v in row] for row in rows]
    return EX.ask(fn, rows2, priority="interactive", **kw)


def vals(reply):
    return [from_bits(r["bits"]) if r.get("type") == "number" else r for r in reply["rows"]]


def meta(reply):
    return {k: reply.get(k) for k in ("excel_build", "excel_file_version", "env_id", "machine")}


def save(name, obj):
    with open(os.path.join(HERE, name), "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=0)


def rn(x, nbits=53):
    """Round a Fraction to nbits significand bits, nearest, ties to even (normal range only)."""
    if x == 0:
        return F(0)
    s = -1 if x < 0 else 1
    x = abs(x)
    e = x.numerator.bit_length() - x.denominator.bit_length()
    k = e - nbits + 1
    q = x / F(2) ** k if k >= 0 else x * F(2) ** (-k)
    while q >= 2 ** nbits:
        k += 1; q /= 2
    while q < 2 ** (nbits - 1):
        k -= 1; q *= 2
    n = q.numerator // q.denominator
    r = q - n
    if r > F(1, 2) or (r == F(1, 2) and n % 2 == 1):
        n += 1
    return s * F(n) * (F(2) ** k if k >= 0 else F(1, 2 ** (-k)))


def r64(x):
    return rn(x, 64)


def st(x):
    """x87 store of a PC64 register value to binary64."""
    return rn(rn(x, 64), 53)


def rand_double(rng, lo_exp, hi_exp, sign=False):
    m = 1 + rng.getrandbits(52) / 2 ** 52
    v = m * 2.0 ** rng.randint(lo_exp, hi_exp)
    return -v if sign and rng.random() < 0.5 else v
