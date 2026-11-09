"""Shared helper for the fable-2 probes: typed questions to the Excel question service (interactive priority),
exact rational rounding to 64- and 53-bit significands (x87 PC64 then binary64), and bit helpers."""
import json, math, os, struct, sys
from fractions import Fraction as F
sys.path.insert(0, r"C:\Work\DnaCalc\DnaCalc\tools\excel-service")
from client import Client, number, boolean, number_array, from_bits  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
EX = Client()


def h(x):
    return struct.pack(">d", float(x)).hex()


def fb(s):
    return struct.unpack(">d", bytes.fromhex(s.removeprefix("0x")))[0]


def ulps(a, b):
    """Signed distance in ulps between two finite doubles of one sign (b - a)."""
    ia = struct.unpack(">q", struct.pack(">d", a))[0]
    ib = struct.unpack(">q", struct.pack(">d", b))[0]
    return ib - ia


def nxt(x, i):
    return struct.unpack(">d", struct.pack(">q", struct.unpack(">q", struct.pack(">d", x))[0] + i))[0]


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


# ---- exact rounding -------------------------------------------------------------------------------
def rn(x, bits=53):
    """Round a Fraction to the nearest binary float with `bits` significand bits, ties to even (normal range only)."""
    if x == 0:
        return F(0)
    s = -1 if x < 0 else 1
    x = abs(x)
    e = x.numerator.bit_length() - x.denominator.bit_length()
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


def st(x):
    """x87 store: round to 64 bits, then to 53 (the double rounding of a PC64 register stored to binary64)."""
    return rn(rn(x, 64), 53)


def d(x):
    """Fraction to float (exact when x is representable)."""
    return float(x)
