"""Shared helper for the fable-1 probes: typed questions to the Excel question service, interactive priority."""
import json, os, struct, sys
sys.path.insert(0, r"C:\Work\DnaCalc\DnaCalc\tools\excel-service")
from client import Client, number, boolean, number_array, from_bits  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
EX = Client()


def h(x):
    return struct.pack(">d", float(x)).hex()


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
