"""N2: offline reading of n1's discovery answers by cashflow count: where do the simplest loops stop matching?"""
from common import *
from n1_npv_shapes import CANDS
d = json.load(open(os.path.join(HERE, sys.argv[1])))
rows = [[struct.unpack(">d", bytes.fromhex(x))[0] for x in r] for r in d["rows"]]
ex = [struct.unpack(">d", bytes.fromhex(x))[0] for x in d["excel"]]
by = {}
for row, e in zip(rows, ex):
    k = len(row) - 1
    r, vs = row[0], row[1:]
    exact = sum(F(v) / (1 + F(r)) ** (i + 1) for i, v in enumerate(vs))
    cr = float(rn(exact))
    a = CANDS["A1000"](r, vs)
    s = by.setdefault(k, {"n": 0, "correctly_rounded": 0, "A1000": 0, "err": []})
    s["n"] += 1; s["correctly_rounded"] += bits(cr) == bits(e); s["A1000"] += bits(a) == bits(e)
    s["err"].append(ulps(cr, e))
for k in sorted(by):
    s = by[k]; print(k, s["n"], "cr", s["correctly_rounded"], "A1000", s["A1000"], "ulp err vs exact:", sorted(set(s["err"])))
# n = 1 detail
for row, e in zip(rows, ex):
    if len(row) == 2:
        r, v = row
        b53 = rn(1 + F(r)); b64 = r64(1 + F(r))
        c = {"v/b53 one rounding": float(rn(F(v) / b53)), "v/b64 st": float(st(F(v) / b64)), "v/b53 st": float(st(F(v) / b53)),
             "v*st(1/b53)": float(st(F(v) * st(1 / b53))), "v*r64(1/b53)": float(st(F(v) * r64(1 / b53))),
             "v*r64(1/b64)": float(st(F(v) * r64(1 / b64)))}
        print(h(r), h(v), {k: ulps(e, x) for k, x in c.items()})
