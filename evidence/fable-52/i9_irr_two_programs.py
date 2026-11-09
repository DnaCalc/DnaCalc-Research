"""I9: offline on I7 (two flows). Every operation stored; the update expression and the linear term varied.
Updates (dx = x1-x0, dy = y1-y0): a x1-(y1*dx)/dy  b x1-y1*(dx/dy)  c x1-y1/(dy/dx)  d (x0*y1-x1*y0)/dy  e x0-(y0*dx)/dy
  f x1-(dx*y1)/dy (=a)  g x1 - dx/(dy/y1)  h x1 - dx*(y1/dy)
Objectives: M f0+f1*v; D f0+f1/(1/v); R (f0*b+f1)/b with b = 1/v.  Stop as landed. Also the same with stop on |y| only."""
from common import *
from i1_irr_objective import A, S, M, D
fh = lambda s: struct.unpack(">d", bytes.fromhex(s))[0]
def load(name):
    fz = json.load(open(os.path.join(HERE, name + "_frozen.json"))); ex = json.load(open(os.path.join(HERE, name + ".json")))["excel"]
    return [([F(fh(x)) for x in fl], fh(e)) for (fl, g), e in zip(fz["cases"], ex) if isinstance(e, str) and len(e) == 16]
OBJ = {"M": lambda f, v: A(f[0], M(f[1], v)), "D": lambda f, v: A(f[0], D(f[1], D(F(1), v))),
       "R": lambda f, v: (lambda b: D(A(M(f[0], b), f[1]), b))(D(F(1), v))}
UPD = {"a": lambda x0, y0, x1, y1: S(x1, D(M(y1, S(x1, x0)), S(y1, y0))), "b": lambda x0, y0, x1, y1: S(x1, M(y1, D(S(x1, x0), S(y1, y0)))),
       "c": lambda x0, y0, x1, y1: S(x1, D(y1, D(S(y1, y0), S(x1, x0)))), "d": lambda x0, y0, x1, y1: D(S(M(x0, y1), M(x1, y0)), S(y1, y0)),
       "e": lambda x0, y0, x1, y1: S(x0, D(M(y0, S(x1, x0)), S(y1, y0))), "g": lambda x0, y0, x1, y1: S(x1, D(S(x1, x0), D(S(y1, y0), y1))),
       "h": lambda x0, y0, x1, y1: S(x1, M(S(x1, x0), D(y1, S(y1, y0))))}
def run(o, u, f, stop):
    tol = F(1e-7); x0 = D(F(1), A(F(1), F(0.1))); y0 = o(f, x0)
    x1 = A(x0, F(0.001) if y0 > 0 else -F(0.001)); y1 = o(f, x1)
    for _ in range(60):
        if S(y1, y0) == 0 or y1 == 0: break
        try: x2 = u(x0, y0, x1, y1)
        except ZeroDivisionError: return None
        step = abs(S(x2, x1)); x0, y0 = x1, y1; x1 = x2; y1 = o(f, x1)
        if (stop == "both" and step < tol and abs(y1) < tol) or (stop == "y" and abs(y1) < tol) or (stop == "step" and step < tol): break
    return float(S(D(F(1), x1), F(1)))
cases = load(sys.argv[1]); res = []
for on, o in OBJ.items():
    for un, u in UPD.items():
        for stop in ("both", "y", "step"):
            res.append((sum(1 for f, e in cases if (p := run(o, u, f, stop)) is not None and bits(p) == bits(e)), on, un, stop))
for r in sorted(res, reverse=True)[:10]: print(r, "/", len(cases))
