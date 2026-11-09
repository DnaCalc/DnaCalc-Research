"""I2: offline on I1's answers. Secant in v with the seq objective; vary the update expression and the final conversion.
  U1: x1 - (y1*(x1-x0))/(y1-y0)   U2: x1 - y1*((x1-x0)/(y1-y0))   U3: x1 - y1/((y1-y0)/(x1-x0))   U4: (x0*y1 - x1*y0)/(y1-y0)
  F1: 1/v - 1    F2: (1 - v)/v
Also reports, per case, the returned v against Excel's implied v when the rate conversion is inverted."""
from common import *
from i1_irr_objective import A, S, M, D, obj
def solve(flows, guess, upd, fin, base_of=lambda v: D(F(1), v), second=F(0.001)):
    flows = [F(x) for x in flows]; g = F(guess); tol = F(1e-7)
    f = lambda v: obj("seq", flows, base_of(v))
    x0 = D(F(1), A(F(1), g)); y0 = f(x0)
    if y0 == 0: return float(g), None
    x1 = A(x0, second if y0 > 0 else -second); y1 = f(x1)
    for it in range(200):
        ch = S(y1, y0)
        if ch == 0: return None, None
        dx = S(x1, x0)
        x2 = {"U1": lambda: S(x1, D(M(y1, dx), ch)), "U2": lambda: S(x1, M(y1, D(dx, ch))),
              "U3": lambda: S(x1, D(y1, D(ch, dx))), "U4": lambda: D(S(M(x0, y1), M(x1, y0)), ch)}[upd]()
        step = abs(S(x2, x1)); x0, y0 = x1, y1; x1 = x2
        if x1 <= 0: return None, None
        y1 = f(x1)
        if step < tol and abs(y1) < tol:
            return float(S(D(F(1), x1), F(1)) if fin == "F1" else D(S(F(1), x1), x1)), (it + 1, float(x1))
    return None, None
if __name__ == "__main__":
    fz = json.load(open(os.path.join(HERE, sys.argv[1] + "_frozen.json"))); ex = json.load(open(os.path.join(HERE, sys.argv[1] + ".json")))["excel"]
    fh = lambda s: struct.unpack(">d", bytes.fromhex(s))[0]
    cases = [([fh(x) for x in fl], fh(g)) for fl, g in fz["cases"]]
    for upd in ("U1", "U2", "U3", "U4"):
        for fin in ("F1", "F2"):
            ok = tot = 0; its = {}
            for (fl, g), e in zip(cases, ex):
                if not (isinstance(e, str) and len(e) == 16): continue
                tot += 1; p, info = solve(fl, g, upd, fin); hit = p is not None and bits(p) == bits(fh(e)); ok += hit
                if info: its.setdefault(info[0], [0, 0]); its[info[0]][0] += hit; its[info[0]][1] += 1
            print(upd, fin, ok, "/", tot, "by iterations (hit, n):", dict(sorted(its.items())))
