"""I5: offline on I1's answers. Objective loops with each variable stored to binary64 or kept in the x87 register.
Shapes: seq (t += flow/p; p *= base, base = 1/v), mul (t += flow*q; q *= v), div (t += flow/p; p /= v), horner.
Mask letters: B base stored, P power stored, Q term stored, T sum stored (capital = stored, lower = register).
Solver steps stay stored (U1, 1/v - 1). Optional argv[2]: comma list of keys to score only those."""
from common import *
import itertools
from i1_irr_objective import A, S, M, D
I_ = lambda x: x
def make(shape, mb, mp, mq, mt):
    kb, kp, kq, kt = (rn if m else I_ for m in (mb, mp, mq, mt))
    def o(flows, v):
        t = F(0)
        if shape == "seq":
            b = kb(r64(1 / v)); p = F(1)
            for f in flows: t = kt(r64(t + kq(r64(f / p)))); p = kp(r64(p * b))
        elif shape == "mul":
            q = F(1)
            for f in flows: t = kt(r64(t + kq(r64(f * q)))); q = kp(r64(q * v))
        elif shape == "div":
            p = F(1)
            for f in flows: t = kt(r64(t + kq(r64(f / p)))); p = kp(r64(p / v))
        else:
            for f in reversed(flows): t = kt(r64(kq(r64(t * v)) + f))
        return rn(t)
    return o
def solve(o, flows, guess, nmax=200):
    flows = [F(x) for x in flows]; g = F(guess); tol = F(1e-7)
    x0 = D(F(1), A(F(1), g)); y0 = o(flows, x0)
    if y0 == 0: return float(g)
    x1 = A(x0, F(0.001) if y0 > 0 else -F(0.001)); y1 = o(flows, x1)
    for _ in range(nmax):
        ch = S(y1, y0)
        if ch == 0: return None
        x2 = S(x1, D(M(y1, S(x1, x0)), ch)); step = abs(S(x2, x1)); x0, y0 = x1, y1; x1 = x2
        if x1 <= 0: return None
        y1 = o(flows, x1)
        if step < tol and abs(y1) < tol: return float(S(D(F(1), x1), F(1)))
    return None
def key(shape, m): return shape + ":" + "".join(c.upper() if b else c for c, b in zip("bpqt", m))
CANDS = {key(s, m): make(s, *m) for s in ("seq", "mul", "div", "horner") for m in itertools.product((0, 1), repeat=4)
         if not (s in ("mul", "div", "horner") and m[0] == 0) and not (s == "horner" and m[1] == 0)}
if __name__ == "__main__":
    fz = json.load(open(os.path.join(HERE, sys.argv[1] + "_frozen.json"))); ex = json.load(open(os.path.join(HERE, sys.argv[1] + ".json")))["excel"]
    fh = lambda s: struct.unpack(">d", bytes.fromhex(s))[0]
    cases = [([fh(x) for x in fl], fh(g), fh(e)) for (fl, g), e in zip(fz["cases"], ex) if isinstance(e, str) and len(e) == 16]
    keys = sys.argv[2].split(",") if len(sys.argv) > 2 else list(CANDS)
    res = {}
    for k in keys:
        ok = 0
        for fl, g, e in cases:
            p = solve(CANDS[k], fl, g); ok += p is not None and bits(p) == bits(e)
        res[k] = ok
    for k, v in sorted(res.items(), key=lambda kv: -kv[1]): print(k, v, "/", len(cases))
