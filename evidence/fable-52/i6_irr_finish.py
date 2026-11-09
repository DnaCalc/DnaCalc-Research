"""I6: offline on I1's answers. The final conversion from v to the rate, with the quotient kept in the x87 register:
  F1: st(st(1/v) - 1)   F3: st(r64(1/v) - 1)   and the objective shapes of I5 with all variables stored."""
from common import *
import i5_irr_masks as I5
from i1_irr_objective import A, S, M, D
def solve(o, flows, guess, fin):
    flows = [F(x) for x in flows]; g = F(guess); tol = F(1e-7)
    x0 = D(F(1), A(F(1), g)); y0 = o(flows, x0)
    if y0 == 0: return float(g)
    x1 = A(x0, F(0.001) if y0 > 0 else -F(0.001)); y1 = o(flows, x1)
    for _ in range(200):
        ch = S(y1, y0)
        if ch == 0: return None
        x2 = S(x1, D(M(y1, S(x1, x0)), ch)); step = abs(S(x2, x1)); x0, y0 = x1, y1; x1 = x2
        if x1 <= 0: return None
        y1 = o(flows, x1)
        if step < tol and abs(y1) < tol:
            return float(st(r64(1 / x1) - 1)) if fin == "F3" else float(S(D(F(1), x1), F(1)))
    return None
if __name__ == "__main__":
    fh = lambda s: struct.unpack(">d", bytes.fromhex(s))[0]
    cases = []
    for name in sys.argv[1].split(","):
        fz = json.load(open(os.path.join(HERE, name + "_frozen.json"))); ex = json.load(open(os.path.join(HERE, name + ".json")))["excel"]
        cases += [([fh(x) for x in fl], fh(g), fh(e)) for (fl, g), e in zip(fz["cases"], ex) if isinstance(e, str) and len(e) == 16]
    keys = sys.argv[2].split(",") if len(sys.argv) > 2 else ["seq:BPQT", "mul:BPQT", "div:BPQT", "horner:BPQT", "seq:bPQT"]
    for k in keys:
        for fin in ("F1", "F3"):
            ok = sum(1 for fl, g, e in cases if (p := solve(I5.CANDS[k], fl, g, fin)) is not None and bits(p) == bits(e))
            print(k, fin, ok, "/", len(cases))
