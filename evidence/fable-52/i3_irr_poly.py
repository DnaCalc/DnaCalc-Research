"""I3: offline on I1's answers. Secant in v, U1 update, 1/v - 1 finish; the objective as a polynomial in v:
  mulseq: q = 1; t += flow*q; q *= v      horner: t = flow_n; t = t*v + flow_k      divseq: p = 1; t += flow/p; p /= v
and the seq objective with base kept in the register (not stored)."""
from common import *
import i1_irr_objective as I
from i2_irr_forms import solve
A, S, M, D = I.A, I.S, I.M, I.D
def mk(kind):
    def o(_, flows, v):
        if kind == "mulseq":
            t = F(0); q = F(1)
            for f in flows: t = A(t, M(f, q)); q = M(q, v)
        elif kind == "horner":
            t = F(0)
            for f in reversed(flows): t = A(M(t, v), f)
        elif kind == "divseq":
            t = F(0); p = F(1)
            for f in flows: t = A(t, D(f, p)); p = D(p, v)
        elif kind == "seq_resident_base":
            b = r64(1 / v); t = F(0); p = F(1)
            for f in flows: t = A(t, D(f, p)); p = M(p, b)
        return t
    return o
fz = json.load(open(os.path.join(HERE, sys.argv[1] + "_frozen.json"))); ex = json.load(open(os.path.join(HERE, sys.argv[1] + ".json")))["excel"]
fh = lambda s: struct.unpack(">d", bytes.fromhex(s))[0]
cases = [([fh(x) for x in fl], fh(g)) for fl, g in fz["cases"]]
for kind in ("mulseq", "horner", "divseq", "seq_resident_base"):
    I.obj = mk(kind); import i2_irr_forms; i2_irr_forms.obj = I.obj
    ok = tot = 0; dist = {}
    for (fl, g), e in zip(cases, ex):
        if not (isinstance(e, str) and len(e) == 16): continue
        tot += 1; p, info = solve(fl, g, "U1", "F1", base_of=lambda v: v)
        hit = p is not None and bits(p) == bits(fh(e)); ok += hit
    print(kind, ok, "/", tot)
