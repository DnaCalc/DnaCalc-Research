"""I8: offline on I7. Which secant iterate (x2, x3, x4, ...) is Excel's answer for two flows, per objective form."""
from common import *
from collections import Counter
from i1_irr_objective import A, S, M, D
fh = lambda s: struct.unpack(">d", bytes.fromhex(s))[0]
fz = json.load(open(os.path.join(HERE, sys.argv[1] + "_frozen.json"))); ex = json.load(open(os.path.join(HERE, sys.argv[1] + ".json")))["excel"]
def objs():
    yield "M:Q", lambda f0, f1, v: st(f0 + st(f1 * v))
    yield "Dv:Q", lambda f0, f1, v: st(f0 + st(f1 / st(1 / v)))
for name, o in objs():
    tally = Counter()
    for (fl, g), e in zip(fz["cases"], ex):
        f0, f1 = (F(fh(x)) for x in fl); e = fh(e)
        x0 = D(F(1), A(F(1), F(0.1))); y0 = o(f0, f1, x0); x1 = A(x0, F(0.001) if y0 > 0 else -F(0.001)); y1 = o(f0, f1, x1)
        hit = []
        for it in range(2, 7):
            ch = S(y1, y0)
            if ch == 0: break
            x2 = S(x1, D(M(y1, S(x1, x0)), ch)); x0, y0 = x1, y1; x1 = x2; y1 = o(f0, f1, x1)
            if bits(float(S(D(F(1), x1), F(1)))) == bits(e): hit.append(it)
        tally[tuple(hit)] += 1
    print(name, sorted(tally.items(), key=lambda kv: -kv[1]))
