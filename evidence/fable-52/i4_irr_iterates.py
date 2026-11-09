"""I4: offline on I1's answers. For each case, the secant iterates (seq objective, U1) continued past calc-func's stop:
which iterate, converted by 1/v - 1, is Excel's answer, and what |step| and |npv| were there. Reads the stop rule."""
from common import *
from i1_irr_objective import A, S, M, D, obj
fz = json.load(open(os.path.join(HERE, sys.argv[1] + "_frozen.json"))); ex = json.load(open(os.path.join(HERE, sys.argv[1] + ".json")))["excel"]
fh = lambda s: struct.unpack(">d", bytes.fromhex(s))[0]
rows = []
for (fl, g), e in zip(fz["cases"], ex):
    if not (isinstance(e, str) and len(e) == 16): continue
    e = fh(e); flows = [F(fh(x)) for x in fl]; g = F(fh(g))
    f = lambda v: obj("seq", flows, D(F(1), v))
    x0 = D(F(1), A(F(1), g)); y0 = f(x0); x1 = A(x0, F(0.001) if y0 > 0 else -F(0.001)); y1 = f(x1)
    hist = []; stop = None
    for it in range(1, 40):
        ch = S(y1, y0)
        if ch == 0: break
        x2 = S(x1, D(M(y1, S(x1, x0)), ch)); step = abs(S(x2, x1)); x0, y0 = x1, y1; x1 = x2; y1 = f(x1)
        rate = float(S(D(F(1), x1), F(1)))
        hist.append((it, rate, float(step), float(abs(y1))))
        if stop is None and step < F(1e-7) and abs(y1) < F(1e-7): stop = it
    best = min(hist, key=lambda t: abs(ulps(t[1], e)) if t[1] * e > 0 else 1 << 62)
    rows.append((stop, best[0], ulps(best[1], e) if best[1] * e > 0 else None, best[2], best[3], len(flows), float(max(abs(x) for x in flows))))
from collections import Counter
print("offset (excel iterate - our stop), residual ulp at that iterate:")
print(Counter((b - s if s else None, min(max(u, -3), 3) if u is not None else None) for s, b, u, *_ in rows).most_common(20))
print("examples where Excel stops later: (our stop, excel iterate, ulp, step at excel iterate, |npv| at excel iterate, n, max|flow|)")
for r in [r for r in rows if r[0] and r[1] > r[0]][:10]: print(r)
print("our-stop step/npv for cases where Excel stopped later vs same:")
def at_stop(case, e):
    fl, g = case; flows = [F(fh(x)) for x in fl]; g = F(fh(g))
    f = lambda v: obj("seq", flows, D(F(1), v))
    x0 = D(F(1), A(F(1), g)); y0 = f(x0); x1 = A(x0, F(0.001) if y0 > 0 else -F(0.001)); y1 = f(x1)
    for it in range(1, 40):
        ch = S(y1, y0); x2 = S(x1, D(M(y1, S(x1, x0)), ch)); step = abs(S(x2, x1)); xp = x1; x0, y0 = x1, y1; x1 = x2; y1 = f(x1)
        if step < F(1e-7) and abs(y1) < F(1e-7):
            r_new = S(D(F(1), x1), F(1)); r_old = S(D(F(1), xp), F(1))
            return float(step), float(abs(y1)), float(abs(r_new - r_old)), float(abs(y0))
later, same = [], []
i = 0
for case, e in zip(fz["cases"], ex):
    if not (isinstance(e, str) and len(e) == 16): continue
    r = rows[i]; i += 1
    if r[0] is None: continue
    (later if r[1] > r[0] else same).append(at_stop(case, fh(e)))
for name, grp in (("excel later", later), ("excel same", same)):
    print(name, len(grp), "step(v) range %.3g..%.3g" % (min(g[0] for g in grp), max(g[0] for g in grp)),
          "|npv| %.3g..%.3g" % (min(g[1] for g in grp), max(g[1] for g in grp)),
          "step(rate) %.3g..%.3g" % (min(g[2] for g in grp), max(g[2] for g in grp)),
          "|npv prev| %.3g..%.3g" % (min(g[3] for g in grp), max(g[3] for g in grp)))
