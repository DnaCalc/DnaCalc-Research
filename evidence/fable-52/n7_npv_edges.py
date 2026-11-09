"""N7: NPV overflow and underflow of the running power p (b = 2^100 and b = 2^-100, n cashflows of 1.0 or 2^900)."""
from common import *
from n6_npv_claim import npv
rows = []; tag = []
for n in (9, 10, 11, 12, 13):
    rows.append([2.0 ** 100] + [1.0] * n); tag.append(("b=2^100 ones", n))
    rows.append([2.0 ** 100] + [2.0 ** 900] * n); tag.append(("b=2^100 big", n))
for n in (9, 10, 11, 12, 13):
    rows.append([2.0 ** -100 - 1.0] + [1.0] * n); tag.append(("b=2^-100 ones", n))
    rows.append([2.0 ** -100 - 1.0] + [2.0 ** -900] * n); tag.append(("b=2^-100 tiny", n))
got = vals(ask("NPV", rows))
out = []
for t, r, g in zip(tag, rows, got):
    try:
        p = npv(r[0], r[1:]); ps = h(p)
    except ZeroDivisionError:
        ps = "python: division by zero"
    gs = h(g) if isinstance(g, float) else g.get("code", g)
    out.append([t, gs, ps]); print(t, "excel", gs, "claim", ps)
save("n7_npv_edges.json", out)
