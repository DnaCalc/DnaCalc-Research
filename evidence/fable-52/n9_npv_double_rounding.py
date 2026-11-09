"""N9: designed separation of plain binary64 from x87 double rounding in NPV's power p *= b.
Pick b in [1, 2) whose square rounds differently in one step (RN53) and in two (RN64 then RN53); rate = b - 1 is exact.
NPV(rate, 0, 0, v) = (v / (b*b)) / b then reads that product. Predictions for both models are frozen before asking."""
from common import *
rng = random.Random(5209); rows = []; pa = []; pb = []
while len(rows) < 80:
    b = 1 + rng.getrandbits(52) / 2 ** 52
    fb_ = F(b); p1, p2 = rn(fb_ * fb_), st(fb_ * fb_)
    if p1 == p2: continue
    v = rand_double(rng, -3, 12, sign=True)
    a = float(rn(rn(F(v) / p1) / fb_)); c = float(rn(rn(F(v) / p2) / fb_))
    if bits(a) == bits(c): continue
    rows.append([b - 1.0, 0.0, 0.0, v]); pa.append(a); pb.append(c)
save("n9_npv_double_rounding_frozen.json", {"rows": [[h(x) for x in r] for r in rows], "plain": [h(x) for x in pa], "x87": [h(x) for x in pb]})
reply = ask("NPV", rows); got = vals(reply)
print(meta(reply), "plain binary64:", sum(bits(g) == bits(a) for g, a in zip(got, pa)), "x87 stored:", sum(bits(g) == bits(c) for g, c in zip(got, pb)), "of", len(rows))
save("n9_npv_double_rounding.json", {"meta": meta(reply), "excel": [h(g) for g in got]})
