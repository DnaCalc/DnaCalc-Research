"""N1: NPV(rate, v1..vn). Enumerate the loops a period developer would type, each operation in an x87 register
(64-bit significand) with every variable either stored to binary64 or kept in the register, and eliminate against
Excel on random inputs. No candidate is ranked by match count: one contradiction kills.
  A: d *= b; t += v / d            (b = 1 + rate)
  B: f = 1 / b; d *= f; t += v * d
  H: backward Horner t = (t + v) / b
Usage: python n1_npv_shapes.py discover SEED N | python n1_npv_shapes.py proof SEED N KEY[,KEY...]"""
from common import *
import itertools


def cand_A(sb, sd, sq, stt):
    def f(r, vs):
        b = r64(1 + F(r)); b = rn(b) if sb else b
        d = F(1); t = F(0)
        for v in vs:
            d = r64(d * b); d = rn(d) if sd else d
            q = r64(F(v) / d); q = rn(q) if sq else q
            t = r64(t + q); t = rn(t) if stt else t
        return float(rn(t))
    return f


def cand_B(sb, sf, sd, sq, stt):
    def f(r, vs):
        b = r64(1 + F(r)); b = rn(b) if sb else b
        fac = r64(1 / b); fac = rn(fac) if sf else fac
        d = F(1); t = F(0)
        for v in vs:
            d = r64(d * fac); d = rn(d) if sd else d
            q = r64(F(v) * d); q = rn(q) if sq else q
            t = r64(t + q); t = rn(t) if stt else t
        return float(rn(t))
    return f


def cand_H(sb, ss, stt):
    def f(r, vs):
        b = r64(1 + F(r)); b = rn(b) if sb else b
        t = F(0)
        for v in reversed(vs):
            s = r64(t + F(v)); s = rn(s) if ss else s
            t = r64(s / b); t = rn(t) if stt else t
        return float(rn(t))
    return f


CANDS = {}
for m in itertools.product((0, 1), repeat=4):
    CANDS["A" + "".join(map(str, m))] = cand_A(*m)
for m in itertools.product((0, 1), repeat=5):
    CANDS["B" + "".join(map(str, m))] = cand_B(*m)
for m in itertools.product((0, 1), repeat=3):
    CANDS["H" + "".join(map(str, m))] = cand_H(*m)


def inputs(seed, n):
    rng = random.Random(seed)
    rows = []
    for _ in range(n):
        kind = rng.random()
        if kind < 0.5:
            r = rand_double(rng, -8, 1)
        elif kind < 0.7:
            r = rand_double(rng, -60, -40)       # 1 + rate inexact in binary64, exact or nearly in 64 bits
        elif kind < 0.85:
            r = -rand_double(rng, -6, -2)
        else:
            r = rand_double(rng, 2, 12)
        k = rng.randint(1, 9)
        vs = [rand_double(rng, -3, 12, sign=True) for _ in range(k)]
        rows.append([r] + vs)
    return rows


if __name__ == "__main__":
    mode, seed, n = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    rows = inputs(seed, n)
    keys = list(CANDS) if mode == "discover" else sys.argv[4].split(",")
    pred = {k: [CANDS[k](row[0], row[1:]) for row in rows] for k in keys}   # frozen before asking
    reply = ask("NPV", rows)
    got = vals(reply)
    out = {"meta": meta(reply), "mode": mode, "seed": seed, "rows": [[h(x) for x in row] for row in rows],
           "excel": [h(g) if isinstance(g, float) else g for g in got], "score": {}}
    for k in keys:
        bad = [i for i, g in enumerate(got) if not (isinstance(g, float) and bits(g) == bits(pred[k][i]))]
        out["score"][k] = {"matched": len(rows) - len(bad), "first_kill": bad[0] if bad else None}
    alive = [k for k in keys if out["score"][k]["first_kill"] is None]
    out["alive"] = alive
    save(f"n1_npv_{mode}_{seed}.json", out)
    print(meta(reply)); print("rows", len(rows), "non-number", sum(not isinstance(g, float) for g in got))
    for k in sorted(keys, key=lambda k: -out["score"][k]["matched"])[:12]:
        print(k, out["score"][k])
    print("alive:", alive)
