"""N3: isolate NPV's discount factor: NPV(r, 0, ..., 0, v) with the cashflow in position k only (k = 2..8).
Asked beside Excel's own POWER(b, k) and POWER(b, -k), b = RN53(1 + r), as oracles for a shared power leaf."""
from common import *
rng = random.Random(int(sys.argv[1]))
rows, pw, pwn, meta_rows = [], [], [], []
for k in range(2, 9):
    for _ in range(40):
        r = rand_double(rng, -7, 1) if rng.random() < 0.8 else -rand_double(rng, -6, -2)
        v = rand_double(rng, -3, 12, sign=True)
        b = float(rn(1 + F(r)))
        rows.append([r] + [0.0] * (k - 1) + [v]); pw.append([b, float(k)]); pwn.append([b, float(-k)]); meta_rows.append((k, r, v, b))
by_k = {}
for k in range(2, 9):
    idx = [i for i, m in enumerate(meta_rows) if m[0] == k]
    by_k[k] = (idx, vals(ask("NPV", [rows[i] for i in idx])))
P = vals(ask("POWER", pw)); PN = vals(ask("POWER", pwn))
out = []
for k in range(2, 9):
    idx, got = by_k[k]
    tally = {}
    for j, i in enumerate(idx):
        _, r, v, b = meta_rows[i]; e = got[j]; fb_, fv = F(b), F(v)
        d64 = F(1); d53 = F(1)
        for _ in range(k):
            d64 = r64(d64 * fb_); d53 = rn(r64(d53 * fb_))
        f64_ = r64(1 / fb_); g64 = F(1); g53 = F(1); f53 = rn(f64_)
        for _ in range(k):
            g64 = r64(g64 * f64_); g53 = rn(r64(g53 * f53))
        c = {"v/seq64": st(fv / d64), "v/seq53": st(fv / d53), "v*recipseq64": st(fv * g64), "v*recipseq53": st(fv * g53),
             "v/POWER(b,k)": st(fv / F(P[i])), "v*POWER(b,-k)": st(fv * F(PN[i])), "exact": rn(fv / fb_ ** k),
             "chain v/b/b..": None}
        t = fv
        for _ in range(k):
            t = r64(t / fb_)
        c["chain v/b/b.."] = rn(t)
        t = fv
        for _ in range(k):
            t = rn(r64(t / fb_))
        c["chain stored"] = t
        for name, x in c.items():
            tally[name] = tally.get(name, 0) + (bits(float(x)) == bits(e))
        out.append({"k": k, "r": h(r), "v": h(v), "excel": h(e), "power": h(P[i]), "power_neg": h(PN[i])})
    print(k, len(idx), tally)
save(f"n3_npv_single_term_{sys.argv[1]}.json", out)
