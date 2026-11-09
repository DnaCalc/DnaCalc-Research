"""N5: offline on n3: split forms. (a) (v/b)/b^(k-1); (b) (v/b^(k-1))/b; power sequential or by squaring, 53 or 64 bits."""
from common import *
d = json.load(open(os.path.join(HERE, sys.argv[1])))
fh = lambda s: F(struct.unpack(">d", bytes.fromhex(s))[0])
def seq(b, n, keep):
    p = F(1)
    for _ in range(n): p = keep(r64(p * b))
    return p
def sq(b, n, keep):          # right-to-left binary powering
    res = F(1); base = b
    while n:
        if n & 1: res = keep(r64(res * base))
        n >>= 1
        if n: base = keep(r64(base * base))
    return res
tally = {}
I = lambda x: x
for row in d:
    k = row["k"]; e = fh(row["excel"]); v = fh(row["v"]); b = rn(1 + fh(row["r"]))
    c = {}
    for pn, pf in (("seq", seq), ("sq", sq)):
        for kn, keep in (("53", rn), ("64", I)):
            p = pf(b, k - 1, keep)
            for qn, qk in (("q53", rn), ("q64", I)):
                c[f"(v/b)/{pn}{kn}^(k-1) {qn}"] = st(qk(r64(v / b)) / p)
                c[f"(v/{pn}{kn}^(k-1))/b {qn}"] = st(qk(r64(v / p)) / b)
            p = pf(b, k, keep)
            c[f"v/{pn}{kn}^k"] = st(v / p)
    for name, x in c.items():
        tally.setdefault(name, {}).setdefault(k, 0); tally[name][k] += (x == e)
for name, t in sorted(tally.items(), key=lambda kv: -sum(kv[1].values()))[:10]:
    print(f"{name:30s}", [t[k] for k in sorted(t)])
