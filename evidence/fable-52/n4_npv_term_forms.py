"""N4: offline, on n3's single-term answers: wider family of ways to form v / b^k.
p-recurrences (p*=b, p/=b, p*=f) at 53 or 64 bits, then v/p or v*p; plus trees for k = 3."""
from common import *
d = json.load(open(os.path.join(HERE, sys.argv[1])))
fh = lambda s: F(struct.unpack(">d", bytes.fromhex(s))[0])
def forms(k, v, b):
    c = {}
    for pname, pb in (("p53", rn), ("p64", lambda x: x)):
        p = F(1)
        for _ in range(k): p = pb(r64(p / b))
        c[f"v*({pname} /=b)"] = st(v * p)
        p = F(1)
        for _ in range(k): p = pb(r64(p * b))
        c[f"v/({pname} *=b)"] = st(v / p)
        c[f"v*(1/({pname} *=b))64"] = st(v * r64(1 / p)); c[f"v*(1/({pname} *=b))53"] = st(v * rn(r64(1 / p)))
        for fname, fb_ in (("f53", rn), ("f64", lambda x: x)):
            f = fb_(r64(1 / b)); p = F(1)
            for _ in range(k): p = pb(r64(p * f))
            c[f"v*({pname} *={fname})"] = st(v * p)
            t = v
            for _ in range(k): t = pb(r64(t * f))
            c[f"chain v*{fname} ({pname})"] = rn(t)
    t = v
    for _ in range(k): t = rn(r64(t / b))
    c["chain v/b stored"] = t
    t = v
    for _ in range(k): t = r64(t / b)
    c["chain v/b resident"] = rn(t)
    return c
tally = {}
for row in d:
    k = row["k"]; e = fh(row["excel"]); c = forms(k, fh(row["v"]), rn(1 + fh(row["r"])))
    for name, x in c.items():
        tally.setdefault(name, {}).setdefault(k, 0)
        tally[name][k] += (x == e)
for name, t in sorted(tally.items(), key=lambda kv: -sum(kv[1].values())):
    print(f"{name:28s}", [t[k] for k in sorted(t)])
