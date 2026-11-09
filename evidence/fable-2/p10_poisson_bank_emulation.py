"""P10: emulate the current calc-func Poisson mass (first route below 1.1 or x <= 1; Loader route above, series or stored
deviance) over every mass row of the bank (OxFunc cache, Excel 20131; x in {1,2,3}), with today's public LN and EXP as leaves,
and list where it fails."""
from common import *
from p9_poisson_series_and_fixed_d import comp_chain, deviance_series, in_series, sq2pix, r64, S
import gzip, json, collections
G = {1: 0.0, 2: fb("3fe62e42fefa39f0"), 3: fb("3ffcab0bfa2a2002")}
bank = []
with gzip.open(r"C:\Work\DnaCalc\CalcFunc-Evidence\cache\oracle\build-20131\POISSON.DIST.jsonl.gz", "rt") as f:
    for line in f:
        rr = json.loads(line); a = rr["args"][0]
        if any(v.get("kind") != "number" for v in a) or rr["outcome"]["kind"] != "number":
            continue
        x, m = a[0]["value"], a[1]["value"]; cum = a[2]["value"] if len(a) > 2 else 0
        if cum or x not in (1.0, 2.0, 3.0) or m <= 0:
            continue
        bank.append((int(x), fb(a[1]["bits_hex"][2:]), fb(rr["outcome"]["bits_hex"][2:])))
print(len(bank), "mass rows")
ms = sorted({m for x, m, p in bank})
LN = dict(zip(ms, vals(ask("LN", [[m] for m in ms]))))
first = [(i, x, m) for i, (x, m, p) in enumerate(bank) if x <= 1 or m < 1.1]
loader = [(i, x, m) for i, (x, m, p) in enumerate(bank) if not (x <= 1 or m < 1.1)]
# first route: stored exponent chain
stored_args = []
for i, x, m in first:
    t = st(F(x) * F(LN[m]))
    n1 = st(t - F(m)); n2 = st(n1 - F(G[x]))
    stored_args.append(d(n2))
EXP1 = dict(zip(sorted(set(stored_args)), vals(ask("EXP", [[a] for a in sorted(set(stored_args))]))))
# loader route
ratios = sorted({d(st(F(x) / F(m))) for i, x, m in loader if not in_series(x, m)})
LNR = dict(zip(ratios, vals(ask("LN", [[q] for q in ratios]))))
loader_args = []
for i, x, m in loader:
    if in_series(x, m):
        dev = deviance_series(x, m)
    else:
        q = d(st(F(x) / F(m)))
        dev = st(st(F(x) * F(LNR[q])) + st(F(m) - x))
    loader_args.append(d(st(-F(S[x]) - dev)))
EXP2 = dict(zip(sorted(set(loader_args)), vals(ask("EXP", [[a] for a in sorted(set(loader_args))]))))
pred = [None] * len(bank)
for (i, x, m), a in zip(first, stored_args):
    pred[i] = comp_chain(d(st(F(x) * F(LN[m]))), [m, G[x]], lambda t: EXP1[a])
for (i, x, m), a in zip(loader, loader_args):
    pred[i] = d(st(F(EXP2[a]) / F(sq2pix(x))))
fails = []
by = collections.Counter(); tot = collections.Counter()
import math
for (x, m, p), v in zip(bank, pred):
    reg = (x, "m<1.1" if m < 1.1 else "series" if in_series(x, m) else "m<x" if m < x else "m>=1.5x")
    dec = math.floor(math.log10(m))
    tot[reg] += 1
    if v != p:
        by[reg] += 1
        fails.append({"x": x, "m": h(m), "excel": h(p), "model": h(v), "ulps": ulps(p, v), "decade": dec})
print("region            rows  fails")
for reg in sorted(tot):
    print(f"{str(reg):18} {tot[reg]:5} {by[reg]:5}")
dec = collections.Counter((f["x"], f["decade"]) for f in fails)
print("fails by (x, decade of m):", sorted(dec.items()))
print("ulps histogram:", sorted(collections.Counter(f["ulps"] for f in fails).items()))
save("p10_poisson_bank_emulation.json", {"rows": len(bank), "fails": fails, "build": meta(ask("LN", [[2.0]]))})
