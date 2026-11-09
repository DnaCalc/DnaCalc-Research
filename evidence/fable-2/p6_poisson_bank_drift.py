"""P6: do the bank's POISSON.DIST mass answers (Excel build 20131, OxFunc cache) equal Excel 20527's answers today?
Re-ask 100 rows per count x in {1, 2, 3}, spread over the mean's range; also GAMMALN(2..4) and LN(2), LN(6)."""
from common import *
import gzip, json, random
random.seed(5)
bank = {1: [], 2: [], 3: []}
with gzip.open(r"C:\Work\DnaCalc\CalcFunc-Evidence\cache\oracle\build-20131\POISSON.DIST.jsonl.gz", "rt") as f:
    for line in f:
        r = json.loads(line); a = r["args"][0]
        if any(v.get("kind") != "number" for v in a) or r["outcome"]["kind"] != "number":
            continue
        x, m = a[0]["value"], a[1]["value"]; cum = a[2]["value"] if len(a) > 2 else 0
        if cum or x not in bank:
            continue
        bank[int(x)].append((a[1]["bits_hex"][2:], r["outcome"]["bits_hex"][2:]))
rows = []
for x in bank:
    bank[x].sort(key=lambda t: fb(t[0]))
    step = max(1, len(bank[x]) // 100)
    rows += [(x, mh, ob) for mh, ob in bank[x][::step]]
ans = ask("POISSON.DIST", [[float(x), number(mh), boolean(False)] for x, mh, ob in rows])["rows"]
by = {}
out = []
for (x, mh, ob), a in zip(rows, ans):
    m = fb(mh)
    reg = "m<1.1" if m < 1.1 else "m<10" if m < 10 else "m<100" if m < 100 else "m>=100"
    same = a.get("bits") == ob
    dd = by.setdefault((x, reg), [0, 0]); dd[0] += 1; dd[1] += same
    out.append({"x": x, "m": mh, "bank_20131": ob, "now": a.get("bits"), "same": same,
                "ulps": ulps(fb(ob), fb(a["bits"])) if "bits" in a else None})
print("x region   n  same")
for (x, reg), (n, s) in sorted(by.items()):
    print(f"{x} {reg:7} {n:3} {s:3}")
g = ask("GAMMALN", [[2.0], [3.0], [4.0]])["rows"]
l = ask("LN", [[2.0], [6.0]])["rows"]
print("GAMMALN(2,3,4):", [r.get("bits") for r in g], " LN(2), LN(6):", [r.get("bits") for r in l])
save("p6_poisson_bank_drift.json", {"rows": out, "gammaln_2_3_4": [r.get("bits") for r in g], "ln_2_6": [r.get("bits") for r in l],
                                    "build": meta(ask("POISSON.DIST", [[1.0, 1.0, boolean(False)]]))})
