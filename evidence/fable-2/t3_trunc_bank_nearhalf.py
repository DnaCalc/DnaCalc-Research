"""T3 (offline on the ROUNDDOWN bank, Excel 20430, OxFunc cache): rows with |n| > 27 whose exact x 10^n lies within 3 units
of 2^-64 V of a half; score the constructions' half decisions against Excel's answers. Adds K9 (low nibble as repeated
x87 products by 10, higher chunks from the table) and K10 (every step a product by 10 or by RN64(0.1))."""
from common import *
import gzip, json, math, struct
from t1_trunc_constructions import rz, r64, chunks, K1, K7, chain
from t2_trunc_designed import chain_round, rhu, decision

def K9(x, n):
    v = F(x); cs = chunks(n)
    for i, c in enumerate(cs):
        if i == 0 and c < 16:
            for _ in range(c):
                v = r64(v * 10) if n >= 0 else r64(v * r64(F(1, 10)))
        else:
            v = r64(v * (r64(F(10) ** c) if n >= 0 else r64(F(1, 10 ** c))))
    return v

def K10(x, n):
    v = F(x)
    for _ in range(abs(n)):
        v = r64(v * 10) if n >= 0 else r64(v * r64(F(1, 10)))
    return v

def K11(x, n):
    """Low first, but the stored chunk constants for n < 0 are RN64(10^-c) built as RN64(1/RN64(10^c))."""
    v = F(x)
    for c in chunks(n):
        f = r64(F(10) ** c) if n >= 0 else r64(1 / r64(F(10) ** c))
        v = r64(v * f)
    return v

cands = {"K1": K1, "K2": lambda x, n: chain(x, n, "low", "mul"), "K2z": lambda x, n: chain_round(x, n, rz),
         "K7": K7, "K9": K9, "K10": K10, "K11": K11, "K3": lambda x, n: chain(x, n, "low", "div")}
rows = []
with gzip.open(r"C:\Work\DnaCalc\CalcFunc-Evidence\cache\oracle-20430\build-20430\ROUNDDOWN.jsonl.gz", "rt") as f:
    for line in f:
        r = json.loads(line); a = r["args"][0]
        if len(a) != 2 or any(v.get("kind") != "number" for v in a) or r["outcome"]["kind"] != "number":
            continue
        x, digits = a[0]["value"], a[1]["value"]
        if x == 0 or not math.isfinite(x) or digits != int(digits):
            continue
        rows.append((x, int(digits), r["outcome"]["value"]))
print(len(rows), "numeric ROUNDDOWN rows")
near = []
for x, digits, ans in rows:
    e10 = math.floor(math.log10(abs(x)))
    n = 14 - e10
    if abs(n) <= 27:
        continue
    # exact fraction of |x| 10^n
    mant, e = math.frexp(abs(x)); m = int(mant * 2 ** 53); e -= 53
    num, den = (m * 10 ** n, 1) if n >= 0 else (m, 10 ** (-n))
    if e >= 0: num <<= e
    else: den <<= -e
    ip = num // den; frac = F(num - ip * den, den); V = F(num, den)
    dist = (frac - F(1, 2)) / (V / 2 ** 64)
    if abs(dist) <= 3:
        # Excel's decision: the reading went up iff the 15-digit integer is ip + 1, visible when the dropped digits are not read
        # ROUNDDOWN(x, digits) truncates the reading at 10^-digits; the reading's integer is recoverable when digits >= n
        sc = F(10) ** n if n >= 0 else F(1, 10 ** (-n))
        read = F(ans) * sc
        if digits < n:
            continue          # the truncation hides the 15th digit
        up = round(read) == ip + 1
        near.append((x, n, dist, up, ans))
print(len(near), "near-half rows with |n| > 27 and the 15th digit visible")
score = {k: 0 for k in cands}
detail = []
for x, n, dist, up, ans in near:
    preds = {k: decision(f(abs(x), n)) for k, f in cands.items()}
    for k in cands:
        score[k] += preds[k] == up
    if preds["K2"] != up:
        detail.append((h(x), n, float(dist), up, {k: v for k, v in preds.items()}))
print(score)
print("K2 misses:")
for t in detail:
    print("  ", t)
save("t3_trunc_bank_nearhalf.json", {"near_rows": len(near), "score": score, "k2_misses": detail})
