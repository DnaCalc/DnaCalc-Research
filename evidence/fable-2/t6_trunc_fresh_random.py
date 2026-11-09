"""T6: fresh random check of the nibble-chain reading (K2, ties up) on ordinary x at scales with |n| > 27: random mantissas,
exponents 1e-300..1e-30 and 1e+30..1e+300, 240 rows each for ROUNDDOWN(x, n), ROUNDUP(x, n), TRUNC(x, n) with n = 14 - e10.
The prediction is the 15-digit integer N (the reading) recovered from the answer; K1 (one product with the nearest 64-bit
constant) is scored beside K2 to show the separation on unselected inputs."""
from common import *
import math, random
from t1_trunc_constructions import K1, chain
random.seed(33)


def k2(x, n):
    return chain(x, n, "low", "mul")


def reading(f, x, n):
    v = f(x, n); ip = v.numerator // v.denominator
    frac = v - ip
    return ip + (1 if frac >= F(1, 2) else 0)


xs = []
for _ in range(240):
    e = random.choice([random.randint(-300, -30), random.randint(30, 300)])
    xs.append(random.uniform(1.0, 9.999) * 10.0 ** e)
rows = []
for x in xs:
    e10 = math.floor(math.log10(abs(x)))
    n = 14 - e10
    rows.append((x, n))
score = {"ROUNDDOWN": [0, 0, 0], "ROUNDUP": [0, 0, 0], "TRUNC": [0, 0, 0]}
recs = []
for fn in score:
    ans = vals(ask(fn, [[x, float(n)] for x, n in rows]))
    for (x, n), a in zip(rows, ans):
        if not isinstance(a, float):
            continue
        sc = F(10) ** n if n >= 0 else F(1, 10 ** (-n))
        N_excel = round(F(a) * sc)
        exact = F(x) * sc; ip = exact.numerator // exact.denominator
        # ROUNDUP adds a unit when the dropped part was not zero; the reading itself is N; for ROUNDUP the published N is
        # N + (dropped digits != 0) where the dropped part is what lies below the 15th digit of the reading: compare readings
        # through ROUNDDOWN/TRUNC only; for ROUNDUP compare N_excel - 1 when the exact x is not the reading.
        n2 = reading(k2, x, n); n1 = reading(K1, x, n)
        if fn == "ROUNDUP":
            n2u = n2 + (1 if F(x) * sc != n2 else 0); n1u = n1 + (1 if F(x) * sc != n1 else 0)
            score[fn][0] += 1; score[fn][1] += N_excel == n2u; score[fn][2] += N_excel == n1u
        else:
            score[fn][0] += 1; score[fn][1] += N_excel == n2; score[fn][2] += N_excel == n1
        recs.append({"fn": fn, "x": h(x), "n": n, "answer": h(a), "k2_reading": str(n2), "k1_reading": str(n1), "excel_reading": str(N_excel)})
for fn, (n, s2, s1) in score.items():
    print(f"{fn:10} rows {n}: K2 {s2}  K1 {s1}")
save("t6_trunc_fresh_random.json", {"rows": recs, "build": meta(ask("TRUNC", [[1.5, 0.0]]))})
