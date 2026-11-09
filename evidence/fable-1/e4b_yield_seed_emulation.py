"""E4b (offline): qualitative replay of YIELD candidates on the round-46 zero-coupon cash-scale map. A candidate must
reproduce the regimes (exact root at tiny cash, deficit linear in pr at small cash, linear in 1/pr at large cash,
the erratic crossover near face amounts, the negative root for pr=80) before any bit comparison."""
import gzip, json, struct
from common import from_bits

fz = json.load(gzip.open("../round46/retained/yield-map-frozen.json.gz"))["cases"]
an = json.load(gzip.open("../round46/retained/yield-map-result.json.gz"))["answers"]
excel = {}
for c, a in zip(fz, an):
    if a.get("type") == "number":
        excel[(c["n"], c["frequency"], c["base"], c["exponent"])] = from_bits(a["bits"])


def price(y, red, N, f):
    return red / (1 + y / f) ** N


def dprice(y, red, N, f):
    return -N / f * red / (1 + y / f) ** (N + 1)


def newton(y0, pr, red, N, f, tol=1e-10, before=True, fd=None, maxit=100):
    y = y0
    for _ in range(maxit):
        P = price(y, red, N, f)
        if fd is None:
            d = dprice(y, red, N, f)
        else:
            d = (price(y + fd, red, N, f) - P) / fd
        if d == 0:
            return None
        step = (P - pr) / d
        if before and abs(step) < tol:
            return y
        y = y - step
        if not before and abs(step) < tol:
            return y
    return y


seeds = {
    "textbook par100: (C+(100-P)/N)/((100+P)/2)": lambda pr, red, N, f: ((100 - pr) / N) / ((100 + pr) / 2),
    "textbook red:    (C+(R-P)/N)/((R+P)/2)": lambda pr, red, N, f: ((red - pr) / N) / ((red + pr) / 2),
    "par100 in years: (C+(100-P)/(N/f))/((100+P)/2)": lambda pr, red, N, f: ((100 - pr) / (N / f)) / ((100 + pr) / 2),
    "constant 0.1": lambda pr, red, N, f: 0.1,
    "(100-P)/(N*P) capped at 1": lambda pr, red, N, f: min(1.0, (100 - pr) / (N * pr)) if pr > 0 else 1.0,
}
designs = [(2, 1, 2.0), (4, 1, 2.0), (2, 2, 2.0), (4, 2, 2.0), (2, 1, 1.5)]
for name, seed in seeds.items():
    print("==", name)
    for (n, f, base) in designs:
        N = n * f
        root = f * (base - 1)
        pr0, red0 = 100.0 / base ** N, 100.0
        line = []
        for e in (-300, -100, -20, -8, -4, 0, 4, 8, 12, 16, 20, 40, 300):
            s = 2.0 ** e
            try:
                y = newton(seed(pr0 * s, red0 * s, N, f), pr0 * s, red0 * s, N, f)
            except (ZeroDivisionError, OverflowError):
                y = None
            ex = excel.get((n, f, base, e))
            line.append(f"e{e:+d}:{(y - root) if y is not None else float('nan'):+.1e}/{(ex - root):+.1e}")
        print(f"  n{n} f{f} base{base}: " + " ".join(line))
