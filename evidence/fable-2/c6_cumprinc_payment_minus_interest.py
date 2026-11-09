"""C6 (offline on c4 and the big-loan scan): the last period's principal as payment less interest with free constants, by exact
preimage intersection (q = -pv exactly at these loans, so the payment is T = st(pv r)):
 H1 P = st(T - st(pv B)), B free (interest = pv x constant)
 H2 P = st(T - st(st(pv B1) r)), B1 free (interest = r x balance, balance = pv x constant)
 H4 P = st(T - st(T B2)), B2 free (interest = payment x constant)
 H6 P = st(st(pv B1) r) ... no: covered by K fits.  H7 P = st(T - st(T - st(pv B3))) (principal = payment - (payment - X))."""
from common import *
import json


def pre(v):
    v = abs(v); lo = nxt(v, -1); hi = nxt(v, 1)
    return (F(v) + F(lo)) / 2, (F(v) + F(hi)) / 2


def doubles_with_st_eq(target, centre):
    """Doubles I near `centre` with st(I) == target (I itself is a double, so I == target); kept for symmetry."""
    return [target]


def inv_sub(T, P):
    """Doubles I with st(T - I) == P (|I| < T): I in (T - hi, T - lo] roughly; enumerate doubles in that range."""
    lo, hi = pre(P)
    a, b = T - hi, T - lo
    out = []
    v = float(a)
    v = nxt(v, -2)
    for _ in range(600):
        if F(v) > b + (b - a):
            break
        if d(st(T - F(v))) == P:
            out.append(F(v))
        v = nxt(v, 1)
    return out


def inv_mul(I, c):
    """Doubles S with st(S c) == I."""
    out = []
    v = nxt(float(I / c), -4)
    for _ in range(9):
        if d(st(F(v) * c)) == float(I):
            out.append(F(v))
        v = nxt(v, 1)
    return out


def intersect(cur, ivs):
    if cur is None:
        return ivs
    out = []
    for a in cur:
        for b in ivs:
            lo, hi = max(a[0], b[0]), min(a[1], b[1])
            if lo <= hi:
                out.append((lo, hi))
    return out


def fit(samples, r, shape):
    cur = None
    for pv, P in samples:
        pv = F(pv); T = st(pv * F(r))
        ivs = []
        for I in inv_sub(T, P):
            if shape == "H1":
                lo, hi = pre(float(I)); ivs.append((lo / pv, hi / pv))
            elif shape == "H2":
                for S in inv_mul(I, F(r)):
                    lo, hi = pre(float(S)); ivs.append((lo / pv, hi / pv))
            elif shape == "H4":
                lo, hi = pre(float(I)); ivs.append((lo / T, hi / T))
            elif shape == "H7":
                for X in inv_sub(T, float(I)):
                    lo, hi = pre(float(X)); ivs.append((lo / pv, hi / pv))
        if not ivs:
            return None
        cur = intersect(cur, ivs)
        if not cur:
            return None
    return cur


data = json.load(open("c4_cumprinc_last_period_bias.json"))
by = {}
for row in data["rows"]:
    by.setdefault(row["rate"], []).append((fb(row["pv"]), abs(fb(row["cumprinc"]))))
by[0.75] += [(fb(a), abs(fb(c))) for a, b, c in (l.split() for l in open("../ld1/pvscan-big1.txt") if l.strip())]
for r, samples in sorted(by.items()):
    exp = {"H1": F(r) * F(r) / (1 + F(r)), "H2": F(r) / (1 + F(r)), "H4": F(r) / (1 + F(r)), "H7": F(r) / (1 + F(r))}
    line = f"rate {r} ({len(samples)}): "
    for shape in ("H1", "H2", "H4", "H7"):
        iv = fit(samples, r, shape)
        if iv is None:
            line += f"{shape} empty; "
        else:
            lo, hi = iv[0]
            line += f"{shape} K in [{float(lo)!r},{float(hi)!r}] rel {float((lo - exp[shape]) / exp[shape] / F(2) ** -53):+.2f}..{float((hi - exp[shape]) / exp[shape] / F(2) ** -53):+.2f}; "
    print(line)
