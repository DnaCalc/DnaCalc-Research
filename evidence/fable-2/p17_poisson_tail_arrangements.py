"""P17 (offline on p16 answers): for a <= m (both the a <= m < 31 and the m >= 31 regions) the cumulative as the public mass at
a = k + 1 (M) times the finite backward series (the NSWC asymptotic form; exact for integer a), in several arrangements:
 A1 st(st(M st(a/m)) (1+summ_nswc))     A2 st(M (a/m) (1+summ)) registers, one store   A3 st(st(st(M a)/m) (1+summ))
 A4 st(st(st(M/m) a) (1+summ))           A5 st(M st(st(a/m) (1+summ)))                   A6 A1 with a plain forward binary64 sum
 A7 A1 with the Horner sum               A8 st(st(M st(a/m)) st(1 + summ_x87)) with the series summed in PC64 stores
 A9 st(M st(a/m)) + st(st(M st(a/m)) summ) stored     A10 st(st(M a) st((1+summ)/m))."""
from common import *
import json, collections
from p15_poisson_cumulative_from_mass_k import nswc_summ


def r64(q):
    return rn(q, 64)


def series_variants(a, m):
    s_nswc = nswc_summ(a, m)
    t = 1.0; s_fwd = 0.0
    for i in range(1, a):
        t = t * (a - i) / m; s_fwd += t
    hh = 0.0
    for j in range(a - 1, 0, -1):
        hh = (1.0 + hh) * j / m
    t = F(1); s_x87 = F(0)
    for i in range(1, a):
        t = st(st(t * (a - i)) / F(m)); s_x87 = st(s_x87 + t)
    return s_nswc, s_fwd, hh, d(s_x87)


rows = json.load(open("p16_poisson_cumulative_prefactor.json"))["rows"]
by = collections.defaultdict(collections.Counter); tot = collections.Counter()
for r in rows:
    k = r["k"]; a = k + 1; m = fb(r["m"]); q = fb(r["cum"]); M = F(fb(r["mass_a"]))
    if m < a:
        continue
    reg = "a<=m<31" if m < 31 else "m>=31"
    tot[reg] += 1
    s1, s2, s3, s4 = series_variants(a, m)
    ax = st(F(a) / F(m))
    cand = {
        "A1": st(st(M * ax) * F(1.0 + s1)),
        "A2": st(r64(M * r64(F(a) / F(m))) * F(1.0 + s1)),
        "A3": st(st(st(M * a) / F(m)) * F(1.0 + s1)),
        "A4": st(st(st(M / F(m)) * a) * F(1.0 + s1)),
        "A5": st(M * st(ax * F(1.0 + s1))),
        "A6": st(st(M * ax) * F(1.0 + s2)),
        "A7": st(st(M * ax) * F(1.0 + s3)),
        "A8": st(st(M * ax) * st(1 + F(s4))),
        "A9": st(st(M * ax) + st(st(M * ax) * F(s1))),
        "A10": st(st(M * a) * st(F(1.0 + s1) / F(m))),
        "A11": st(st(M * ax) * F(1.0 + s1)) if False else st(r64(st(M * ax) * F(1.0 + s1))),
        "A12": st(st(M * ax) * (1 + F(s2))),
        "A13": st(r64(r64(M * r64(F(a) / F(m))) * (1 + F(s2)))),
    }
    for n_, v in cand.items():
        by[reg][n_] += d(v) == q
for reg in sorted(tot):
    print(f"{reg:8} n={tot[reg]:3}: " + "  ".join(f"{k} {v}" for k, v in sorted(by[reg].items(), key=lambda t: -t[1])))
