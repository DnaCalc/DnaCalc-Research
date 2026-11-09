"""P19: for a <= m, the cumulative's prefactor as the mass's FIRST route (exp(a LN m - m - GAMMALN(a+1)), compensated stores, or
plain) even where the public mass takes the Loader route, with the asymptotic, continued-fraction and finite-series tails."""
from common import *
import json, collections
from p9_poisson_series_and_fixed_d import comp_chain
from p15_poisson_cumulative_from_mass_k import nswc_summ
from p14_poisson_cf import gamma_fraction
rows = json.load(open("p16_poisson_cumulative_prefactor.json"))["rows"]
G = dict(zip(range(2, 21), vals(ask("GAMMALN", [[float(a)] for a in range(2, 21)]))))
sel = [r for r in rows if fb(r["m"]) >= r["k"] + 1]
args = []
for r in sel:
    a = r["k"] + 1; m = fb(r["m"]); l = fb(r["ln"])
    t = st(F(a) * F(l)); n1 = st(t - F(m)); n2 = st(n1 - F(G[a]))
    args.append(d(n2))
EX = dict(zip(sorted(set(args)), vals(ask("EXP", [[e] for e in sorted(set(args))]))))
by = collections.defaultdict(collections.Counter); tot = collections.Counter()
for r, e2 in zip(sel, args):
    a = r["k"] + 1; m = fb(r["m"]); l = fb(r["ln"]); q = fb(r["cum"])
    reg = "a<=m<31" if m < 31 else "m>=31"
    tot[reg] += 1
    pf7 = F(comp_chain(d(st(F(a) * F(l))), [m, G[a]], lambda t: EX[e2]))      # compensated first route
    pf8 = F(EX[e2])                                                            # plain stored exponent
    summ = nswc_summ(a, m); cf = gamma_fraction(a, m)
    for n_, pf in (("PF7", pf7), ("PF8", pf8)):
        ax = st(F(a) / F(m))
        cand = {"asy": st(st(pf * ax) * F(1.0 + summ)), "asy regs": st(rn(pf * rn(F(a) / F(m), 64), 64) * F(1.0 + summ)),
                "cf": st(st(pf * a) * cf), "cf regs": st(rn(pf * a, 64) * cf), "asy /m": st(st(st(pf * a) / F(m)) * F(1.0 + summ))}
        for cn, v in cand.items():
            by[reg][n_ + " " + cn] += d(v) == q
for reg in sorted(tot):
    print(f"{reg:8} n={tot[reg]:3}: " + "  ".join(f"{k} {v}" for k, v in sorted(by[reg].items(), key=lambda t: -t[1])))
save("p19_poisson_first_route_prefactor.json", {"gammaln_2_20": {str(a): h(v) for a, v in G.items()}, "exp_args": {h(k): h(v) for k, v in EX.items()}, "build": meta(ask("LN", [[2.0]]))})
