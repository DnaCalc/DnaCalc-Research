"""I7: IRR with two flows: P(v) = f0 + f1*v is linear, the secant lands in two steps, and the last step is driven by the
rounding residue of P at the previous iterate. That exposes how the linear term is formed:
  M:  f1 * v        Dv: f1 / (1/v)       each with the term stored (Q) or kept in the register (q).
Predictions frozen before asking."""
from common import *
from client import number_array
from i1_irr_objective import A, S, M, D
from i6_irr_finish import solve
def mk(kind, stored):
    kq = rn if stored else (lambda x: x)
    def o(flows, v):
        f0, f1 = flows
        term = r64(f1 * v) if kind == "M" else r64(f1 / rn(r64(1 / v)))
        return rn(r64(f0 + kq(term)))
    return o
C = {"M:Q": mk("M", 1), "M:q": mk("M", 0), "Dv:Q": mk("Dv", 1), "Dv:q": mk("Dv", 0)}
seed = int(sys.argv[1]); rng = random.Random(seed); cases = []
for _ in range(int(sys.argv[2])):
    f0 = -rand_double(rng, 0, 12); ratio = 1 + rng.random() * rng.choice([0.05, 0.3, 1.5])
    cases.append(([f0, -f0 * ratio * (1 + rng.random() * 1e-3)], 0.1))
pred = {k: [solve(o, fl, g, "F1") for fl, g in cases] for k, o in C.items()}
save(f"i7_irr_two_{seed}_frozen.json", {"cases": [[[h(x) for x in fl], h(g)] for fl, g in cases],
                                         "pred": {k: [h(x) if x is not None else None for x in p] for k, p in pred.items()}})
reply = EX.ask("IRR", [[number_array([[x] for x in fl])] for fl, g in cases], priority="interactive"); got = vals(reply)
save(f"i7_irr_two_{seed}.json", {"meta": meta(reply), "excel": [h(g) if isinstance(g, float) else g for g in got]})
nums = [i for i, g in enumerate(got) if isinstance(g, float)]
exact = [i for i in nums if bits(got[i]) == bits(cases[i][0][1] / -cases[i][0][0] - 1)]
print(meta(reply), "numeric", len(nums), "equal to float(f1/-f0 - 1):", len(exact))
for k, p in pred.items():
    print(k, sum(1 for i in nums if p[i] is not None and bits(p[i]) == bits(got[i])), "/", len(nums),
          "distinct from M:Q on", sum(1 for i in nums if p[i] != pred["M:Q"][i]))
