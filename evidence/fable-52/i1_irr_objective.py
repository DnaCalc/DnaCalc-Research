"""I1: IRR. calc-func's solver (secant in v = 1/(1+rate), second point 0.001 away on the side of the sign of the NPV,
stop when |step| < 1e-7 and |npv| < 1e-7, 200 steps), every operation an x87 store, with the objective varied:
  seq:  NPV's loop (N6): p = 1; t += flow / p; p *= base            base = 1 / v
  l2r / r2l: flow / base^k with base^k by binary powering (left-to-right, right-to-left)
and the solver variable varied: v (as landed) or rate (base = 1 + rate, second point rate +- 0.001).
Usage: i1_irr_objective.py SEED N   (predictions are written before Excel is asked)."""
from common import *
from client import number_array
A = lambda a, b: st(a + b); S = lambda a, b: st(a - b); M = lambda a, b: st(a * b); D = lambda a, b: st(a / b)
def pw_l2r(b, k):
    if k == 0: return F(1)
    r = b
    for bit in bin(k)[3:]:
        r = M(r, r)
        if bit == "1": r = M(r, b)
    return r
def pw_r2l(b, k):
    res = F(1); base = b
    while k:
        if k & 1: res = M(res, base)
        k >>= 1
        if k: base = M(base, base)
    return res
def obj(kind, flows, base):
    t = F(0); p = F(1)
    for k, f in enumerate(flows):
        if kind == "seq":
            t = A(t, D(f, p)); p = M(p, base)
        else:
            t = A(t, D(f, (pw_l2r if kind == "l2r" else pw_r2l)(base, k)))
    return t
def solve(kind, var, flows, guess):
    flows = [F(x) for x in flows]; g = F(guess); e3 = F(0.001); tol = F(1e-7)
    if var == "v":
        f = lambda v: obj(kind, flows, D(F(1), v)); x0 = D(F(1), A(F(1), g)); out = lambda v: S(D(F(1), v), F(1))
    else:
        f = lambda r: obj(kind, flows, A(F(1), r)); x0 = g; out = lambda r: r
    y0 = f(x0)
    if y0 == 0: return float(g)
    side = 1 if y0 > 0 else -1
    x1 = A(x0, e3 * side) if var == "v" else A(x0, e3 * side)
    y1 = f(x1)
    for _ in range(200):
        ch = S(y1, y0)
        if ch == 0: return None
        x2 = S(x1, D(M(y1, S(x1, x0)), ch))
        step = abs(S(x2, x1)); x0, y0 = x1, y1; x1 = x2
        if x1 <= 0 and var == "v": return None
        y1 = f(x1)
        if step < tol and abs(y1) < tol: return float(out(x1))
    return None
KEYS = [(k, v) for k in ("seq", "l2r", "r2l") for v in ("v", "rate")]
if __name__ == "__main__":
    seed, n = int(sys.argv[1]), int(sys.argv[2]); rng = random.Random(seed); cases = []
    for _ in range(n):
        k = rng.randint(3, 9)
        flows = [-rand_double(rng, 6, 10)] + [rand_double(rng, 3, 9) * (1 if rng.random() < 0.9 else -1) for _ in range(k - 1)]
        guess = rng.choice([0.1, 0.1, 0.05, 0.2, 0.5, rand_double(rng, -5, -1)])
        cases.append((flows, guess))
    pred = {f"{k}/{v}": [solve(k, v, fl, g) for fl, g in cases] for k, v in KEYS}
    save(f"i1_irr_{seed}_frozen.json", {"cases": [[[h(x) for x in fl], h(g)] for fl, g in cases],
                                         "pred": {k: [h(x) if x is not None else None for x in p] for k, p in pred.items()}})
    reply = EX.ask("IRR", [[number_array([[x] for x in fl]), number(g)] for fl, g in cases], priority="interactive")
    got = vals(reply); nums = [i for i, g in enumerate(got) if isinstance(g, float)]
    save(f"i1_irr_{seed}.json", {"meta": meta(reply), "excel": [h(g) if isinstance(g, float) else g for g in got]})
    print(meta(reply), "cases", n, "numeric", len(nums))
    for key, p in pred.items():
        ok = sum(1 for i in nums if p[i] is not None and bits(p[i]) == bits(got[i]))
        errs = sum(1 for i, g in enumerate(got) if not isinstance(g, float) and p[i] is None)
        d = [abs(ulps(p[i], got[i])) for i in nums if p[i] is not None and p[i] * got[i] > 0]
        print(f"{key:9s} exact {ok}/{len(nums)}  error cells agreed {errs}/{n - len(nums)}  max ulp {max(d) if d else None}")
