"""N6: the NPV claim, in plain binary64 (Python floats are IEEE binary64, one rounding per operation):
    b = 1 + rate;  p = 1;  t = 0;  for v in values: t += v / p; p *= b;   NPV = t / b
i.e. the flows are discounted as if the first one were at time 0, and the sum is divided by (1 + rate) once at the end.
Modes: replay FILE (offline, on earlier answers) | proof SEED N (predictions frozen to disk before Excel is asked)."""
from common import *
from n1_npv_shapes import inputs


def npv(rate, values):
    b = 1.0 + rate; p = 1.0; t = 0.0
    for v in values:
        t += v / p; p *= b
    return t / b


def edge_inputs(seed, n):
    rng = random.Random(seed); rows = []
    for _ in range(n):
        kind = rng.randrange(5)
        r = [rand_double(rng, -30, -10), -rand_double(rng, -1, -1) , rand_double(rng, 5, 40), -1 - rand_double(rng, -3, 3),
             rng.choice([0.05, 0.1, 0.08, 0.125, 0.07, 1e-3, 0.0])][kind]
        k = rng.randint(1, 30)
        rows.append([r] + [rand_double(rng, -10, 30, sign=True) if rng.random() < 0.9 else 0.0 for _ in range(k)])
    return rows


if __name__ == "__main__":
    if sys.argv[1] == "replay":
        d = json.load(open(os.path.join(HERE, sys.argv[2])))
        f = lambda s: struct.unpack(">d", bytes.fromhex(s))[0]
        ok = sum(bits(npv(f(r[0]), [f(x) for x in r[1:]])) == bits(f(e)) for r, e in zip(d["rows"], d["excel"]))
        print("replay", sys.argv[2], ok, "/", len(d["rows"]))
    else:
        seed, n = int(sys.argv[2]), int(sys.argv[3])
        rows = inputs(seed, n) + edge_inputs(seed + 1, n)
        pred = [npv(r[0], r[1:]) for r in rows]
        frozen = {"claim": __doc__, "seed": seed, "rows": [[h(x) for x in r] for r in rows], "predicted": [h(p) for p in pred]}
        save(f"n6_npv_proof_{seed}_frozen.json", frozen)                      # written before the question is sent
        reply = ask("NPV", rows); got = vals(reply)
        miss = [i for i, g in enumerate(got) if not (isinstance(g, float) and bits(g) == bits(pred[i]))]
        frozen.update(meta=meta(reply), excel=[h(g) if isinstance(g, float) else g for g in got], misses=miss)
        save(f"n6_npv_proof_{seed}.json", frozen)
        print(meta(reply)); print("fresh proof:", len(rows) - len(miss), "/", len(rows))
        for i in miss[:12]:
            g = got[i]
            print(i, len(rows[i]) - 1, h(rows[i][0]), g if not isinstance(g, float) else ulps(pred[i], g), h(pred[i]))
