"""R43 fresh proof for the a = 1, x < 0.1 arm of the regularized beta (frozen before asking).
usage: beta_proof.py freeze | ask"""
import sys, os, json, struct, math, hashlib, time
import numpy as np
sys.path.insert(0, "/home/ubuntu/projects/DnaCalc/DnaCalc/tools/excel-service")
sys.path.insert(0, "/home/ubuntu/codex-runs/r43")
import racer
QF = "/home/ubuntu/codex-runs/r43/q/BETA.a1-proof.json"
def h(x): return "%016x" % struct.unpack("<Q", struct.pack("<d", float(x)))[0]
if sys.argv[1] == "freeze":
    rng = np.random.default_rng(4340)
    Q = []
    for _ in range(2000):
        x = float(np.exp(rng.uniform(math.log(1e-300), math.log(0.0999999))))
        b = float(np.exp(rng.uniform(math.log(1e-3), math.log(1e6))))
        Q.append(("H1", "BETA.DIST", [x, 1.0, b, "T"]))
    x = 0.1
    for k in range(100):
        x = float(np.nextafter(x, 0)); b = float(np.exp(rng.uniform(math.log(0.1), math.log(1e4))))
        Q.append(("H2", "BETA.DIST", [x, 1.0, b, "T"]))
    for _ in range(500):
        n = float(int(np.exp(rng.uniform(math.log(2), math.log(1e6)))))
        p = 1.0 - float(np.exp(rng.uniform(math.log(1e-12), math.log(0.1))))
        Q.append(("H3", "BINOM.DIST", [n - 1.0, n, p, "T"]))
    for _ in range(300):
        x = float(np.exp(rng.uniform(math.log(1e-30), math.log(0.0999))))
        b = float(np.exp(rng.uniform(math.log(1e-2), math.log(1e5))))
        Q.append(("H4", "BETADIST", [x, 1.0, b]))
    preds = {}
    for name, binary, env in [("base", "racer-beta", "0"), ("final", "racer-final", "")]:
        racer.RACER = os.path.expanduser("~/codex-runs/r43/" + binary)
        os.environ["R43_BETA_A1"] = env
        out = racer.invoke([(fn, args) for tag, fn, args in Q])
        preds[name] = [o if isinstance(o, str) else h(o) for o in out]
    json.dump(dict(questions=Q, predictions=preds), open(QF, "w"))
    print(len(Q), "separating", sum(a != b for a, b in zip(preds["base"], preds["final"])),
          hashlib.sha256(open(QF, "rb").read()).hexdigest())
else:
    from client import Client, number, boolean
    d = json.load(open(QF)); Q = d["questions"]
    c = Client(timeout=300)
    def arg(v): return boolean(True) if v == "T" else number(v)
    answers = [None] * len(Q); t0 = time.time()
    for fn in ["BETA.DIST", "BINOM.DIST", "BETADIST"]:
        idx = [i for i, q in enumerate(Q) if q[1] == fn]
        r = c.ask(fn, [[arg(v) for v in Q[i][2]] for i in idx], priority="bulk")
        for i, row in zip(idx, r["rows"]):
            answers[i] = row.get("bits") if row.get("type") == "number" else "error:" + str(row.get("code"))
    dt = time.time() - t0
    json.dump(dict(answers=answers, seconds=dt, build=r["excel_build"], env=r["env_id"]),
              open("/home/ubuntu/codex-runs/r43/fresh/BETA.a1-proof.replies.json", "w"))
    import collections
    c2 = collections.Counter()
    for q, a, pb, pf in zip(Q, answers, d["predictions"]["base"], d["predictions"]["final"]):
        c2[q[0] + ":n"] += 1; c2[q[0] + ":base"] += a == pb; c2[q[0] + ":final"] += a == pf
    print("asked in %.2fs" % dt, dict(sorted(c2.items())))
