"""R43 fresh proof for BINOM.DIST's k = n - 1 arm (frozen before asking).
usage: binom_proof.py freeze | ask"""
import sys, os, json, struct, math, hashlib, time, collections
import numpy as np
sys.path.insert(0, "/home/ubuntu/projects/DnaCalc/DnaCalc/tools/excel-service")
sys.path.insert(0, "/home/ubuntu/codex-runs/r43")
import racer
QF = "/home/ubuntu/codex-runs/r43/q/BINOM.n1-proof.json"
def h(x): return "%016x" % struct.unpack("<Q", struct.pack("<d", float(x)))[0]
if sys.argv[1] == "freeze":
    rng = np.random.default_rng(4350)
    Q = []
    for _ in range(1000):
        n = float(int(np.exp(rng.uniform(math.log(2), math.log(1e7)))))
        p = 1.0 - float(np.exp(rng.uniform(math.log(1e-14), math.log(0.1))))
        Q.append(("H5", "BINOM.DIST", [n - 1.0, n, p, "T"]))
    for _ in range(300):
        n = float(int(np.exp(rng.uniform(math.log(2), math.log(1e5)))))
        p = float(rng.uniform(0.5, 0.9))
        Q.append(("H6", "BINOM.DIST", [n - 1.0, n, p, "T"]))
    for _ in range(200):
        n = float(int(np.exp(rng.uniform(math.log(2), math.log(1e6)))))
        p = 1.0 - float(np.exp(rng.uniform(math.log(1e-12), math.log(0.1))))
        Q.append(("H7", "BINOM.DIST.RANGE", [n, p, 0.0, n - 1.0]))
    p = 0.9
    for k in range(50):
        p = float(np.nextafter(p, 2)); n = float(int(np.exp(rng.uniform(math.log(2), math.log(1e4)))))
        Q.append(("H8", "BINOM.DIST", [n - 1.0, n, p, "T"]))
    p = 0.9
    for k in range(50):
        p = float(np.nextafter(p, 0)); n = float(int(np.exp(rng.uniform(math.log(2), math.log(1e4)))))
        Q.append(("H8", "BINOM.DIST", [n - 1.0, n, p, "T"]))
    preds = {}
    for name, binary in [("base", "racer-base2"), ("final", "racer-final2")]:
        racer.RACER = os.path.expanduser("~/codex-runs/r43/" + binary)
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
    for fn in ["BINOM.DIST", "BINOM.DIST.RANGE"]:
        idx = [i for i, q in enumerate(Q) if q[1] == fn]
        r = c.ask(fn, [[arg(v) for v in Q[i][2]] for i in idx], priority="bulk")
        for i, row in zip(idx, r["rows"]):
            answers[i] = row.get("bits") if row.get("type") == "number" else "error:" + str(row.get("code"))
    dt = time.time() - t0
    json.dump(dict(answers=answers, seconds=dt, build=r["excel_build"], env=r["env_id"]),
              open("/home/ubuntu/codex-runs/r43/fresh/BINOM.n1-proof.replies.json", "w"))
    c2 = collections.Counter()
    for q, a, pb, pf in zip(Q, answers, d["predictions"]["base"], d["predictions"]["final"]):
        c2[q[0] + ":n"] += 1; c2[q[0] + ":base"] += a == pb; c2[q[0] + ":final"] += a == pf
    print("asked in %.2fs" % dt, dict(sorted(c2.items())))
