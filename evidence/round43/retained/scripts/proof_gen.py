"""R43 fresh proof set (frozen before asking): ones column stored last and swapped in,
pivot remainders clipped at zero with the earlier column winning ties, and the collinearity
test on vv when columns remain. New seed; LINEST and LOGEST.
G1 zero or tiny third column beside an exact or perturbed pair (all layouts);
G2 three columns with permuted entries (exact pivot ties), with and without a constant;
G3 perturbed pair near the 2^-41 switch with an independent column before or after;
G4 two zero or tiny columns beside a pair; G5 full-rank controls."""
import json, random, struct, sys, math
def hx(x): return '0x%016x' % struct.unpack('<Q', struct.pack('<d', x))[0]
rng = random.Random(4311)
cases = []
def add(tag, cols, const):
    n = len(cols[0])
    cases.append(dict(tag=tag, n=n, cols=cols, const=const))
def pair(n, lo=18, hi=27):
    a = [float(rng.randint(-40, 40)) for _ in range(n)]
    while len(set(a)) < 3: a = [float(rng.randint(-40, 40)) for _ in range(n)]
    m = rng.choice([2.0, 4.0, 0.5, -1.0, 3.0, -2.0])
    b = [m * v for v in a]
    if rng.random() < 0.6:
        i = rng.randrange(n); b[i] += 2.0 ** -rng.randint(lo, hi) * max(1.0, abs(b[i]))
    return a, b
def third(n):
    k = rng.random()
    if k < 0.5: return [0.0] * n
    s = 2.0 ** -rng.choice([200, 300, 600])
    return [s * rng.randint(-9, 9) for _ in range(n)]
for t in range(96):
    n = rng.choice([4, 5, 6, 8, 12, 16, 32]); a, b = pair(n); z = third(n)
    add('G1', rng.choice([[z, a, b], [a, z, b], [a, b, z], [z, b, a], [b, a, z], [b, z, a]]), rng.random() < 0.5)
for t in range(64):
    n = rng.choice([4, 5, 6, 8, 9]); base = [float(rng.randint(-9, 9)) for _ in range(n)]
    while len(set(base)) < 3: base = [float(rng.randint(-9, 9)) for _ in range(n)]
    p1 = base[:]; rng.shuffle(p1); p2 = base[:]; rng.shuffle(p2)
    s = 2.0 ** rng.randint(-12, 12)
    cols = [[s * v for v in base], [s * v for v in p1], [s * v for v in p2]]
    if rng.random() < 0.4: cols = cols[:2]
    rng.shuffle(cols)
    add('G2', cols, rng.random() < 0.6)
for t in range(64):
    n = rng.choice([5, 6, 8, 12, 16]); a, b = pair(n, 19, 23)
    w = [float(rng.randint(-40, 40)) for _ in range(n)]
    add('G3', rng.choice([[w, a, b], [a, w, b], [a, b, w], [b, a, w], [w, b, a]]), rng.random() < 0.5)
for t in range(48):
    n = rng.choice([5, 6, 8, 12]); a, b = pair(n); z1 = third(n); z2 = third(n)
    cols = [a, b, z1, z2]; rng.shuffle(cols)
    add('G4', cols, rng.random() < 0.5)
for t in range(32):
    n = rng.choice([5, 6, 8, 12]); k = rng.choice([2, 3])
    cols = [[float(rng.randint(-30, 30)) for _ in range(n)] for _ in range(k)]
    add('G5', cols, rng.random() < 0.5)
for fn in ['LINEST', 'LOGEST']:
    ws = []
    r2 = random.Random(4312 if fn == 'LINEST' else 4313)
    for ci, c in enumerate(cases):
        n = c['n']; k = len(c['cols'])
        ys = [float(r2.randint(-20, 20)) for _ in range(n)] if fn == 'LINEST' else [float(r2.randint(1, 60)) for _ in range(n)]
        args = [[[hx(v)] for v in ys], [[hx(c['cols'][j][i]) for j in range(k)] for i in range(n)], c['const'], True]
        sels = [(1, cc) for cc in range(1, k + 2)] + [(4, 2), (5, 2), (3, 2), (2, 1)]
        for r, cc in sels:
            ws.append(dict(id='G%d:%d:%d' % (ci, r, cc), tag=c['tag'], args=args, result_index=[r, cc]))
    json.dump(dict(function=fn, witnesses=ws), open(sys.argv[1] + '/%s.proof-questions.json' % fn, 'w'))
    print(fn, len(ws))
