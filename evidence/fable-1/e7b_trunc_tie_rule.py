"""E7b: the tie rule of the fifteen-digit reading, isolated from scaling (n = 0: no product). x = 16 significant
digits ending in .5, exactly representable. ROUNDDOWN(x,0) returns the reading's integer; the 15th digit's parity
separates half-up from half-even."""
from common import *
xs = [123456789012344.5, 123456789012345.5, 123456789012346.5, 123456789012347.5, 999999999999998.5, 100000000000000.5]
for fn in ("ROUNDDOWN", "TRUNC", "ROUNDUP", "ROUND"):
    v = vals(ask(fn, [[x, 0.0] for x in xs]))
    print(fn, [f"{x!r}->{a!r}" for x, a in zip(xs, v)])
save("e7b_trunc_tie_rule.json", {"x": [h(x) for x in xs],
     "answers": {fn: [h(a) for a in vals(ask(fn, [[x, 0.0] for x in xs]))] for fn in ("ROUNDDOWN", "TRUNC", "ROUNDUP", "ROUND")}})
