"""N8: offline. Same NPV shape with every operation done in an x87 register and stored (double rounding), against the
plain binary64 claim, on all NPV answers so far. Separates SSE-style arithmetic from x87 stores for this kernel."""
from common import *
f = lambda s: struct.unpack(">d", bytes.fromhex(s))[0]
def npv_st(rate, values):
    b = st(1 + F(rate)); p = F(1); t = F(0)
    for v in values:
        t = st(t + st(F(v) / p)); p = st(p * b)
    return float(st(t / b))
def npv_plain(rate, values):
    b = 1.0 + rate; p = 1.0; t = 0.0
    for v in values:
        t += v / p; p *= b
    return t / b
tot = plain = dbl = differ = 0
for name in sys.argv[1:]:
    d = json.load(open(os.path.join(HERE, name)))
    for r, e in zip(d["rows"], d["excel"]):
        if not (isinstance(e, str) and len(e) == 16): continue
        r = [f(x) for x in r]; e = f(e)
        if abs(r[0]) > 1e30 or len(r) > 12: continue        # keep the exact-rational replay cheap and in range
        a, b_ = npv_plain(r[0], r[1:]), npv_st(r[0], r[1:])
        tot += 1; plain += bits(a) == bits(e); dbl += bits(b_) == bits(e); differ += bits(a) != bits(b_)
print("rows", tot, "plain binary64", plain, "x87 stored", dbl, "rows where the two differ", differ)
