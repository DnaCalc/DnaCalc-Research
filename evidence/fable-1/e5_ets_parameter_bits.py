"""E5: what do FORECAST.ETS.STAT's reported alpha/beta/gamma look like at the bit level? float32-representable?
short decimals plus epsilon? The answer constrains the search grid and the parameter storage."""
import random, struct
from common import *

random.seed(51)


def f32_exact(v):
    return struct.unpack(">f", struct.pack(">f", v))[0] == v


def decimal_digits(v):
    for d in range(1, 10):
        if abs(round(v, d) - v) < 1e-15:
            return d
    return None


out = []
for i in range(14):
    n = random.choice([12, 16, 20, 24, 30])
    level, trend = random.uniform(50, 500), random.uniform(-5, 5)
    series = []
    for t in range(n):
        level += trend + random.gauss(0, random.uniform(1, 30))
        series.append(round(level, 3))
    timeline = [float(t + 1) for t in range(n)]
    v = number_array([[s] for s in series])
    tl = number_array([[t] for t in timeline])
    for seas in (0.0, 1.0):
        rows = [[v, tl, float(stat), seas, 1.0, 1.0] for stat in (1, 2, 3, 8)]
        r = ask("FORECAST.ETS.STAT", rows)
        got = vals(r)
        desc = []
        for stat, g in zip((1, 2, 3, 8), got):
            if isinstance(g, float):
                desc.append(f"stat{stat}={g!r} ({h(g)}) f32={f32_exact(g)} dec={decimal_digits(g)}")
            else:
                desc.append(f"stat{stat}={g}")
        print(f"hist {i} n={n} seas={int(seas)}: " + " | ".join(desc))
        out.append({"series": [h(s) for s in series], "seasonality": seas, "stats": [h(g) if isinstance(g, float) else g for g in got]})
save("e5_ets_parameter_bits.json", {"rows": out, "build": meta(r)})
