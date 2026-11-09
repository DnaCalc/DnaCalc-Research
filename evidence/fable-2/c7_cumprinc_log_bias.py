"""C7 (offline on c4): is the rate-dependent bias of the last-period constant the rounding of the stored logarithm L = LN(1+r)
(a chain e^-L in the registers has relative error -(L - ln(1+r)))? Units of 2^-53."""
from common import *
from decimal import Decimal, getcontext
import json
getcontext().prec = 60
data = json.load(open("c4_cumprinc_last_period_bias.json"))
meas = {0.75: -2.26, 0.5: -0.40, 0.25: +0.18, 0.1: -0.96, 0.05: +1.16, 1.0: 0.0, 3.0: 0.0, 1.5: +1.78, 0.3: +0.10}
print("rate   measured C_eff dev   -(L - ln(1+r))   u=EXP(-L) dev   RN(1/(1+r)) dev    (units 2^-53, relative)")
for rs, lh in data["ln_1_plus_r"].items():
    r = float(rs); L = fb(lh); u = fb(data["exp_neg_ln"][rs])
    dl = -(F(L) - F(Decimal(1 + F(r)).ln() if False else Decimal(1) + Decimal(r) if False else F(Decimal(1 + F(r)).ln()) if False else F((Decimal(1) + Decimal(r)).ln())))
    v = 1 / (1 + F(r))
    print(f"{r:5}  {meas[r]:+6.2f}             {float(dl / F(2) ** -53):+7.3f}        {float((F(u) - v) / v / F(2) ** -53):+7.3f}        {float((rn(v) - v) / v / F(2) ** -53):+7.3f}")
