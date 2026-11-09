fable-1: cheap probes behind design/fresh-ideas-fable-1.md
Claude Fable 5.1, Windows GOVERT-OMEN, 2026-11-09. Excel 16.0.20527 (file 16.0.20527.20150),
env f3269992527c, x64 EXCEL.EXE, through tools/excel-service/client.py (touchstone socket,
priority interactive). Each e*.py writes the e*.json beside it: inputs and Excel's bits as hex,
plus the build record (meta). common.py is the helper. No application code changed.

e1   GAMMALN vs LN(GAMMA), EXP(GAMMALN) vs GAMMA, one-step stored recurrences up/down on [0.7,8);
     bisection of the upward-recurrence join: last holds 3fe6666666666664, first fails 3fe6666666666666.
e1b  [3,8): GAMMALN(x) = RN(GAMMALN(x-k) + LN(product from the lowest factor)) 600/600; GAMMA forms.
e1c  core [0.7,3) error against math.lgamma (double truth; mpmath absent).
e1d  core is not RN(ln(binary64 G)): 74/558 on grid, 75.5 expected by chance.
e1e  offline: absolute error of the core binned in x (flat, 2-5e-16, not vanishing at the roots).
e2   NORM.S.INV mirror by region (exact except (0.925, 15/16]), windows at 1/16 and 15/16,
     NORM.INV(p,m,s) == RN(m + RN(s*NORM.S.INV(p))) 400/400.
e3   BETA.DIST / GAMMA.DIST cdf vs pdf*x/a at tiny x; pdf vs public-leaf exponent forms.
e3b  dyadic scaling of BETA.DIST pdf and cdf for integer a (both exponent calculations).
e4   YIELD at joint cash scales 2^-500..2^300, two-period bonds (seven designs).
e4b  offline: regime replay of candidate seeds on the round-46 exponent map (no match; harness).
e5   FORECAST.ETS.STAT alpha/beta/gamma/step bits for 14 histories, seasonality 0 and 1.
e6   offline: exact interval intersection over a hidden constant for CUMPRINC two-stage shapes,
     on evidence/ld1/pvscan-big1.txt and pvscan-40.txt (all six shapes empty).
e7   ROUNDDOWN at designed near-half fifteen-digit readings, four scale classes (threshold at 0,
     bias about +1 unit at n=314, about -0.3 unit at n=-286).
e7b  tie rule without scaling: directed functions read the half up; ROUND reads …344.5 as …344.
