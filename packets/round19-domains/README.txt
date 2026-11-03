R19 adaptive Bessel domain packet. Python 3 and Cargo are required on Windows.
Run from any directory, after pulling the published DnaCalc master:
  python domain_bisection.py run --checkout C:\Work\DnaCalc\DnaCalc-win --evidence C:\Work\DnaCalc\OxFunc-Evidence --output C:\Work\DnaCalc\codex-runs\excel-results-r19-domains

Use a new output directory outside both checkouts and this packet. The script
uses only cargo xtask witness read --inputs; no --new-only and no new oracle.
It rechecks all endpoints, bisects positive binary64 magnitude coordinates
independently per function/order/sign, then asks both endpoints and two exterior
neighbors. I's adjacent 709.7 boundary is separately rechecked through order10000.
Initial controls also surround an unverified2^27-1 private-phase guard; only
their observed answers can shorten a bracket. The full bisection remains valid
if that hypothesis is wrong. No predicted answer is submitted to the reader.
Every requested/submitted/Value2 operand and raw result must agree before a
reading changes a bracket. Predictions never count as answers. Missing results,
unexpected errors, operand changes, conflicting repeated results and observed
nonmonotonicity abort. A batch failure preserves the partial attempt.

The bisection assumes one numeric-to-NUM transition inside each starting bracket;
the packet cannot prove global monotonicity or identify Excel's hidden program.
Initial high endpoints lacking an R16 NUM answer expand by powers of two only
if Excel returns a number. A completed report identifies adjacent local endpoints.
Return the whole output directory: packet, requests, canonical WitnessSets,
worker requests/results, frozen predictions, step journal and file hashes.
The original canonical readings also stay in OxFunc-Evidence/readings.
Do not submit plan.json, summaries or predictions to the witness reader.
