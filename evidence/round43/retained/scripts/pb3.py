p = "/home/ubuntu/projects/DnaCalc/DnaCalc/crates/calc-func/src/ported/discrete_dist.rs"
s = open(p).read()
old = """fn binom_cumulative(k: f64, n: f64, p: f64, q: f64) -> f64 {
    if k == n {
        1.0
    } else if k == 0.0 {
        binom_mass(k, n, p, q)
    } else {"""
new = """fn binom_cumulative(k: f64, n: f64, p: f64, q: f64) -> f64 {
    if k == n {
        1.0
    } else if k == 0.0 {
        binom_mass(k, n, p, q)
    } else if n - k == 1.0 && q < 0.1 {
        // One short of n with q below 0.1: -expm1(n * log1p(-q)) with the product stored in
        // binary64, not BETA.DIST's x87 chain at a = 1 (R43: Excel's BETA.DIST(q, 1, n) and
        // BINOM.DIST(n - 1, n, p, TRUE) differ on 91 of 500 fresh calls; this matches 499 of
        // the binomial's). Which EXPM1 graph runs is undetermined: none of the asked
        // arguments separates the internal and spilled ones.
        use crate::ported::excel_numeric::{excel_expm1_internal, excel_log1p};
        -excel_expm1_internal(n * excel_log1p(-q))
    } else {"""
assert old in s
s = s.replace(old, new)
open(p, "w").write(s)
print("ok")
