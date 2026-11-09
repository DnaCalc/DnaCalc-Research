p = "/home/ubuntu/projects/DnaCalc/DnaCalc/crates/calc-func/src/ported/regression_forecast_family.rs"
s = open(p).read()
old = """        let collinear = scale[j].1 > 0.0 && matrix[j].0.is_some() && {
            let before = scale[j].1;
            let remainder = matrix[j].1[..j].iter().fold(before, |rest, r| rest - r * r);
            remainder < before * COLLINEAR_REMAINDER
        };"""
new = """        let collinear = scale[j].1 > 0.0 && matrix[j].0.is_some() && {
            let before = scale[j].1;
            if j + 1 == matrix.len() {
                let remainder = matrix[j].1[..j].iter().fold(before, |rest, r| rest - r * r);
                remainder < before * COLLINEAR_REMAINDER
            } else {
                // With columns still to come the test reads the reflector's normaliser
                // `vv = rest + (x0 - alpha)^2` (R43: near-collinear pairs beside a zero or tiny
                // third column switch at vv = 2^-41 times the initial sum of squares, n = 4 to
                // 32, with and without a constant; the last column switches on its remainder).
                let norm = squares.sqrt();
                let alpha = if x0 > 0.0 { -norm } else { norm };
                let v0 = excel_numeric::excel_x87_sub(x0, alpha);
                let vv = excel_numeric::excel_x87_add(rest, excel_numeric::excel_x87_mul(v0, v0));
                vv < before * COLLINEAR_REMAINDER
            }
        };"""
assert old in s
s = s.replace(old, new)
open(p, "w").write(s)
print("ok")
