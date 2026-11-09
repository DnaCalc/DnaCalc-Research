# research-only switch (not committed): a = 1, x < 0.1 candidates for BETA.DIST
p = "/home/ubuntu/projects/DnaCalc/DnaCalc/crates/calc-func/src/ported/special_math_common.rs"
s = open(p).read()
old = """    if a == 1.0 && x >= 0.1 {
        let w1 = crate::ported::excel_numeric::excel_exp(b * crt_log(y));
        return (0.5 + (0.5 - w1), w1);
    }
"""
new = old + """    if a == 1.0 && x < 0.1 {
        use crate::ported::excel_numeric as en;
        let v: u32 = std::env::var("R43_BETA_A1").ok().and_then(|s| s.parse().ok()).unwrap_or(0);
        let w = match v {
            1 => -en::excel_expm1(b * en::excel_log1p(-x)),
            2 => -en::excel_expm1(en::excel_x87_mul(b, en::excel_log1p(-x))),
            3 => -en::excel_expm1(en::excel_log1p_times(b, -x)),
            4 => -en::excel_expm1_internal(b * en::excel_log1p(-x)),
            5 => -en::excel_expm1_spilled(b * en::excel_log1p(-x)),
            6 => -en::excel_expm1_spilled(en::excel_log1p_times(b, -x)),
            7 => -en::excel_expm1_internal(en::excel_log1p_times(b, -x)),
            _ => f64::NAN,
        };
        if v != 0 {
            return (w, 0.5 + (0.5 - w));
        }
    }
"""
assert old in s
s = s.replace(old, new)
open(p, "w").write(s)
print("ok")
