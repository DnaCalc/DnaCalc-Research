import subprocess
p = "/home/ubuntu/projects/DnaCalc/DnaCalc/crates/calc-func/src/ported/special_math_common.rs"
subprocess.run(["git", "-C", "/home/ubuntu/projects/DnaCalc/DnaCalc", "checkout", "--", p], check=True)
s = open(p).read()
old = """    // I_{1-x}(1, a). With a = 1 and x >= 0.1 it is that complement's own arm,
    // 0.5 + (0.5 - exp(b * log y)); just below 0.1 it is not.
    if b == 1.0 && x <= 0.9 {
        let w = crate::ported::excel_numeric::excel_exp(a * crt_log(x));
        return (w, 0.5 + (0.5 - w));
    }
    if a == 1.0 && x >= 0.1 {
        let w1 = crate::ported::excel_numeric::excel_exp(b * crt_log(y));
        return (0.5 + (0.5 - w1), w1);
    }
"""
new = """    // I_{1-x}(1, a). With a = 1 and x >= 0.1 it is that complement's own arm,
    // 0.5 + (0.5 - exp(b * log y)); below 0.1 it is -expm1(b * log1p(-x)) with the
    // product taken on the x87 chain (R43: 3,000 of 3,000 public calls, b from 0.3 to 1,000
    // and x from 1e-9 to 0.0999; the spilled and internal EXPM1 graphs agree on all 20,000
    // sampled arguments, so which of them Excel runs is undetermined).
    if b == 1.0 && x <= 0.9 {
        let w = crate::ported::excel_numeric::excel_exp(a * crt_log(x));
        return (w, 0.5 + (0.5 - w));
    }
    if a == 1.0 && x >= 0.1 {
        let w1 = crate::ported::excel_numeric::excel_exp(b * crt_log(y));
        return (0.5 + (0.5 - w1), w1);
    }
    if a == 1.0 {
        use crate::ported::excel_numeric::{excel_expm1_spilled, excel_log1p_times};
        let w = -excel_expm1_spilled(excel_log1p_times(b, -x));
        return (w, 0.5 + (0.5 - w));
    }
"""
assert old in s, "anchor"
s = s.replace(old, new)
open(p, "w").write(s)
print("ok")
