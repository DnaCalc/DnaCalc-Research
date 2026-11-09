# research-only switch (never committed): x87 store variants in bpser
p = "/home/ubuntu/projects/DnaCalc/DnaCalc/crates/calc-func/src/ported/special_math_common.rs"
s = open(p).read()
old1 = """    let a0 = a.min(b);
    if a0 >= 1.0 {
        let z = a * crt_log(x) - bratio_betaln(a, b);
        bp = crate::ported::excel_numeric::excel_exp(z) / a;
    } else {"""
new1 = """    let a0 = a.min(b);
    let sw: u32 = std::env::var("R43_BPSER").ok().and_then(|s| s.parse().ok()).unwrap_or(0);
    use crate::ported::excel_numeric::{excel_x87_add as xa, excel_x87_div as xd, excel_x87_mul as xm, excel_x87_sub as xs};
    if a0 >= 1.0 {
        let z = if sw & 1 != 0 { xs(xm(a, crt_log(x)), bratio_betaln(a, b)) } else { a * crt_log(x) - bratio_betaln(a, b) };
        bp = if sw & 8 != 0 { xd(crate::ported::excel_numeric::excel_exp(z), a) } else { crate::ported::excel_numeric::excel_exp(z) / a };
    } else {"""
assert old1 in s
s = s.replace(old1, new1)
old2 = """    loop {
        n += 1.0;
        c = c * (0.5 + (0.5 - b / n)) * x;
        let w = c / (a + n);
        summ += w;
        if w.abs() <= tol {
            break;
        }
    }
    bp * (1.0 + a * summ)
}"""
new2 = """    loop {
        n += 1.0;
        let w;
        if sw & 2 != 0 {
            c = xm(xm(c, xa(0.5, xs(0.5, xd(b, n)))), x);
            w = xd(c, a + n);
            summ = xa(summ, w);
        } else {
            c = c * (0.5 + (0.5 - b / n)) * x;
            w = c / (a + n);
            summ += w;
        }
        if w.abs() <= tol {
            break;
        }
    }
    if sw & 4 != 0 { xm(bp, xa(1.0, xm(a, summ))) } else { bp * (1.0 + a * summ) }
}"""
assert old2 in s
s = s.replace(old2, new2)
open(p, "w").write(s)
print("ok")
