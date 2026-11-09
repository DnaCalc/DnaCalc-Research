# research-only switch (never committed): BINOM CDF routes for p < 0.5
p = "/home/ubuntu/projects/DnaCalc/DnaCalc/crates/calc-func/src/ported/discrete_dist.rs"
s = open(p).read()
old = """    } else {
        crate::ported::special_math_common::regularized_beta(q, n - k, k + 1.0)
    }
}"""
new = """    } else {
        let v: u32 = std::env::var("R43_BINOM").ok().and_then(|s| s.parse().ok()).unwrap_or(0);
        match v {
            1 => 1.0 - crate::ported::special_math_common::regularized_beta(p, k + 1.0, n - k),
            2 => crate::ported::special_math_common::bratio(k + 1.0, n - k, p, q).1,
            3 => crate::ported::special_math_common::bratio(n - k, k + 1.0, q, p).0,
            4 => 0.5 + (0.5 - crate::ported::special_math_common::regularized_beta(p, k + 1.0, n - k)),
            _ => crate::ported::special_math_common::regularized_beta(q, n - k, k + 1.0),
        }
    }
}"""
assert old in s, "anchor"
s = s.replace(old, new)
open(p, "w").write(s)
print("ok")
