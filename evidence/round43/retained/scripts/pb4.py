p1 = "/home/ubuntu/projects/DnaCalc/DnaCalc/crates/calc-func/src/ported/beta_gamma_stats_family.rs"
s = open(p1).read()
anchor = """    #[cfg_attr(not(target_arch = "wasm32"), test)]
    #[cfg_attr(target_arch = "wasm32", wasm_bindgen_test::wasm_bindgen_test)]
    fn gamma_dist_half_shape_scale_two_matches_erf_sqrt_half_x() {"""
test = """    /// BETA.DIST(x, 1, b, TRUE) below x = 0.1 is -expm1(b * log1p(-x)) on the x87 chain
    /// (R43). Excel 16.0 build 20527, the question service, 2026-11-09.
    #[cfg_attr(not(target_arch = "wasm32"), test)]
    #[cfg_attr(target_arch = "wasm32", wasm_bindgen_test::wasm_bindgen_test)]
    fn beta_dist_unit_first_shape_below_a_tenth_matches_excel() {
        let rows: [(f64, f64, u64); 6] = [
            (0.05, 3.0, 0x3fc24189374bc6a8),
            (0.001, 2.5, 0x3f6476f2d098e78f),
            (0.01, 10.0, 0x3fb87a6a94bfa203),
            (0.07, 0.5, 0x3fa23ebda17f7999),
            (0.003, 7.0, 0x3f954fb9cbf27f4d),
            (0.02, 1.5, 0x3f9e90ddc821a601),
        ];
        for (x, b, expected) in rows {
            let got = beta_dist_kernel(x, 1.0, b, true, 0.0, 1.0).unwrap();
            assert_eq!(got.to_bits(), expected, "BETA.DIST({x},1,{b},TRUE)");
        }
    }

"""
assert anchor in s
s = s.replace(anchor, test + anchor)
open(p1, "w").write(s)

p2 = "/home/ubuntu/projects/DnaCalc/DnaCalc/crates/calc-func/src/ported/discrete_dist.rs"
s = open(p2).read()
anchor2 = """    #[cfg_attr(not(target_arch = "wasm32"), test)]
    #[cfg_attr(target_arch = "wasm32", wasm_bindgen_test::wasm_bindgen_test)]
    fn argument_rules_match_excel() {"""
test2 = """    /// BINOM.DIST(n - 1, n, p, TRUE) with 1 - p below 0.1 is -expm1(n * log1p(p - 1)) with a
    /// binary64 product, not BETA.DIST(1 - p, 1, n): at n = 10, p = 0.99 the two Excel answers
    /// differ by five units. Excel 16.0 build 20527, the question service, 2026-11-09 (R43).
    #[cfg_attr(not(target_arch = "wasm32"), test)]
    #[cfg_attr(target_arch = "wasm32", wasm_bindgen_test::wasm_bindgen_test)]
    fn binomial_cdf_one_short_of_n_matches_excel() {
        let rows: [(f64, f64, f64, u64); 6] = [
            (19.0, 20.0, 0.95, 0x3fe4874887e5b9b4),
            (99.0, 100.0, 0.999, 0x3fb85f8ab6781056),
            (9.0, 10.0, 0.99, 0x3fb87a6a94bfa208),
            (49.0, 50.0, 0.97, 0x3fe9059bc30fc602),
            (4.0, 5.0, 0.93, 0x3fd379d77d69880c),
            (29.0, 30.0, 0.995, 0x3fc1deee48dd912d),
        ];
        for (k, n, p, expected) in rows {
            let got = binom_dist_kernel(k, n, p, true).unwrap();
            assert_eq!(got.to_bits(), expected, "BINOM.DIST({k},{n},{p},TRUE)");
        }
    }

"""
assert anchor2 in s
s = s.replace(anchor2, test2 + anchor2)
open(p2, "w").write(s)
print("ok")
