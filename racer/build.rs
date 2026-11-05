use std::{env, fs, path::PathBuf};

fn main() {
    // Old units keep their original guard; a new unit explicitly freezes its
    // own source revision before any research build uses it.
    println!("cargo:rerun-if-env-changed=DNA_RESEARCH_UNIT");
    println!("cargo:rustc-check-cfg=cfg(bessel_high_j_audit)");
    println!("cargo:rustc-check-cfg=cfg(bessel_high_groups_audit)");
    let unit_name = env::var("DNA_RESEARCH_UNIT").unwrap_or_else(|_| "pmt-r6-unit.json".into());
    let unit_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join(unit_name);
    let unit: serde_json::Value = serde_json::from_slice(&fs::read(&unit_path).unwrap()).unwrap();
    let revision = unit["application_revision"].as_str().unwrap();
    let checkout = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../DnaCalc");
    let mut guarded_paths = vec![
        "crates/calc-func/src",
        "crates/calc-value/src",
        "tools/calc-verdict/src",
    ];
    if unit["route_audit"].as_bool() == Some(true) {
        guarded_paths.extend([
            "apps/dna-calc-worker/src",
            "apps/dna-calc-model/src",
            "apps/dna-calc-protocol/src",
            "apps/dna-formula-ui/src",
            "apps/dna-strand/src",
            "crates/calc-host/src",
            "crates/calc-engine/src",
            "crates/calc-fml/src",
            "crates/calc-doc/src",
            "crates/calc-doc-model/src",
            "crates/calc-syntax/src",
        ]);
    }
    let untracked = std::process::Command::new("git")
        .current_dir(&checkout)
        .args(["ls-files", "--others", "--exclude-standard", "--"])
        .args(&guarded_paths)
        .output().unwrap();
    assert!(untracked.status.success());
    let untracked = std::str::from_utf8(&untracked.stdout).unwrap();
    if let Some(blobs) = unit["application_working_blobs"].as_object() {
        // A pending product change can be studied before its commit only if
        // the unit binds every changed guarded file to its exact Git blob.
        // Historical revision-only units retain the original stricter mode.
        let diff = std::process::Command::new("git")
            .current_dir(&checkout)
            .args(["diff", "--name-only", revision, "--"])
            .args(&guarded_paths)
            .output().unwrap();
        assert!(diff.status.success());
        let mut changed: std::collections::BTreeSet<_> = std::str::from_utf8(&diff.stdout).unwrap()
            .lines().collect();
        changed.extend(untracked.lines());
        let declared: std::collections::BTreeSet<_> = blobs.keys().map(String::as_str).collect();
        assert_eq!(changed, declared, "unbound changed application sources");
        for (path, expected) in blobs {
            assert!(path.starts_with("crates/") || path.starts_with("apps/") || path.starts_with("tools/"));
            let blob = std::process::Command::new("git").current_dir(&checkout)
                .args(["hash-object", "--", path]).output().unwrap();
            assert!(blob.status.success());
            assert_eq!(std::str::from_utf8(&blob.stdout).unwrap().trim(), expected.as_str().unwrap(),
                "application source differs from frozen working blob: {path}");
            println!("cargo:rerun-if-changed={}", checkout.join(path).display());
        }
    } else {
        assert!(untracked.is_empty(), "unbound untracked application sources");
        assert!(std::process::Command::new("git")
            .current_dir(&checkout)
            .args(["diff", "--quiet", revision, "--"])
            .args(guarded_paths)
            .status()
            .unwrap()
            .success(),
        "DnaCalc sources differ from the research unit's frozen revision"
        );
    }
    println!("cargo:rerun-if-changed={}", unit_path.display());
    let source = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../DnaCalc/crates/calc-func/src/ported/excel_numeric")
        .canonicalize()
        .expect("DnaCalc sibling checkout");
    let mut text = fs::read_to_string(source.join("mod.rs")).unwrap();
    let output = PathBuf::from(env::var_os("OUT_DIR").unwrap());
    // Freeze the admitted unit inside the executable. Packet metadata must
    // describe its compiled application sources rather than the current HEAD.
    let mut binding=unit.clone();
    binding["compiled_features"]=serde_json::json!({
        "soft_x87":env::var_os("CARGO_FEATURE_SOFT_X87").is_some()
    });
    fs::write(output.join("research-unit.json"),serde_json::to_vec_pretty(&binding).unwrap()).unwrap();
    for module in [
        "tables",
        "x87",
        "soft87",
        "x87_chains",
        "crt_pow",
        "crt_libm",
        "gammaln",
    ] {
        let mut path = source.join(format!("{module}.rs"));
        if module == "x87_chains" {
            let chain = fs::read_to_string(&path).unwrap();
            println!("cargo:rerun-if-changed={}", path.display());
            path = output.join("chains-research.rs");
            fs::write(
                &path,
                chain
                    .replace("fn exp_ext<", "pub(super) fn exp_ext<")
                    .replace("fn pow2<", "pub(super) fn pow2<"),
            )
            .unwrap();
        }
        text = text.replace(
            &format!("mod {module};"),
            &format!("#[path = {:?}] mod {module};", path),
        );
        println!("cargo:rerun-if-changed={}", path.display());
    }
    // Keep the existing numeric module's tests usable from this generated
    // location too. These replacements change paths, never test assertions.
    let cpu_tests = source
        .join("../../../tests/f1_support/x87_cpu.rs")
        .canonicalize()
        .unwrap();
    text = text.replace(
        "#[path = \"../../../tests/f1_support/x87_cpu.rs\"]",
        &format!("#[path = {:?}]", cpu_tests),
    );
    let ground = source.join("x87_excel_ground_truth.tsv");
    text = text.replace(
        "include_str!(\"x87_excel_ground_truth.tsv\")",
        &format!("include_str!({:?})", ground),
    );
    println!("cargo:rerun-if-changed={}", cpu_tests.display());
    println!("cargo:rerun-if-changed={}", ground.display());
    // Only visibility and module locations change; all numeric operations are the
    // current cleared DnaCalc sources, compiled as research without editing them.
    text.push_str("\npub(crate) use soft87::{Ext80 as Extended, ext_from_f64 as research_load};\n");
    text.push_str("\npub(crate) fn research_exp<F: Fpu>(x: &Extended) -> Extended { x87_chains::exp_ext::<F>(x) }\n");
    text.push_str("\npub(crate) fn research_log<F: Fpu>(x: &Extended) -> Extended { F::fyl2x(&ext_ln2(), x, CW) }\n");
    text.push_str("\npub(crate) fn research_pow<F: Fpu>(b: f64, e: f64) -> Extended { x87_chains::pow2::<F>(&F::fyl2x(&load(e), &load(b), CW)) }\n");
    text.push_str("\npub(crate) fn research_pow2<F: Fpu>(t: &Extended) -> Extended { x87_chains::pow2::<F>(t) }\n");
    text.push_str("\npub(crate) fn research_ln2() -> Extended { ext_ln2() }\n");
    text.push_str("\npub(crate) fn research_l2e() -> Extended { ext_l2e() }\n");
    println!("cargo:rerun-if-changed={}", source.join("mod.rs").display());
    fs::write(output.join("numeric.rs"), text).unwrap();
    fs::write(
        output.join("bridge.rs"),
        format!(
            "#[path = {:?}] pub mod excel_numeric;",
            output.join("numeric.rs")
        ),
    )
    .unwrap();
    // Audit the current cleared bodies without changing their calculations.
    // Old units retain their original bridge; no historical executable is inferred.
    if unit["route_audit"].as_bool() == Some(true) {
        use std::io::Write;
        let mut bridge = fs::OpenOptions::new()
            .append(true)
            .open(output.join("bridge.rs"))
            .unwrap();
        for module in [
            "power_fn",
            "financial_time_value_family",
            "real_policy",
            "exp_fn",
            "complex_decimal",
            "complex_family",
            "bessel",
        ] {
            let path = source.parent().unwrap().join(format!("{module}.rs"));
            println!("cargo:rerun-if-changed={}", path.display());
            writeln!(bridge, "#[path = {:?}] pub mod {module};", path).unwrap();
        }
        if unit["bessel_rational_audit"].as_bool() == Some(true) {
            let path = source.parent().unwrap().join("bessel_rational.rs");
            println!("cargo:rerun-if-changed={}", path.display());
            writeln!(bridge, "#[path = {:?}] pub mod bessel_rational;", path).unwrap();
        }
        if unit["bessel_high_j_audit"].as_bool() == Some(true) {
            let original=source.parent().unwrap().join("bessel.rs");
            let mut text=fs::read_to_string(&original).unwrap();
            for (old,new) in [
                ("libm::j0(x.0)","crate::bessel_rational::inferred_high_j_seed(0,x.0)"),
                ("libm::j1(x.0)","crate::bessel_rational::inferred_high_j_seed(1,x.0)"),
            ] {
                assert_eq!(text.matches(old).count(),1,"changed high-J seed call");
                text=text.replace(old,new);
            }
            let path=output.join("bessel-high-j-research.rs");fs::write(&path,text).unwrap();
            writeln!(bridge,"#[path = {:?}] pub mod bessel_high_j;",path).unwrap();
            println!("cargo:rerun-if-changed={}",original.display());
            println!("cargo:rustc-cfg=bessel_high_j_audit");
        }
        if unit["bessel_high_groups_audit"].as_bool() == Some(true) {
            let original=source.parent().unwrap().join("bessel_rational.rs");
            let mut text=fs::read_to_string(&original).unwrap();
            for (trig,count) in [("sin",1),("cos",2)] {
                let old=format!("((excel_{trig}(phase) * correction) * 8.0) / x");
                let new=format!("((8.0 / x) * excel_{trig}(phase)) * correction");
                assert_eq!(text.matches(&old).count(),count,"changed high-region tail");
                text=text.replace(&old,&new);
            }
            let path=output.join("bessel-grouped-rational.rs");fs::write(&path,text).unwrap();
            writeln!(bridge,"#[path = {:?}] pub mod bessel_grouped_rational;",path).unwrap();
            println!("cargo:rerun-if-changed={}",original.display());
            let original=source.parent().unwrap().join("bessel.rs");
            let text=fs::read_to_string(&original).unwrap();
            assert!(text.contains("super::bessel_rational::"));
            let text=text.replace("super::bessel_rational::","super::bessel_grouped_rational::");
            let path=output.join("bessel-grouped.rs");fs::write(&path,text).unwrap();
            writeln!(bridge,"#[path = {:?}] pub mod bessel_grouped;",path).unwrap();
            println!("cargo:rerun-if-changed={}",original.display());
            println!("cargo:rustc-cfg=bessel_high_groups_audit");
        }
        if unit["bond_solver_audit"].as_bool() == Some(true) {
            let day_count = source.parent().unwrap().join("day_count_common.rs");
            writeln!(bridge, "#[path = {:?}] pub mod day_count_common;", day_count).unwrap();
            println!("cargo:rerun-if-changed={}", day_count.display());
            let original = source.parent().unwrap().join("bond_core_family.rs");
            let mut bond = fs::read_to_string(&original).unwrap();
            // A separately named diagnostic price copy changes one date-factor
            // producer only. The original registered body remains untouched.
            let begin=bond.find("fn pcomp_disc(").unwrap();
            let end=bond[begin..].find("\nfn solve(").unwrap()+begin;
            let alternative=bond[begin..end].replace("fn pcomp_disc(","fn research_pcomp_days_num(")
                .replace("let a = dd(p.prev, c.settlement, c.basis)?;",
                         "let a = days_between_num(p.prev, c.settlement, c.basis)?;");
            assert!(alternative.contains("let a = days_between_num(p.prev, c.settlement, c.basis)?;"));
            bond.push_str(&alternative);
            // The current price/date bodies are unchanged. This appended
            // wrapper exposes their existing private schedule inside the same
            // module, solely for the research solver comparison.
            bond.push_str(r#"
pub(crate) struct ResearchPrice { context: Ctx, coupons: Period, coupon_rate: f64, redemption: f64 }
pub(crate) fn research_mat_factors(a:&[f64],numerator_days:bool)->Result<(f64,f64,f64),BondCoreEvalError> {
    let settlement=dser(a[0])?;let maturity=dser(a[1])?;let issue=dser(a[2])?;
    let basis=basis(a.get(5).copied().unwrap_or(0.0))?;
    let b=days_in_year_for_mat(issue,settlement,basis)?;
    let days=if numerator_days {days_between_num} else {dd};
    Ok((b,days(issue,maturity,basis)?,days(issue,settlement,basis)?))
}
impl ResearchPrice {
    pub(crate) fn new(a: &[f64]) -> Result<Self, BondCoreEvalError> {
        let context = ctx(a[0], a[1], a[5], a.get(6).copied())?;
        Ok(Self { context, coupons: period(context)?, coupon_rate: rate(a[2])?, redemption: pos(a[4])? })
    }
    pub(crate) fn count(&self) -> i64 { self.coupons.n }
    pub(crate) fn value(&self, y: f64, staged: bool) -> Result<f64, BondCoreEvalError> {
        Ok(pcomp_disc(self.coupon_rate, y, self.redemption, self.context, self.coupons, staged)?.0)
    }
    pub(crate) fn is_par(&self, target: f64) -> bool {
        self.coupons.n > 1 && self.context.settlement == self.coupons.prev
            && target == 100.0 && self.redemption == 100.0
            && (100.0*self.coupon_rate/self.context.frequency as f64).is_finite()
    }
    pub(crate) fn original_root(&self, target: f64, staged: bool, numerator_days: bool) -> Result<f64, BondCoreEvalError> {
        solve(target,self.context.frequency,|y| {
            let price=if numerator_days {research_pcomp_days_num} else {pcomp_disc};
            Ok(price(self.coupon_rate,y,self.redemption,self.context,self.coupons,staged)?.0)
        })
    }
}
"#);
            let path = output.join("bond-research.rs");
            fs::write(&path, bond).unwrap();
            writeln!(bridge, "#[path = {:?}] pub mod bond_core_family;", path).unwrap();
            println!("cargo:rerun-if-changed={}", original.display());
            println!("cargo:rustc-cfg=bond_solver_audit");
        }
        if unit["closed_security_audit"].as_bool() == Some(true) {
            let original=source.parent().unwrap().join("discount_bill_yearfrac_family.rs");
            let mut body=fs::read_to_string(&original).unwrap();
            body.push_str(r#"
pub(crate) fn research_security_days(a:&[f64],bill:bool)->Result<(f64,f64),DiscountBillYearfracEvalError> {
    if bill {Ok((tbill_days(a[0],a[1])?,360.0))}
    else {security_days_and_year(a[0],a[1],a.get(4).copied())}
}
"#);
            let path=output.join("security-research.rs");fs::write(&path,body).unwrap();
            writeln!(bridge,"#[path = {:?}] pub mod discount_bill_yearfrac_family;",path).unwrap();
            println!("cargo:rerun-if-changed={}",original.display());
            println!("cargo:rustc-cfg=closed_security_audit");
        }
        if unit["bessel_k_audit"].as_bool() == Some(true) {
            let original=source.parent().unwrap().join("bessel.rs");
            let original_body=fs::read_to_string(&original).unwrap();
            // Numerical research templates only. The registered module and
            // its full application tests remain the original source above.
            let (body,_)=original_body.split_once("\n#[cfg(test)]\nmod tests {").unwrap();
            let start=body.find("macro_rules! flushing_op").unwrap();
            let end=body.find("\nimpl Neg for D").unwrap();
            let replacement=r#"macro_rules! flushing_op {
    ($tr:ident, $f:ident, $helper:ident) => {
        impl $tr for D {
            type Output = D;
            fn $f(self, rhs: D) -> D {
                D::new(super::excel_numeric::$helper(D::new(self.0).0, D::new(rhs.0).0))
            }
        }
    };
}
flushing_op!(Add, add, excel_x87_add);
flushing_op!(Sub, sub, excel_x87_sub);
flushing_op!(Mul, mul, excel_x87_mul);
flushing_op!(Div, div, excel_x87_div);
"#;
            let body=format!("{}{}{}",&body[..start],replacement,&body[end..]);
            for (name,staged_sqrt) in [("bessel_k_pc64",false),("bessel_k_pc64_sqrt",true)] {
                let mut text=body.clone();
                if staged_sqrt {
                    assert_eq!(text.matches("D::new(self.0.sqrt())").count(),1);
                    text=text.replace("D::new(self.0.sqrt())","D::new(super::excel_numeric::excel_sqrt_staged(self.0))");
                }
                let path=output.join(format!("{name}.rs"));fs::write(&path,text).unwrap();
                writeln!(bridge,"#[path = {:?}] pub mod {name};",path).unwrap();
            }
            println!("cargo:rerun-if-changed={}",original.display());
            println!("cargo:rustc-cfg=bessel_k_audit");
        }
        println!("cargo:rustc-cfg=route_audit");
    }
    println!("cargo:rustc-check-cfg=cfg(route_audit)");
    println!("cargo:rustc-check-cfg=cfg(bond_solver_audit)");
    println!("cargo:rustc-check-cfg=cfg(closed_security_audit)");
    println!("cargo:rustc-check-cfg=cfg(bessel_k_audit)");
    if unit["discrete_audit"].as_bool() == Some(true) {
        use std::io::Write;
        let original = source.parent().unwrap().join("discrete_dist.rs");
        let text = fs::read_to_string(&original).unwrap();
        // Expose current cleared primitives; the copied body and constants
        // are unchanged. Hypotheses live separately in discrete.rs.
        let text = text.replace("fn stirling_error(", "pub(crate) fn stirling_error(")
            .replace("fn deviance(", "pub(crate) fn deviance(")
            .replace("fn binom_mass(", "pub(crate) fn binom_mass(");
        let path = output.join("discrete-research.rs");
        fs::write(&path, text).unwrap();
        let mut bridge = fs::OpenOptions::new().append(true).open(output.join("bridge.rs")).unwrap();
        writeln!(bridge, "#[path = {:?}] pub mod discrete_dist;", path).unwrap();
        let special = source.parent().unwrap().join("special_math_common.rs");
        // Visibility only: the original kept TOMS708 arithmetic is unchanged.
        let special_text=fs::read_to_string(&special).unwrap();
        assert_eq!(special_text.matches("fn bratio_brcomp(").count(),1);
        let generated_special=output.join("special-research.rs");
        let mut special_text=special_text.replace("fn bratio_brcomp(","pub(crate) fn bratio_brcomp(");
        if unit["gamma_exponent_audit"].as_bool()==Some(true) {
            assert_eq!(special_text.matches("a * crt_log(x) - x").count(),3,"changed gamma exponent sites");
            special_text=special_text.replace("a * crt_log(x) - x","crate::discrete::gamma_exponent(a,x)");
            let finish="fn gratio_temme_finish(l: f64, c: f64, w: f64, t: f64, a: f64) -> (f64, f64) {";
            assert_eq!(special_text.matches(finish).count(),1,"changed Temme finish site");
            special_text=special_text.replace(finish,&format!("{finish}\n    if let Some(result)=crate::discrete::gamma_temme_finish(l,c,w,t,a,GRATIO_RT2PIN) {{ return result; }}"));
            println!("cargo:rustc-cfg=gamma_exponent_audit");
        }
        fs::write(&generated_special,special_text).unwrap();
        writeln!(bridge, "#[path = {:?}] pub mod special_math_common;", generated_special).unwrap();
        println!("cargo:rerun-if-changed={}", original.display());
        println!("cargo:rerun-if-changed={}", special.display());
        println!("cargo:rustc-cfg=discrete_audit");
    }
    println!("cargo:rustc-check-cfg=cfg(discrete_audit)");
    println!("cargo:rustc-check-cfg=cfg(gamma_exponent_audit)");
    if unit["regression_audit"].as_bool() == Some(true) {
        use std::io::Write;
        let original=source.parent().unwrap().join("regression_forecast_family.rs");
        let mut text=fs::read_to_string(&original).unwrap();
        // Observe only final standard-error operands. Every original arithmetic
        // expression is retained, and the collector checks the selected result
        // against the registered invocation before admitting any captured leaf.
        for (old,new) in [
            ("(ss_resid * covariance(position, position) / df).sqrt()",
             "crate::route_audit::regression_leaf(ss_resid, covariance(position, position), df, 2, fit.model.coefficients.len()-*index, 0, (ss_resid * covariance(position, position) / df).sqrt())"),
            ("(ss_resid / df * variance).sqrt()",
             "crate::route_audit::regression_leaf(ss_resid, variance, df, 2, first.len(), 1, (ss_resid / df * variance).sqrt())"),
            ("let se_b = if exact",
             "let variance = crate::route_audit::regression_factors(&w, &inverse, n, variance); let se_b = if exact"),
        ] {
            assert_eq!(text.matches(old).count(),1,"changed regression observation site");
            text=text.replace(old,new);
        }
        if unit["regression_inverse_audit"].as_bool() == Some(true) {
            let site = "    // (R'R)^-1 = R^-1 R^-T.";
            assert_eq!(text.matches(site).count(), 1, "changed inverse-R observation site");
            text = text.replace(site, &format!(
                "    crate::route_audit::regression_triangular(&fit.r, nz, &inverse, &fit.order, fit.model.coefficients.len());\n{site}"));
        }
        if unit["regression_outer_audit"].as_bool() == Some(true) {
            // Six coherent operation groups surrounding the unchanged reflector.
            // No rank threshold, coefficient, pivot tie or data-dependent selector changes.
            for (old, new) in [
                ("values.iter().fold(0.0, |sum, v| sum + v) / n as f64", "crate::route_audit::regression_outer_op(0, 3, values.iter().fold(0.0, |sum, v| crate::route_audit::regression_outer_op(0, 0, sum, *v)), n as f64)"),
                ("column.iter().map(|v| v - means[index]).collect()", "column.iter().map(|v| crate::route_audit::regression_outer_op(1, 1, *v, means[index])).collect()"),
                ("y.iter().map(|v| v - y_mean).collect()", "y.iter().map(|v| crate::route_audit::regression_outer_op(1, 1, *v, y_mean)).collect()"),
                ("values[from..].iter().fold(0.0, |sum, v| sum + v * v)", "values[from..].iter().fold(0.0, |sum, v| crate::route_audit::regression_outer_op(2, 0, sum, crate::route_audit::regression_outer_op(2, 2, *v, *v)))"),
                ("let norm = squares.sqrt();", "let norm = crate::route_audit::regression_outer_op(3, 4, squares, 0.0);"),
                ("((j + 1)..rank).fold(0.0, |sum, k| sum + matrix[k].1[j] * solution[k])", "((j + 1)..rank).fold(0.0, |sum, k| crate::route_audit::regression_outer_op(4, 0, sum, crate::route_audit::regression_outer_op(4, 2, matrix[k].1[j], solution[k])))"),
                ("(target[j] - d) / diagonal[j]", "crate::route_audit::regression_outer_op(4, 3, crate::route_audit::regression_outer_op(4, 1, target[j], d), diagonal[j])"),
                (".fold(y_mean, |a, (b, m)| a - b * m)", ".fold(y_mean, |a, (b, m)| crate::route_audit::regression_outer_op(5, 1, a, crate::route_audit::regression_outer_op(5, 2, *b, *m)))"),
            ] {
                assert_eq!(text.matches(old).count(), 1, "changed outer QR audit site: {old}");
                text = text.replace(old, new);
            }
            for (old,new) in [
                ("values.iter().fold(0.0, |sum, v| crate::route_audit::regression_outer_op(0, 0, sum, *v))", "crate::route_audit::regression_outer_fold(0, 0.0, values.iter().map(|v| (*v,1.0)))"),
                ("values[from..].iter().fold(0.0, |sum, v| crate::route_audit::regression_outer_op(2, 0, sum, crate::route_audit::regression_outer_op(2, 2, *v, *v)))", "crate::route_audit::regression_outer_fold(2, 0.0, values[from..].iter().map(|v| (*v,*v)))"),
                ("((j + 1)..rank).fold(0.0, |sum, k| crate::route_audit::regression_outer_op(4, 0, sum, crate::route_audit::regression_outer_op(4, 2, matrix[k].1[j], solution[k])))", "crate::route_audit::regression_outer_fold(4, 0.0, ((j+1)..rank).map(|k| (matrix[k].1[j],solution[k])))"),
                ("coefficients\n            .iter()\n            .zip(means.iter())\n            .fold(y_mean, |a, (b, m)| crate::route_audit::regression_outer_op(5, 1, a, crate::route_audit::regression_outer_op(5, 2, *b, *m)))", "crate::route_audit::regression_outer_fold(5, y_mean, coefficients.iter().zip(means.iter()).map(|(b,m)| (*b,*m)))"),
            ] {
                assert_eq!(text.matches(old).count(),1,"changed QR fold site: {old}");
                text=text.replace(old,new);
            }
            println!("cargo:rustc-cfg=regression_outer_audit");
        }
        if unit["regression_householder_audit"].as_bool() == Some(true) {
            let first="        let norm = squares.sqrt();";
            let last="        diagonal[j] = alpha;";
            assert_eq!(text.matches(first).count(),1,"changed reflector start");
            assert_eq!(text.matches(last).count(),1,"changed reflector end");
            let start=text.find(first).unwrap();let end=text.find(last).unwrap();assert!(start<end);
            text.replace_range(start..end,"        let alpha = crate::route_audit::regression_householder_step(&mut matrix, &mut target, j, x0, rest, squares);\n");
            println!("cargo:rustc-cfg=regression_householder_audit");
        }
        let path=output.join("regression-research.rs");fs::write(&path,text).unwrap();
        let mut bridge=fs::OpenOptions::new().append(true).open(output.join("bridge.rs")).unwrap();
        writeln!(bridge,"#[path = {:?}] pub mod regression_forecast_family;",path).unwrap();
        println!("cargo:rerun-if-changed={}",original.display());
        println!("cargo:rustc-cfg=regression_audit");
    }
    println!("cargo:rustc-check-cfg=cfg(regression_audit)");
    println!("cargo:rustc-check-cfg=cfg(regression_householder_audit)");
    println!("cargo:rustc-check-cfg=cfg(regression_outer_audit)");
}
