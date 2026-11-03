use std::{env, fs, path::PathBuf};

fn main() {
    // Old units keep their original guard; a new unit explicitly freezes its
    // own source revision before any research build uses it.
    println!("cargo:rerun-if-env-changed=DNA_RESEARCH_UNIT");
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
    assert!(
        std::process::Command::new("git")
            .current_dir(&checkout)
            .args(["diff", "--quiet", revision, "--"])
            .args(guarded_paths)
            .status()
            .unwrap()
            .success(),
        "DnaCalc sources differ from the research unit's frozen revision"
    );
    println!("cargo:rerun-if-changed={}", unit_path.display());
    let source = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../DnaCalc/crates/calc-func/src/ported/excel_numeric")
        .canonicalize()
        .expect("DnaCalc sibling checkout");
    let mut text = fs::read_to_string(source.join("mod.rs")).unwrap();
    let output = PathBuf::from(env::var_os("OUT_DIR").unwrap());
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
        if unit["bond_solver_audit"].as_bool() == Some(true) {
            let day_count = source.parent().unwrap().join("day_count_common.rs");
            writeln!(bridge, "#[path = {:?}] pub mod day_count_common;", day_count).unwrap();
            println!("cargo:rerun-if-changed={}", day_count.display());
            let original = source.parent().unwrap().join("bond_core_family.rs");
            let mut bond = fs::read_to_string(&original).unwrap();
            // The current price/date bodies are unchanged. This appended
            // wrapper exposes their existing private schedule inside the same
            // module, solely for the research solver comparison.
            bond.push_str(r#"
pub(crate) struct ResearchPrice { context: Ctx, coupons: Period, coupon_rate: f64, redemption: f64 }
impl ResearchPrice {
    pub(crate) fn new(a: &[f64]) -> Result<Self, BondCoreEvalError> {
        let context = ctx(a[0], a[1], a[5], a.get(6).copied())?;
        Ok(Self { context, coupons: period(context)?, coupon_rate: rate(a[2])?, redemption: pos(a[4])? })
    }
    pub(crate) fn count(&self) -> i64 { self.coupons.n }
    pub(crate) fn value(&self, y: f64, staged: bool) -> Result<f64, BondCoreEvalError> {
        Ok(pcomp_disc(self.coupon_rate, y, self.redemption, self.context, self.coupons, staged)?.0)
    }
}
"#);
            let path = output.join("bond-research.rs");
            fs::write(&path, bond).unwrap();
            writeln!(bridge, "#[path = {:?}] pub mod bond_core_family;", path).unwrap();
            println!("cargo:rerun-if-changed={}", original.display());
            println!("cargo:rustc-cfg=bond_solver_audit");
        }
        println!("cargo:rustc-cfg=route_audit");
    }
    println!("cargo:rustc-check-cfg=cfg(route_audit)");
    println!("cargo:rustc-check-cfg=cfg(bond_solver_audit)");
}
