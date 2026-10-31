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
    assert!(
        std::process::Command::new("git")
            .current_dir(&checkout)
            .args([
                "diff",
                "--quiet",
                revision,
                "--",
                "crates/calc-func/src",
                "crates/calc-value/src",
                "tools/calc-verdict/src"
            ])
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
    // Only visibility and module locations change; all numeric operations are the
    // current cleared DnaCalc sources, compiled as research without editing them.
    text.push_str("\npub(crate) use soft87::{Ext80 as Extended, ext_from_f64 as research_load};\n");
    text.push_str("\npub(crate) fn research_exp<F: Fpu>(x: &Extended) -> Extended { x87_chains::exp_ext::<F>(x) }\n");
    text.push_str("\npub(crate) fn research_log<F: Fpu>(x: &Extended) -> Extended { F::fyl2x(&ext_ln2(), x, CW) }\n");
    text.push_str("\npub(crate) fn research_pow<F: Fpu>(b: f64, e: f64) -> Extended { x87_chains::pow2::<F>(&F::fyl2x(&load(e), &load(b), CW)) }\n");
    text.push_str("\npub(crate) fn research_pow2<F: Fpu>(t: &Extended) -> Extended { x87_chains::pow2::<F>(t) }\n");
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
}
