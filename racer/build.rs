use std::{env, fs, path::PathBuf};

fn main() {
    let unit_path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../pmt-r6-unit.json");
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
    for module in [
        "tables",
        "x87",
        "soft87",
        "x87_chains",
        "crt_pow",
        "crt_libm",
        "gammaln",
    ] {
        let path = source.join(format!("{module}.rs"));
        text = text.replace(
            &format!("mod {module};"),
            &format!("#[path = {:?}] mod {module};", path),
        );
        println!("cargo:rerun-if-changed={}", path.display());
    }
    // Only visibility and module locations change; all numeric operations are the
    // current cleared DnaCalc sources, compiled as research without editing them.
    text.push_str("\npub(crate) use soft87::{Ext80 as Extended, ext_from_f64 as research_load};\n");
    println!("cargo:rerun-if-changed={}", source.join("mod.rs").display());
    let output = PathBuf::from(env::var_os("OUT_DIR").unwrap());
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
