use std::{env, fs, path::PathBuf};
fn main() {
    let source = PathBuf::from("/home/koderbot/projects/DnaCalc/DnaCalc/crates/calc-func/src/ported/excel_numeric");
    let out = PathBuf::from(env::var_os("OUT_DIR").unwrap());
    let mut text = fs::read_to_string(source.join("mod.rs")).unwrap();
    for module in ["tables","x87","soft87","x87_chains","crt_pow","crt_libm","gammaln"] {
        let mut path = source.join(format!("{module}.rs"));
        println!("cargo:rerun-if-changed={}", path.display());
        if module == "x87_chains" {
            let chain = fs::read_to_string(&path).unwrap();
            path = out.join("chains.rs");
            fs::write(&path, chain.replace("fn exp_ext<", "pub(super) fn exp_ext<").replace("fn pow2<", "pub(super) fn pow2<")).unwrap();
        }
        text = text.replace(&format!("mod {module};"), &format!("#[path = {:?}] mod {module};", path));
    }
    text = text.replace("#[path = \"../../../tests/f1_support/x87_cpu.rs\"]", &format!("#[path = {:?}]", source.join("../../../tests/f1_support/x87_cpu.rs")));
    text = text.replace("include_str!(\"x87_excel_ground_truth.tsv\")", &format!("include_str!({:?})", source.join("x87_excel_ground_truth.tsv")));
    text.push_str("\npub(crate) use soft87::{Ext80 as XExt, ext_from_f64 as xload, ext_one as xone, ext_chs as xchs, ext_ln2 as xln2, ext_l2e as xl2e, CW_PC64_RN as XCW};\n");
    text.push_str("pub(crate) fn xexp(x: &Ext80) -> Ext80 { x87_chains::exp_ext::<Native>(x) }\n");
    text.push_str("pub(crate) fn xexpc(t: &Ext80, cw: u16) -> Ext80 { let p2 = |t: &Ext80| { let w = Native::rndint(t, cw); let fr = Native::sub(t, &w, cw); let m = Native::add(&Native::f2xm1(&fr, cw), &ext_one(), cw); Native::scale(&m, &w, cw) }; let y = Native::mul(t, &ext_l2e(), cw); if y.0[9] & 0x80 != 0 { Native::div(&ext_one(), &p2(&ext_chs(&y, cw)), cw) } else { p2(&y) } }\n");
    text.push_str("pub(crate) fn xaddc(a:&Ext80,b:&Ext80,cw:u16)->Ext80{Native::add(a,b,cw)}\npub(crate) fn xstorec(x: &Ext80, cw: u16) -> f64 { Native::store(x, cw) }\npub(crate) fn xlogc(x:&Ext80,cw:u16)->Ext80{Native::fyl2x(&ext_ln2(), x, cw)}\n");
    text.push_str("pub(crate) fn xstore(x: &Ext80) -> f64 { Native::store(x, CW) }
pub(crate) fn xstore_rz(x: &Ext80) -> f64 { Native::store(x, CW | 0x0C00) }
pub(crate) fn xmulc(a:&Ext80,b:&Ext80,cw:u16)->Ext80{Native::mul(a,b,cw)}
pub(crate) fn xdivc(a:&Ext80,b:&Ext80,cw:u16)->Ext80{Native::div(a,b,cw)}\n");
    text.push_str("pub(crate) fn xadd(a:&Ext80,b:&Ext80)->Ext80{Native::add(a,b,CW)}\npub(crate) fn xsub(a:&Ext80,b:&Ext80)->Ext80{Native::sub(a,b,CW)}\npub(crate) fn xmul(a:&Ext80,b:&Ext80)->Ext80{Native::mul(a,b,CW)}\npub(crate) fn xdiv(a:&Ext80,b:&Ext80)->Ext80{Native::div(a,b,CW)}\npub(crate) fn xlog(x:&Ext80)->Ext80{Native::fyl2x(&ext_ln2(), x, CW)}\n");
    println!("cargo:rerun-if-changed={}", source.join("mod.rs").display());
    fs::write(out.join("numeric.rs"), text).unwrap(); fs::create_dir_all("/home/koderbot/codex-runs/ld1/explore/gen").unwrap(); fs::write("/home/koderbot/codex-runs/ld1/explore/gen/numeric.rs", fs::read_to_string(out.join("numeric.rs")).unwrap()).unwrap();
}
