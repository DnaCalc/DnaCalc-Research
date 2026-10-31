//! Earlier PMT structure: continuous log1p, coupled growth/expm1 components.
//! Implementer; Intel SDM FYL2XP1, DLMF4.6.6, annuity algebra, cleared Fpu.
//! FYL2XP1 is a measured hardware leaf even in the soft-arithmetic replay.
use crate::{
    Row,
    boundary::{self, store},
    dsl::*,
    eval::*,
    json_write,
    ported::excel_numeric::*,
};
use serde_json::json;
use std::{
    fs,
    io::{BufWriter, Write},
    path::Path,
};

fn op<F: Fpu>(kind: Op, a: Val, b: Val) -> Val {
    binop::<F>(kind, EvalModel::ResidentPc64, a, b)
}
fn add<F: Fpu>(a: Val, b: Val) -> Val {
    op::<F>(Op::Add(0, 1), a, b)
}
fn mul<F: Fpu>(a: Val, b: Val) -> Val {
    op::<F>(Op::Mul(0, 1), a, b)
}
fn div<F: Fpu>(a: Val, b: Val) -> Val {
    op::<F>(Op::Div(0, 1), a, b)
}

fn hardware_log1p(rate: f64) -> Val {
    assert!(rate.abs() <= 0.25 && rate.is_finite());
    let ln2 = research_ln2();
    let mut result = Extended([0; 10]);
    let cw = CW;
    let mut saved = 0u16;
    // SAFETY: all pointers are live locals; two loads, one FYL2XP1 pop and
    // one store/pop balance the x87 stack; restore the caller's control word.
    unsafe {
        core::arch::asm!(
            "fnstcw word ptr [{saved}]","fldcw word ptr [{cw}]",
            "fld tbyte ptr [{ln2}]","fld qword ptr [{rate}]","fyl2xp1",
            "fstp tbyte ptr [{result}]","fldcw word ptr [{saved}]",
            saved=in(reg)&mut saved,cw=in(reg)&cw,ln2=in(reg)&ln2,rate=in(reg)&rate,result=in(reg)&mut result,
            options(nostack,preserves_flags),
        );
    }
    Val::Ext(result)
}
fn series_log1p<F: Fpu>(rate: f64) -> Val {
    // DLMF4.6.6, a=1: log(1+r)=2*sum(z^(2k+1)/(2k+1)), z=r/(2+r).
    // PC64 floating summation stops when the stored register value is unchanged.
    assert!(rate.abs() <= 0.25 && rate.is_finite());
    let z = div::<F>(Val::f64(rate), add::<F>(Val::f64(2.0), Val::f64(rate)));
    let square = mul::<F>(z, z);
    let mut term = z;
    let mut sum = z;
    for k in (3..=1023).step_by(2) {
        term = mul::<F>(term, square);
        let next = add::<F>(sum, div::<F>(term, Val::f64(k as f64)));
        if next == sum {
            return mul::<F>(Val::f64(2.0), sum);
        }
        sum = next;
    }
    panic!("log1p floating series did not reach its declared unchanged sum")
}
fn logs<F: Fpu>(rate: f64, frozen: u64) -> [Val; 5] {
    let hw = hardware_log1p(rate);
    let series = series_log1p::<F>(rate);
    [
        Val::F64(frozen),
        hw,
        store::<F>(hw, true),
        series,
        store::<F>(series, true),
    ]
}
fn nper(row: &Row) -> f64 {
    let arg = calc_verdict::decode_argument(&row.raw["args"][1], "nper bits").unwrap();
    match arg.core() {
        crate::value::CoreValue::Number(x) => *x,
        _ => panic!("numeric period scope"),
    }
}

pub fn controls(input: &str, out: &str) {
    let rows = crate::lines::<Row>(input);
    let hs = boundary::helpers();
    let mut checks = 0usize;
    let mut differing = 0usize;
    let mut f = BufWriter::new(fs::File::create(out).unwrap());
    for row in rows.iter().filter(|r| r.args()[0].abs() <= 0.25).take(64) {
        let r = row.args()[0];
        let n = nper(row);
        let hard = logs::<Hard>(r, row.log1p);
        let soft = logs::<Soft>(r, row.log1p);
        assert_eq!(hard[1].exact_encoding(), soft[1].exact_encoding());
        assert_eq!(hard[3].exact_encoding(), soft[3].exact_encoding());
        for h in &hs {
            let (_, _, em) = boundary::helper_components::<Hard>(h, n, Val::F64(row.log1p));
            assert_eq!(
                em.exact_encoding(),
                boundary::helper::<Hard>(h, r, n, f64::from_bits(row.log1p)).exact_encoding()
            );
            checks += 1;
        }
        differing += (hard[1].to_f64::<Hard>().to_bits() != row.log1p) as usize;
        serde_json::to_writer(&mut f,&json!({"witness":row.raw,"native_log_values":hard.map(Val::exact_encoding),"soft_arithmetic_log_values":soft.map(Val::exact_encoding)})).unwrap();
        writeln!(f).unwrap();
    }
    println!(
        "{}",
        json!({"helper_wrapper_checks":checks,"sample_stored_hardware_vs_frozen_differences":differing,"FYL2XP1":"native hardware leaf on both arithmetic backends; no software instruction emulation claim"})
    );
}

fn scoped_rows(input: &str, root: &Path) -> (Vec<Row>, std::collections::BTreeMap<String, usize>) {
    let original = crate::lines::<Row>(input);
    let mut rows = vec![];
    let mut counts = std::collections::BTreeMap::<String, usize>::new();
    let mut excluded = BufWriter::new(fs::File::create(root.join("excluded.jsonl")).unwrap());
    for row in original {
        let a = row.args();
        let tau = nper(&row) * f64::from_bits(row.log1p);
        let reason = if a[2] != 0.0 {
            "FV_nonzero_control"
        } else if a[0] <= 0.0 || a[0] > 0.25 {
            "outside_FYL2XP1_safe_domain"
        } else if !(tau > 0.0 && tau < 1.0) {
            "outside_small_exponent_arm"
        } else {
            "scoped"
        };
        *counts.entry(reason.into()).or_default() += 1;
        if reason == "scoped" {
            rows.push(row)
        } else {
            serde_json::to_writer(&mut excluded, &json!({"reason":reason,"row":row})).unwrap();
            writeln!(excluded).unwrap();
        }
    }
    (rows, counts)
}
fn direct_components<F: Fpu>(id: usize, n: f64, log: Val) -> (Val, Val, Val) {
    let mut code = id / 2;
    let mut next = || { let m = MODELS[code % 3]; code /= 3; m };
    let t = binop::<F>(Op::Mul(0, 1), next(), Val::f64(n), log);
    let x = binop::<F>(Op::Mul(0, 1), next(), t, Val::Ext(research_l2e()));
    let k = Val::Ext(F::rndint(&x.to_ext(), CW));
    let fraction = binop::<F>(Op::Add(0, 1), next(), x, boundary::neg(k));
    // Intel's F2XM1 argument interval is [-1,1]; nearest integer reduction
    // confines this fractional argument to [-1/2,1/2]. No answer-derived cut.
    assert!(fraction.to_f64::<F>().abs() <= 0.5);
    let f = store::<F>(Val::Ext(F::f2xm1(&fraction.to_ext(), CW)), id % 2 != 0);
    let scaled = Val::Ext(F::scale(&f.to_ext(), &k.to_ext(), CW));
    let two_k = Val::Ext(F::scale(&Val::f64(1.0).to_ext(), &k.to_ext(), CW));
    let integer_part = binop::<F>(Op::Add(0, 1), EvalModel::ResidentPc64, two_k, Val::f64(-1.0));
    let em = binop::<F>(Op::Add(0, 1), next(), scaled, integer_part);
    assert_eq!(code, 0);
    (f, k, em)
}

pub fn direct_controls(input: &str, out: &str) {
    let rows = crate::lines::<Row>(input);
    let mut f = BufWriter::new(fs::File::create(out).unwrap());
    let mut count = 0;
    for row in rows.iter().filter(|r| r.args()[0]>0.0 && r.args()[0]<=0.25 && nper(r)*f64::from_bits(r.log1p)>0.0 && nper(r)*f64::from_bits(r.log1p)<1.0).take(64) {
        let lh=logs::<Hard>(row.args()[0],row.log1p);
        let ls=logs::<Soft>(row.args()[0],row.log1p);
        for lp in 0..5 {
            for helper in 0..162 {
                let (hf,hk,he)=direct_components::<Hard>(helper,nper(row),lh[lp]);
                let (sf,sk,se)=direct_components::<Soft>(helper,nper(row),ls[lp]);
                serde_json::to_writer(&mut f,&json!({"witness":row.raw,"log_provider":lp,"helper":helper,
                    "native":[hf.exact_encoding(),hk.exact_encoding(),he.exact_encoding()],
                    "hybrid":[sf.exact_encoding(),sk.exact_encoding(),se.exact_encoding()]})).unwrap();
                writeln!(f).unwrap();count+=1;
            }
        }
    }
    println!("{}",json!({"paired_controls":count,"policy":"Raw direct F2XM1 components on original inputs; differences are retained, never asserted away. FYL2XP1 uses native hardware in both runs."}));
}

pub fn run<F: Fpu>(input: &str, out: &str, direct: bool) {
    fs::create_dir_all(out).unwrap();
    let root = Path::new(out);
    let hs = if direct { vec![] } else { boundary::helpers() };
    let helper_count = if direct {162} else {hs.len()};
    let (rows, counts) = scoped_rows(input, root);
    let selectors = if direct {9} else {5};
    let log_count = 5;
    let per_helper = selectors * 3 * 243;
    let count = log_count * helper_count * per_helper;
    if direct {
        json_write(root.join("grammar.json"),&json!({
            "log_providers":["frozen binary64 log1p","resident native FYL2XP1-ROMln2","stored native FYL2XP1-ROMln2","resident DLMF4.6.E6 PC64 unchanged-sum","stored DLMF4.6.E6"],
            "helpers":162,"helper_layout":"id=2*little_endian_base3(t_product,x_product,fraction_subtraction,em_sum)+f2xm1_binary64_store",
            "helper_programme":"t=nper*log1p; x=t*ROMlog2e; k=nearest_integer(x); f=F2XM1(x-k); em=scale(f,k)+(scale(1,k)-1)",
            "growth_selectors":"0..5: scale(round_model(1+f),k) with three models and final binary64 store boolean;6..8: round_model(1+em)",
            "models":MODELS,"sites":["PV*growth","tf addition","combine op1","combine op2","combine op3"],
            "combines":["((N/em)/tf)*rate","(N*rate)/(em*tf)","(N/tf)*(rate/em)"],
            "layout":"id=(((log_provider*162+helper)*9+growth_selector)*3+combine)*243+little_endian_base3_root_site_code",
            "programs":count,"scope_counts":counts,"instruction_policy":"FYL2XP1 is a native measured leaf in both modes; all other operations including F2XM1 use declared backend",
            "source":"Intel SDM F2XM1/FSCALE/FRNDINT/FLDL2E public ISA; annuity algebra; DLMF4.6.E6; current cleared Fpu",
            "constant_policy":"Existing64-bit ROM log2e/ln2, exact0/1/2, mathematical nearest-integer reduction; no fitted constants or boundary",
            "arithmetic_backend":std::any::type_name::<F>()}));
    } else {
    json_write(
        root.join("grammar.json"),
        &json!({"positive_expm1_helpers":hs,"log_providers":["frozen binary64 log1p","resident FYL2XP1 with64-bit ROM ln2","stored FYL2XP1","resident DLMF4.6.6 PC64 unchanged-sum loop","stored DLMF4.6.6 PC64 unchanged-sum loop"],
        "growth_selectors":["helper exponential u","binary64 store of u","strict 1+em","stored-PC64 1+em","resident-PC64 1+em"],
        "helper_components":"t=nper*log1p, u=exp(t), em=declared Kahan/difference programme; u/em share the same computed t/u",
        "sites":["PV*growth","tf addition","combine op1","combine op2","combine op3"],"models":MODELS,
        "combines":["((N/em)/tf)*rate","(N*rate)/(em*tf)","(N/tf)*(rate/em)"],
        "layout":"id=(((log_provider*1980+helper)*5+growth_selector)*3+combine)*243+little_endian_base3_site_code",
        "programs":count,"scope_counts":counts,"constant_policy":"Only0,1,2, exact positive odd integers and existing64-bit ROM ln2; no fitted coefficients or threshold",
        "arithmetic_backend":std::any::type_name::<F>(),"instruction_backend":"FYL2XP1 leaf measured on this AMD hardware for both modes; remaining operations use declared Fpu backend",
        "sources":["Intel SDM vol2A FYL2XP1 public ISA","https://dlmf.nist.gov/4.6.E6","annuity algebra","cleared Fpu"],
        "claim":"Tests continuous-log and coupled growth/annuity structure on PMT only; no claim that its growth helper is shared with the independently constrained PV/FV program"}),
    );
    }
    let mut alive: Vec<usize> = (0..count).collect();
    let mut evaluations = 0u64;
    let mut kills = BufWriter::new(fs::File::create(root.join("first-kills.jsonl")).unwrap());
    let mut stats = BufWriter::new(fs::File::create(root.join("witness-kills.jsonl")).unwrap());
    let mut killing = vec![];
    for (idx, row) in rows.iter().enumerate() {
        let before = alive.len();
        let mut rejected = std::collections::BTreeMap::<String, Vec<usize>>::new();
        if before > 0 {
            let log = logs::<F>(row.args()[0], row.log1p);
            let n = nper(row);
            let components: Vec<Vec<_>> = log
                .iter()
                .map(|&l| {
                    (0..helper_count)
                        .map(|h| if direct {direct_components::<F>(h,n,l)} else {boundary::helper_components::<F>(&hs[h], -n, l)})
                        .collect()
                })
                .collect();
            alive.retain(|&id| {
                let code = id % 243;
                let combine = (id / 243) % 3;
                let selector = (id / (243 * 3)) % selectors;
                let helper = (id / per_helper) % helper_count;
                let lp = id / (per_helper * helper_count);
                let (a, u, em) = components[lp][helper];
                let g = if direct {
                    if selector < 6 {
                        let one_f=binop::<F>(Op::Add(0,1),MODELS[selector/2],Val::f64(1.0),a);
                        store::<F>(Val::Ext(F::scale(&one_f.to_ext(),&u.to_ext(),CW)),selector%2!=0)
                    } else {
                        binop::<F>(Op::Add(0,1),MODELS[selector-6],Val::f64(1.0),em)
                    }
                } else { match selector {
                    0 => u,
                    1 => store::<F>(u, true),
                    2..=4 => binop::<F>(Op::Add(0, 1), MODELS[selector - 2], Val::f64(1.0), em),
                    _ => unreachable!(),
                }};
                let got = outcome::<F>(crate::annuity::pmt_growth_value::<F>(
                    row, g, em, combine, code,
                ));
                evaluations += 1;
                if got == row.expected() {
                    true
                } else {
                    rejected.entry(got).or_default().push(id);
                    false
                }
            });
            for (got, ids) in rejected {
                serde_json::to_writer(&mut kills,&json!({"candidate_ids":ids,"case":idx,"witness":row,"expected":row.expected(),"got":got})).unwrap();
                writeln!(kills).unwrap();
            }
        }
        let rec = json!({"case":idx,"id":row.raw["id"],"killed":before-alive.len(),"remaining":alive.len()});
        serde_json::to_writer(&mut stats, &rec).unwrap();
        writeln!(stats).unwrap();
        if before != alive.len() {
            println!("{rec}");
            killing.push(rec);
            kills.flush().unwrap();
            stats.flush().unwrap()
        }
    }
    json_write(root.join("survivors.json"), &alive);
    let summary = json!({"programs":count,"survivors":alive.len(),"scoped_rows":rows.len(),"scope_counts":counts,"evaluations":evaluations,"first_kill_rows":killing,"direct_F2XM1":direct,"complete":true,"claim":"Exact continuous-log/coupled-growth PMT exclusions or scoped survivors, not product matches or identification","FYL2XP1_policy":"Native measured leaf in both arithmetic backends; not a pure software replay of that instruction"});
    json_write(root.join("summary.json"), &summary);
    println!("{summary}");
}

fn nper_profile<F: Fpu>(a: [f64; 5], cash_model: EvalModel) -> String {
    let [r, p, pv, fv, t] = a;
    let stored = EvalModel::StoredPc64;
    if a.iter().any(|x| !x.is_finite()) {
        return "error:Value".into();
    }
    let bin = |op, a, b| binop::<F>(op, stored, a, b);
    if r == 0.0 {
        if p == 0.0 {
            return "error:Div0".into();
        }
        return outcome::<F>(bin(
            Op::Div(0, 1),
            boundary::neg(bin(Op::Add(0, 1), Val::f64(fv), Val::f64(pv))),
            Val::f64(p),
        ));
    }
    let tf = if t.trunc() == 0.0 {
        Val::f64(1.0)
    } else {
        bin(Op::Add(0, 1), Val::f64(1.0), Val::f64(r))
    };
    let tfp = bin(Op::Mul(0, 1), tf, Val::f64(p));
    let fr = binop::<F>(Op::Mul(0, 1), cash_model, Val::f64(fv), Val::f64(r));
    let pr = binop::<F>(Op::Mul(0, 1), cash_model, Val::f64(pv), Val::f64(r));
    let numerator = bin(Op::Add(0, 1), tfp, boundary::neg(fr));
    let denominator = bin(Op::Add(0, 1), tfp, pr);
    let ratio = bin(Op::Div(0, 1), numerator, denominator);
    let rr = ratio.to_f64::<F>();
    if rr <= 0.0 || !rr.is_finite() {
        return "error:Num".into();
    }
    let base = bin(Op::Add(0, 1), Val::f64(1.0), Val::f64(r));
    if base.to_f64::<F>() <= 0.0 {
        return "error:Num".into();
    }
    let lnbase = store::<F>(Val::Ext(research_log::<F>(&base.to_ext())), true);
    if lnbase.to_f64::<F>() == 0.0 {
        return "error:Div0".into();
    }
    let lnratio = store::<F>(Val::Ext(research_log::<F>(&ratio.to_ext())), true);
    outcome::<F>(bin(Op::Div(0, 1), lnratio, lnbase))
}

pub fn nper_product_profiles(input: &str, out: &str, full_five: bool) {
    let rows = crate::lines::<serde_json::Value>(input);
    let arithmetic = if full_five {
        vec![
            EvalModel::Strict,
            EvalModel::StoredPc64,
            EvalModel::ResidentPc64,
            EvalModel::StoredPc53,
            EvalModel::ResidentPc53,
        ]
    } else {
        MODELS.to_vec()
    };
    let mut alive: Vec<usize> = (0..arithmetic.len()).collect();
    let mut first = vec![];
    let mut checks = vec![0usize; arithmetic.len()];
    for (idx, row) in rows.iter().enumerate() {
        assert_eq!(row["function"], "NPER");
        assert_eq!(row["args"].as_array().unwrap().len(), 5);
        let a: [f64; 5] = std::array::from_fn(|j| {
            let v = calc_verdict::decode_argument(&row["args"][j], "NPER raw bits").unwrap();
            match v.core() {
                crate::value::CoreValue::Number(x) => *x,
                _ => panic!("NPER typed control outside declared numeric bank"),
            }
        });
        let got = nper_profile::<Hard>(a, EvalModel::StoredPc64);
        assert_eq!(
            got,
            row["expected"].as_str().unwrap(),
            "stored baseline must reproduce every original application reading, {}",
            row["id"]
        );
        alive.retain(|&id|{checks[id]+=1;let got=nper_profile::<Hard>(a,arithmetic[id]);
            if got==row["expected"].as_str().unwrap(){true}else{first.push(json!({"model":arithmetic[id],"case":idx,"witness":row,"expected":row["expected"],"got":got}));false}
        });
    }
    let allowed: Vec<_> = alive.iter().map(|&i| arithmetic[i]).collect();
    let result = json!({"rows":rows.len(),"baseline_application_trace_matches":rows.len(),"cash_product_models":arithmetic,"allowed_model_indices":alive,"allowed_models":allowed,"first_kills":first,"evaluations_per_model":checks,"backend":"native Fpu, including existing x87 logarithms","sharing":"Only PV*rate and FV*rate change model; all other NPER sites retain the existing registered graph. The same interest-product model is imposed on the PMT I=PV*rate hypothesis.","claim":"Exact conditional joint constraint, not uniqueness of the NPER programme or a product gain"});
    json_write(out, &result);
    println!("{result}");
}

fn interest_value<F: Fpu>(row: &Row, em: Val, mut code: usize) -> Val {
    let m: [EvalModel; 5] = std::array::from_fn(|_| {
        let v = MODELS[code % 3];
        code /= 3;
        v
    });
    let [rate, pv, _, timing, _] = row.args();
    let i = binop::<F>(Op::Mul(0, 1), m[0], Val::f64(pv), Val::f64(rate));
    let tf = if timing == 0.0 {
        Val::f64(1.0)
    } else {
        binop::<F>(Op::Add(0, 1), m[1], Val::f64(1.0), Val::f64(rate))
    };
    let amort = binop::<F>(Op::Div(0, 1), m[2], i, em);
    let total = binop::<F>(Op::Add(0, 1), m[3], i, amort);
    boundary::neg(binop::<F>(Op::Div(0, 1), m[4], total, tf))
}
pub fn interest_search<F: Fpu>(input: &str, nper_profiles: &str, out: &str) {
    fs::create_dir_all(out).unwrap();
    let root = Path::new(out);
    let hs = boundary::helpers();
    let (rows, counts) = scoped_rows(input, root);
    let precondition: serde_json::Value =
        serde_json::from_slice(&fs::read(nper_profiles).unwrap()).unwrap();
    let allowed: Vec<usize> =
        serde_json::from_value(precondition["allowed_model_indices"].clone()).unwrap();
    assert!(!allowed.is_empty());
    let mut alive: Vec<usize> = (0..5 * hs.len() * 243)
        .filter(|id| allowed.contains(&(id % 3)))
        .collect();
    let count = alive.len();
    let mut evaluations = 0u64;
    json_write(
        root.join("grammar.json"),
        &json!({"log_providers":["frozen_binary64_log1p","residentFYL2XP1","storedFYL2XP1","residentDLMF4.6.6_PC64_unchanged_sum","storedDLMF4.6.6_PC64_unchanged_sum"],"positive_expm1_helpers":hs,
        "NPER_precondition":precondition,"shared_product":"I=PV*rate, constrained by the same rate-weighted-cashflow product model in NPER",
        "programme":"I=PV*rate; amort=I/em; total=I+amort; output=-total/tf","sharing":"I is computed once and reused unchanged on both edges; implicit strict input stores remain local to each consumer",
        "sites":["I product","tf addition","I/em","I+amort","total/tf"],"models":MODELS,
        "layout":"id=(log_provider*1980+helper)*243+little_endian_base3_site_code; only NPER-allowed first digits admitted",
        "programs":count,"scope_counts":counts,"FYL2XP1_policy":"native leaf in both arithmetic modes, no pure software instruction emulation claim",
        "source":"Annuity decomposition -(PV*rate+(PV+FV)*rate/expm1(n*log1p(rate)))/tf; FV=0 arm uses shared I. This explains the observed opposing-balance identity without selecting answer residuals."}),
    );
    let mut kills = BufWriter::new(fs::File::create(root.join("first-kills.jsonl")).unwrap());
    let mut stats = BufWriter::new(fs::File::create(root.join("witness-kills.jsonl")).unwrap());
    let mut killing = vec![];
    for (idx, row) in rows.iter().enumerate() {
        let before = alive.len();
        let mut rejected = std::collections::BTreeMap::<String, Vec<usize>>::new();
        if before > 0 {
            let l = logs::<F>(row.args()[0], row.log1p);
            let n = nper(row);
            let em: Vec<Vec<_>> = l
                .iter()
                .map(|&log| {
                    hs.iter()
                        .map(|h| boundary::helper_components::<F>(h, -n, log).2)
                        .collect()
                })
                .collect();
            alive.retain(|&id| {
                let lp = id / (hs.len() * 243);
                let h = (id / 243) % hs.len();
                let got = outcome::<F>(interest_value::<F>(row, em[lp][h], id % 243));
                evaluations += 1;
                if got == row.expected() {
                    true
                } else {
                    rejected.entry(got).or_default().push(id);
                    false
                }
            });
            for (got, ids) in rejected {
                serde_json::to_writer(&mut kills,&json!({"candidate_ids":ids,"case":idx,"witness":row,"expected":row.expected(),"got":got})).unwrap();
                writeln!(kills).unwrap();
            }
        }
        let rec = json!({"case":idx,"id":row.raw["id"],"killed":before-alive.len(),"remaining":alive.len()});
        serde_json::to_writer(&mut stats, &rec).unwrap();
        writeln!(stats).unwrap();
        if before != alive.len() {
            println!("{rec}");
            killing.push(rec);
            kills.flush().unwrap();
            stats.flush().unwrap()
        }
    }
    json_write(root.join("survivors.json"), &alive);
    let summary = json!({"programs":count,"survivors":alive.len(),"scoped_rows":rows.len(),"scope_counts":counts,"evaluations":evaluations,"first_kill_rows":killing,"complete":true,"claim":"Exact NPER/PMT interest-product joint constraint and FV=0 decomposition search, not identification or a product pass count"});
    json_write(root.join("summary.json"), &summary);
    println!("{summary}");
}
