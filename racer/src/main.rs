#![allow(dead_code, unused_imports)]
#[cfg(target_arch = "wasm32")]
compile_error!("symbolic-racer is native-only");

pub use calc_func::value;
mod ported {
    include!(concat!(env!("OUT_DIR"), "/bridge.rs"));
}
mod boundary;
mod dsl;
mod enumerate;
mod eval;
mod identity;
mod poisson;
mod saddle;
mod scheduler;

use ported::excel_numeric::{Fpu, Hard, Soft};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::{
    fs,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};

#[derive(Serialize, Deserialize, Clone)]
pub struct Row {
    pub trace_index: usize,
    pub raw: Value,
    pub leaves: [u64; 5],
    pub log1p: u64,
    pub neg_log: u64,
    pub prior_gap: bool,
}
impl Row {
    pub fn args(&self) -> [f64; 5] {
        self.leaves.map(f64::from_bits)
    }
    pub fn expected(&self) -> &str {
        self.raw["expected"].as_str().unwrap()
    }
}

fn json_write(path: impl AsRef<Path>, data: &impl Serialize) {
    let mut f = BufWriter::new(fs::File::create(path).unwrap());
    serde_json::to_writer_pretty(&mut f, data).unwrap();
    writeln!(f).unwrap();
}
fn lines<T: for<'de> Deserialize<'de>>(path: impl AsRef<Path>) -> Vec<T> {
    BufReader::new(fs::File::open(path).unwrap())
        .lines()
        .map(|l| serde_json::from_str(&l.unwrap()).unwrap())
        .collect()
}
fn count(c: &mut std::collections::BTreeMap<String, usize>, key: &str) {
    *c.entry(key.into()).or_default() += 1;
}
struct Given(Vec<value::CalcValue>);
impl calc_func::PreparedArguments for Given {
    fn count(&self) -> usize {
        self.0.len()
    }
    fn argument(&mut self, index: usize) -> value::CalcValue {
        self.0[index].clone()
    }
}
fn invoke_control(function: &str, args: Vec<value::CalcValue>) -> (String, Option<f64>) {
    let target = calc_func::catalogue()
        .resolve_surface(calc_func::Notation::Call, function)
        .unwrap();
    let mode = value::NumericMode::Fast;
    match calc_func::invoke(
        target,
        &mut Given(args),
        &mut calc_func::EvalServices::new(&mode),
        &mut calc_func::FunctionCallScratch::new(),
    ) {
        Ok(v) => {
            let n = match v.core() {
                value::CoreValue::Number(n) => Some(*n),
                _ => None,
            };
            (calc_verdict::encode_result(&v).unwrap(), n)
        }
        Err(f) => {
            let s = format!("{f:?}");
            (
                format!("fault:{}", s.split([' ', '{', '(']).next().unwrap()),
                None,
            )
        }
    }
}

fn prepare(trace: &str, gaps: &str, out: &str) {
    fs::create_dir_all(out).unwrap();
    let gaps: std::collections::HashSet<String> = lines::<Value>(gaps)
        .into_iter()
        .map(|r| r["id"].as_str().unwrap().into())
        .collect();
    let mut counts = std::collections::BTreeMap::new();
    let mut rows = Vec::new();
    let mut exclusions =
        BufWriter::new(fs::File::create(Path::new(out).join("excluded.jsonl")).unwrap());
    let mut baseline = enumerate::template(0, 0);
    baseline.nodes[12].model = dsl::EvalModel::StoredPc64;
    baseline.nodes[13].model = dsl::EvalModel::StoredPc64;
    for (i, line) in BufReader::new(fs::File::open(trace).unwrap())
        .lines()
        .enumerate()
    {
        let raw: Value = serde_json::from_str(&line.unwrap()).unwrap();
        count(&mut counts, "recorded_rows");
        assert_eq!(raw["function"], "PMT");
        let a = raw["args"].as_array().unwrap();
        let decoded: Result<Vec<_>, _> = a
            .iter()
            .enumerate()
            .map(|(j, v)| calc_verdict::decode_argument(v, &format!("args[{j}]")))
            .collect();
        let numbers: Option<Vec<f64>> = decoded.ok().and_then(|v| {
            v.into_iter()
                .map(|x| match x.core() {
                    value::CoreValue::Number(n) => Some(*n),
                    _ => None,
                })
                .collect()
        });
        let mut leaf = None;
        let reason = if let Some(a) = numbers {
            if a.len() < 3 || a.len() > 5 || !raw["result_index"].is_null() {
                "unsupported_input_shape"
            } else if a.iter().any(|x| !x.is_finite() || x.is_subnormal()) {
                "input_admission_refusal"
            } else {
                let r = a[0];
                let n = a[1];
                let pv = a[2];
                let fv = *a.get(3).unwrap_or(&0.0);
                let timing = if a.get(4).unwrap_or(&0.0).trunc() == 0.0 {
                    0.0
                } else {
                    1.0
                };
                if r == 0.0 {
                    "zero_rate"
                } else if n == 0.0 || 1.0 + r <= 0.0 {
                    "domain_error"
                } else {
                    let log = ported::excel_numeric::excel_log1p(r);
                    let t = -(n * log);
                    let em = ported::excel_numeric::excel_expm1_spilled(t);
                    if !t.is_finite() || !em.is_finite() {
                        "nonfinite_exponential"
                    } else if em == 0.0 {
                        "zero_expm1"
                    } else if pv != 0.0 && fv == -pv {
                        "opposing_balances"
                    } else if t >= 0.0 {
                        "growth_control"
                    } else {
                        leaf = Some(([r, pv, fv, timing, em], log, t));
                        "scoped"
                    }
                }
            }
        } else {
            "unsupported_argument_type"
        };
        count(&mut counts, reason);
        if let Some((args, log, t)) = leaf {
            // Unknown expected types are an explicit failed load, never skipped.
            calc_verdict::decode_result(raw["expected"].as_str().unwrap())
                .expect("scoped expected outcome must decode");
            let actual = eval::outcome::<Hard>(eval::eval_graph::<Hard>(&baseline, &args));
            assert_eq!(
                actual,
                raw["actual"].as_str().unwrap(),
                "baseline route mismatch on {}",
                raw["id"]
            );
            rows.push(Row {
                trace_index: i,
                prior_gap: gaps.contains(raw["id"].as_str().unwrap()),
                raw,
                leaves: args.map(f64::to_bits),
                log1p: log.to_bits(),
                neg_log: t.to_bits(),
            });
        } else {
            serde_json::to_writer(
                &mut exclusions,
                &json!({"trace_index":i,"reason":reason,"raw":raw}),
            )
            .unwrap();
            writeln!(exclusions).unwrap();
        }
    }
    rows.sort_by_key(|r| (!r.prior_gap, r.trace_index));
    assert_eq!(
        rows.iter().filter(|r| r.prior_gap).count(),
        734,
        "all prior gaps must be covered"
    );
    let mut f = BufWriter::new(fs::File::create(Path::new(out).join("scoped.jsonl")).unwrap());
    for r in &rows {
        serde_json::to_writer(&mut f, r).unwrap();
        writeln!(f).unwrap();
    }
    counts.insert("baseline_exact_route_replays".into(), rows.len());
    json_write(
        Path::new(out).join("coverage.json"),
        &json!({"counts":counts,"prior_gaps_first":734,
        "discovery_rows":32,"order":"prior-gap boolean then original trace index; first 32 are discovery, all are previously inspected evidence",
        "unsupported_or_unread_are_passes":false,"fresh_validation_rows":0}),
    );
    println!("{}", json!(counts));
}

fn self_test(rows: &[Row]) {
    use dsl::{MODELS, Op};
    use eval::{Val, binop};
    let target = calc_func::catalogue()
        .resolve_surface(calc_func::Notation::Call, "PMT")
        .unwrap();
    let mode = value::NumericMode::Fast;
    for (present, future, periods, raw) in [
        (f64::MIN_POSITIVE, 0.0, 4.0, -f64::MIN_POSITIVE / 4.0),
        (0.0, 0.0, 1.0, -0.0),
        (f64::MAX, f64::MAX, 1.0, f64::NEG_INFINITY),
    ] {
        let actual = calc_func::invoke(
            target,
            &mut Given(
                [0.0, periods, present, future, 0.0]
                    .into_iter()
                    .map(value::CalcValue::number)
                    .collect(),
            ),
            &mut calc_func::EvalServices::new(&mode),
            &mut calc_func::FunctionCallScratch::new(),
        );
        let encoded = match actual {
            Ok(value) => calc_verdict::encode_result(&value).unwrap(),
            Err(fault) => {
                let text = format!("{fault:?}");
                format!("fault:{}", text.split([' ', '{', '(']).next().unwrap())
            }
        };
        assert_eq!(
            encoded,
            eval::outcome::<Hard>(Val::f64(raw)),
            "application publication"
        );
        assert_eq!(
            encoded,
            eval::bits_outcome(eval::publication_bits::<Hard>(Val::f64(raw)))
        );
    }
    let family = enumerate::family();
    let mut state = 0x504d545f5236u64;
    let mut checks = 0;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    for _ in 0..20000 {
        let a = f64::from_bits((next() & 0x800fffffffffffff) | (((next() % 2046) + 1) << 52));
        let b = f64::from_bits((next() & 0x800fffffffffffff) | (((next() % 2046) + 1) << 52));
        for op in [Op::Add(0, 1), Op::Mul(0, 1), Op::Div(0, 1)] {
            for model in MODELS {
                let av = Val::f64(a);
                let bv = Val::f64(b);
                let h = binop::<Hard>(op, model, av, bv);
                let s = binop::<Soft>(op, model, av, bv);
                if h.to_f64::<Hard>().is_nan() {
                    assert!(s.to_f64::<Soft>().is_nan());
                } else {
                    assert_eq!(h, s, "arithmetic parity {op:?} {model:?} {a:?} {b:?}");
                }
                // A resident operation feeding a second operation must retain its
                // full bits, not the binary64 view used for final publication.
                let h = binop::<Hard>(Op::Add(0, 1), model, h, Val::f64(-a));
                let s = binop::<Soft>(Op::Add(0, 1), model, s, Val::f64(-a));
                if h.to_f64::<Hard>().is_nan() {
                    assert!(s.to_f64::<Soft>().is_nan());
                } else {
                    assert_eq!(h, s, "resident continuation parity");
                }
                checks += 2;
            }
        }
    }
    // Exact set propagation is independently checked against every serialized
    // graph, including typed overflow/refusal and signed-zero publication.
    let sample: Vec<_> = rows
        .iter()
        .take(32)
        .chain((0..32).map(|i| &rows[i * (rows.len() - 1) / 31]))
        .collect();
    for r in &sample {
        let mut h = std::collections::BTreeSet::new();
        let mut s = std::collections::BTreeSet::new();
        let mut hr = std::collections::BTreeSet::new();
        let mut sr = std::collections::BTreeSet::new();
        for c in &family {
            let hv = eval::eval_graph::<Hard>(&c.graph, &r.args());
            let sv = eval::eval_graph::<Soft>(&c.graph, &r.args());
            h.insert(eval::publication_bits::<Hard>(hv));
            s.insert(eval::publication_bits::<Soft>(sv));
            hr.insert(hv.exact_encoding());
            sr.insert(sv.exact_encoding());
        }
        assert_eq!(h, s, "full family backend parity");
        assert_eq!(
            h.into_iter().collect::<Vec<_>>(),
            eval::reachable::<Hard>(&r.args())
        );
        assert_eq!(
            s.into_iter().collect::<Vec<_>>(),
            eval::reachable::<Soft>(&r.args())
        );
        assert_eq!(hr, sr, "full resident values, both backends");
        assert_eq!(
            hr,
            eval::root_values::<Hard>(&r.args())
                .into_iter()
                .map(|v| v.exact_encoding())
                .collect()
        );
        assert_eq!(
            sr,
            eval::root_values::<Soft>(&r.args())
                .into_iter()
                .map(|v| v.exact_encoding())
                .collect()
        );
    }
    println!(
        "self-test passed: 3 actual invoke-route publication controls; {checks} raw arithmetic/continuation parity checks; {} serialized programs x {} rows both backends; exact tagged root and publication set equality",
        family.len(),
        sample.len()
    );
}

fn main() {
    let a: Vec<_> = std::env::args().collect();
    match a.get(1).map(String::as_str) {
        Some("poisson-publication-test") => poisson::publication_controls(),
        Some("normalizer-matrix") => {
            let rows = lines::<poisson::MassRow>(&a[2]);
            match a.get(5).map(String::as_str) {
                Some("soft") => saddle::normalizer_matrix::<Soft>(&rows, &a[3], &a[4]),
                _ => saddle::normalizer_matrix::<Hard>(&rows, &a[3], &a[4]),
            }
        }
        Some("questions") => saddle::questions(
            &lines::<poisson::MassRow>(&a[2]),
            &a[3],
            &a[4],
            &a[5],
            &a[6],
        ),
        Some("count-one-builds") => {
            let rows = lines::<poisson::MassRow>(&a[2]);
            match a.get(5).map(String::as_str) {
                Some("soft") => saddle::build_cohorts::<Soft>(&rows, &a[3], &a[4]),
                _ => saddle::build_cohorts::<Hard>(&rows, &a[3], &a[4]),
            }
        }
        Some("precision-test") => eval::precision_controls(),
        Some("count-one-full") => {
            let rows = lines::<poisson::MassRow>(&a[2]);
            match a.get(5).map(String::as_str) {
                Some("soft") => saddle::run_full::<Soft>(&rows, &a[3], &a[4]),
                _ => saddle::run_full::<Hard>(&rows, &a[3], &a[4]),
            }
        }
        Some("count-one-full-test") => saddle::test_full(&lines::<poisson::MassRow>(&a[2]), &a[3]),
        Some("count-one-refined") => {
            let rows = lines::<poisson::MassRow>(&a[2]);
            match a.get(5).map(String::as_str) {
                Some("soft") => saddle::run_refined::<Soft>(&rows, &a[3], &a[4]),
                _ => saddle::run_refined::<Hard>(&rows, &a[3], &a[4]),
            }
        }
        Some("count-one-refined-test") => {
            saddle::test_refined(&lines::<poisson::MassRow>(&a[2]), &a[3])
        }
        Some("boundary-wide") => {
            let rows = lines::<Row>(&a[2]);
            match a.get(4).map(String::as_str) {
                Some("soft") => boundary::run_wide::<Soft>(&rows, &a[3]),
                _ => boundary::run_wide::<Hard>(&rows, &a[3]),
            }
        }
        Some("boundary-wide-test") => boundary::test_wide(&lines::<Row>(&a[2])),
        Some("count-one-expanded") => {
            let rows = lines::<poisson::MassRow>(&a[2]);
            match a.get(5).map(String::as_str) {
                Some("soft") => saddle::run_expanded::<Soft>(&rows, &a[3], &a[4]),
                _ => saddle::run_expanded::<Hard>(&rows, &a[3], &a[4]),
            }
        }
        Some("count-one-expanded-test") => {
            saddle::test_expanded(&lines::<poisson::MassRow>(&a[2]), &a[3])
        }
        Some("count-one") => {
            let rows = lines::<poisson::MassRow>(&a[2]);
            match a.get(5).map(String::as_str) {
                Some("soft") => saddle::run::<Soft>(&rows, &a[3], &a[4]),
                _ => saddle::run::<Hard>(&rows, &a[3], &a[4]),
            }
        }
        Some("count-one-test") => saddle::self_test(&lines::<poisson::MassRow>(&a[2]), &a[3]),
        Some("inverse-controls") => boundary::inverse_controls(&lines::<Row>(&a[2]), &a[3]),
        Some("mass-prepare") => poisson::prepare(&a[2], &a[3]),
        Some("mass-test") => poisson::self_test(&lines::<poisson::MassRow>(&a[2])),
        Some("mass") => {
            let rows = lines::<poisson::MassRow>(&a[2]);
            match a.get(4).map(String::as_str) {
                Some("soft") => poisson::run::<Soft>(&rows, &a[3]),
                _ => poisson::run::<Hard>(&rows, &a[3]),
            }
        }
        Some("identity") => {
            let rows = lines::<Row>(&a[2]);
            match a.get(4).map(String::as_str) {
                Some("soft") => identity::run::<Soft>(&rows, &a[3]),
                _ => identity::run::<Hard>(&rows, &a[3]),
            }
        }
        Some("boundary") => {
            let rows = lines::<Row>(&a[2]);
            match a.get(4).map(String::as_str) {
                Some("soft") => boundary::run::<Soft>(&rows, &a[3]),
                _ => boundary::run::<Hard>(&rows, &a[3]),
            }
        }
        Some("boundary-test") => boundary::self_test(&lines::<Row>(&a[2])),
        Some("prepare") => prepare(&a[2], &a[3], &a[4]),
        Some("self-test") => self_test(&lines::<Row>(&a[2])),
        Some("run") => {
            let rows = lines::<Row>(&a[2]);
            let output = &a[3];
            match a.get(4).map(String::as_str) {
                Some("soft") => scheduler::run::<Soft>(&rows, output),
                _ => scheduler::run::<Hard>(&rows, output),
            }
        }
        _ => panic!("prepare TRACE GAPS OUT | self-test SCOPED | run SCOPED OUT [hard|soft]"),
    }
}
