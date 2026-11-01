//! R11: monotone balance recurrences with every latent payment in its publication cell.
//! This identifies constraints on a consumer, never a payment-producing program.
use crate::{dsl::*, eval::*, ported::excel_numeric::Fpu};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    io::{BufWriter, Write},
    path::Path,
    time::Instant,
};

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct Profile {
    id: usize,
    update: usize,
    payment_store: bool,
    interest_model: EvalModel,
    principal_model: EvalModel,
    balance_model: EvalModel,
    sum_order: usize,
    sum_model: EvalModel,
}

fn profiles() -> Vec<Profile> {
    let mut result = Vec::new();
    for update in 0..2 {
        for payment_store in [false, true] {
            for interest_model in WIDE_MODELS {
                for principal_model in WIDE_MODELS {
                    for balance_model in WIDE_MODELS {
                        for sum_order in 0..3 {
                            for sum_model in WIDE_MODELS {
                                result.push(Profile {
                                    id: result.len(),
                                    update,
                                    payment_store,
                                    interest_model,
                                    principal_model,
                                    balance_model,
                                    sum_order,
                                    sum_model,
                                });
                            }
                        }
                    }
                }
            }
        }
    }
    assert_eq!(result.len(), 7500);
    result
}

pub fn graph(p: &Profile, start: usize, end: usize, beginning: bool) -> Graph {
    graph_variant(p, start, end, beginning, false)
}

fn graph_variant(p: &Profile, start: usize, end: usize, beginning: bool, interest_sum: bool) -> Graph {
    assert!(1 <= start && start <= end && end <= 512);
    let mut g = GraphBuilder::default();
    let rate = g.strict(Op::Arg(0));
    let pv = g.strict(Op::Arg(1));
    let latent = g.strict(Op::Arg(2));
    let one = g.strict(Op::Const {
        bits: 1f64.to_bits(),
    });
    let zero = g.strict(Op::Const { bits: 0 });
    let payment = if p.payment_store {
        g.strict(Op::Mul(latent, one))
    } else {
        latent
    };
    let q = if p.update == 1 {
        g.push(Op::Add(one, rate), p.balance_model)
    } else {
        one
    };
    let mut balance = pv;
    let mut terms = vec![];
    for period in 1..=end {
        let interest = if beginning && period == 1 { zero }
            else { g.push(Op::Mul(balance, rate), p.interest_model) };
        let principal = if beginning && period == 1 {
            payment
        } else {
            g.push(Op::Add(payment, interest), p.principal_model)
        };
        if period >= start {
            terms.push(if interest_sum { interest } else { principal });
        }
        if period < end {
            balance = if p.update == 0 || beginning && period == 1 {
                g.push(Op::Add(balance, principal), p.balance_model)
            } else {
                let grown = g.push(Op::Mul(balance, q), p.balance_model);
                g.push(Op::Add(grown, payment), p.balance_model)
            };
        }
    }
    let mut output = sum_graph(&mut g, &terms, p.sum_order, p.sum_model, zero);
    if interest_sum {
        let count = g.strict(Op::Const { bits: ((end-start+1) as f64).to_bits() });
        let payments = g.push(Op::Mul(count, payment), p.principal_model);
        output = g.push(Op::Add(payments, output), p.principal_model);
    }
    g.finish(output)
}

fn line(out: &mut impl Write, row: &impl Serialize) {
    serde_json::to_writer(&mut *out, row).unwrap();
    writeln!(out).unwrap();
}

fn input(row: &Value) -> (f64, f64, usize, usize, bool) {
    assert_eq!(row["function"], "CUMPRINC");
    assert!(row["result_index"].is_null());
    let a = row["args"].as_array().unwrap();
    assert_eq!(a.len(), 6);
    let a: Vec<_> = a.iter().map(crate::annuity::number).collect();
    assert!(
        a.iter()
            .all(|x| x.is_finite() && (*x == 0.0 || x.is_normal()))
    );
    assert!(a[0] > 0.0 && a[2] > 0.0 && matches!(a[5], 0.0 | 1.0));
    assert!(a[1].fract() == 0.0 && a[3].fract() == 0.0 && a[4].fract() == 0.0);
    assert!(1.0 <= a[3] && a[3] <= a[4] && a[4] <= a[1] && a[1] <= 512.0);
    (a[0], a[2], a[3] as usize, a[4] as usize, a[5] != 0.0)
}

struct Loan {
    raw: Value,
    indices: Vec<usize>,
    latent: Vec<Val>,
}

pub fn run<F: Fpu>(trace: &str, structure: &str, out: &str) {
    run_with_domains::<F>(trace, structure, out, None)
}

pub fn run_with_domains<F: Fpu>(trace: &str, structure: &str, out: &str, domains: Option<&str>) {
    run_variant::<F>(trace, structure, out, domains, false)
}

pub fn run_variant<F: Fpu>(trace: &str, structure: &str, out: &str, domains: Option<&str>, interest_sum: bool) {
    let clock = Instant::now();
    let rows = crate::lines::<Value>(trace);
    assert_eq!(rows.len(), 540);
    for row in &rows {
        let _ = input(row);
    }
    let original: Value = serde_json::from_slice(&fs::read(structure).unwrap()).unwrap();
    let groups = original["groups"].as_array().unwrap();
    assert_eq!(groups.len(), 60);
    let free_domains: Option<Value> = domains.map(|path| {
        serde_json::from_slice(&fs::read(path).unwrap()).unwrap()
    });
    let mut domain_map = BTreeMap::new();
    if let Some(doc) = &free_domains {
        for d in doc["domains"].as_array().unwrap() {
            let mut parts = vec![d["loan_index"].clone(), d["payment_store"].clone(),
                d["interest_model"].clone(), d["principal_model"].clone(), d["sum_model"].clone()];
            if interest_sum { parts.push(d["sum_order"].clone()); }
            let key = serde_json::to_string(&parts).unwrap();
            assert!(domain_map.insert(key, d).is_none());
        }
        assert_eq!(domain_map.len(), if interest_sum {45000} else {15000});
        assert_eq!(doc["interest_sum"].as_bool().unwrap_or(false), interest_sum);
    }
    let mut loans = vec![];
    let mut seen = BTreeSet::new();
    fs::create_dir_all(out).unwrap();
    let root = Path::new(out);
    if let Some(doc) = &free_domains {
        crate::json_write(root.join("initial-domains.json"), doc);
    }
    let mut lf = BufWriter::new(fs::File::create(root.join("loans.jsonl")).unwrap());
    for (loan_index, group) in groups.iter().enumerate() {
        let raw = group["PMT"].clone();
        assert_eq!(raw["function"], "PMT");
        assert!(raw["result_index"].is_null());
        let pa = raw["args"].as_array().unwrap();
        assert_eq!(pa.len(), 5);
        assert_eq!(crate::annuity::number(&pa[3]), 0.0);
        let indices: Vec<_> = group["row_indices"]
            .as_array()
            .unwrap()
            .iter()
            .map(|i| i.as_u64().unwrap() as usize)
            .collect();
        assert_eq!(indices.len(), 9);
        for &i in &indices {
            assert!(seen.insert(i));
            let a = rows[i]["args"].as_array().unwrap();
            assert_eq!(&a[..3], &pa[..3]);
            assert_eq!(a[5], pa[4]);
            assert_eq!(rows[i]["reading"], raw["reading"]);
        }
        let bits = u64::from_str_radix(
            raw["expected"].as_str().unwrap().trim_start_matches("0x"),
            16,
        )
        .unwrap();
        let latent = crate::annuity_early::publication_candidates::<F>(bits);
        line(
            &mut lf,
            &json!({"loan_index":loan_index,"payment":raw,"row_indices":indices,
            "latent_payment_values":latent.iter().map(|p|p.exact_encoding()).collect::<Vec<_>>(),
            "meaning":"complete existential Ext80 publication cell, not a fitted or implemented payment"}),
        );
        loans.push(Loan {
            raw,
            indices,
            latent,
        });
    }
    assert_eq!(
        seen.into_iter().collect::<Vec<_>>(),
        (0..540).collect::<Vec<_>>()
    );
    lf.flush().unwrap();
    let mut pf = BufWriter::new(fs::File::create(root.join("profiles.jsonl")).unwrap());
    let mut cf = BufWriter::new(fs::File::create(root.join("constraints.jsonl")).unwrap());
    let mut gf = BufWriter::new(fs::File::create(root.join("graphs.jsonl")).unwrap());
    let mut tf = BufWriter::new(fs::File::create(root.join("boundary-traces.jsonl")).unwrap());
    let mut sf = BufWriter::new(fs::File::create(root.join("survivors.jsonl")).unwrap());
    let mut completed = 0;
    let mut killed = 0;
    let mut surviving = 0;
    let mut evaluations = 0u64;
    let mut trace_count = 0usize;
    let mut histogram = BTreeMap::<usize, usize>::new();
    for p in profiles() {
        // A stopped budget produces an explicit partial result, never zero survivors.
        if clock.elapsed().as_secs() > 5400 {
            break;
        }
        line(&mut pf, &p);
        let mut first_kill = None;
        let mut feasible_loans = vec![];
        for (loan_index, loan) in loans.iter().enumerate() {
            let free_latent = free_domains.as_ref().map(|_| {
                let mut parts = vec![json!(loan_index), json!(p.payment_store),
                    json!(p.interest_model), json!(p.principal_model), json!(p.sum_model)];
                if interest_sum { parts.push(json!(p.sum_order)); }
                let key = serde_json::to_string(&parts).unwrap();
                let d = domain_map[&key];
                assert_eq!(d["row_index"].as_u64().unwrap() as usize, loan.indices[0]);
                if d["count"].as_u64().unwrap() == 0 {
                    assert!(d["lower"].is_null() && d["upper"].is_null());
                    return Vec::new();
                }
                fn extended(raw: &Value) -> crate::ported::excel_numeric::Extended {
                    let s = raw.as_str().unwrap().strip_prefix("ext80-le:").unwrap();
                    assert_eq!(s.len(), 20);
                    let mut bytes = [0u8; 10];
                    for (i,b) in bytes.iter_mut().enumerate() {
                        *b = u8::from_str_radix(&s[2*i..2*i+2],16).unwrap();
                    }
                    crate::ported::excel_numeric::Extended(bytes)
                }
                let lower = extended(&d["lower"]);
                let upper = extended(&d["upper"]);
                assert_eq!(&lower.0[8..], &upper.0[8..]);
                assert_ne!(lower.0[9] & 128, 0);
                let a = u64::from_le_bytes(lower.0[..8].try_into().unwrap());
                let b = u64::from_le_bytes(upper.0[..8].try_into().unwrap());
                let count = d["count"].as_u64().unwrap();
                assert!(a >= b && a-b+1 == count && count <= 100000);
                (b..=a).rev().map(|sig| {
                    let mut value = lower;
                    value.0[..8].copy_from_slice(&sig.to_le_bytes());
                    Val::Ext(value)
                }).collect::<Vec<_>>()
            });
            let latent = free_latent.as_ref().unwrap_or(&loan.latent);
            if latent.is_empty() {
                first_kill = Some((loan_index, loan.indices[0]));
                break;
            }
            let mut low = 0;
            let mut high = latent.len();
            for &row_index in &loan.indices {
                let (rate, pv, start, end, beginning) = input(&rows[row_index]);
                let g = graph_variant(&p, start, end, beginning, interest_sum);
                let target = crate::annuity::number(&rows[row_index]["expected"]);
                let before = (low, high);
                let mut cache = BTreeMap::<usize, f64>::new();
                let mut value = |index: usize| -> f64 {
                    *cache.entry(index).or_insert_with(|| {
                        evaluations += 1;
                        let args = [
                            Val::f64(rate),
                            Val::f64(pv),
                            latent[index],
                            Val::f64(0.0),
                            Val::f64(0.0),
                        ];
                        let got = eval_graph_typed::<F>(&g, &args).to_f64::<F>();
                        assert!(got.is_finite());
                        got
                    })
                };
                let bound = |strict: bool,
                             mut a: usize,
                             mut b: usize,
                             value: &mut dyn FnMut(usize) -> f64| {
                    while a < b {
                        let middle = a + (b - a) / 2;
                        let y = value(middle);
                        if if strict { y > target } else { y >= target } {
                            b = middle;
                        } else {
                            a = middle + 1;
                        }
                    }
                    a
                };
                let left = bound(false, low, high, &mut value);
                let right = bound(true, left, high, &mut value);
                low = left;
                high = right;
                let mut corners = BTreeSet::new();
                if left > before.0 {
                    corners.insert(left - 1);
                }
                if left < before.1 {
                    corners.insert(left);
                }
                if right > left {
                    corners.insert(right - 1);
                }
                if right < before.1 {
                    corners.insert(right);
                }
                assert!(!corners.is_empty());
                let graph_id = [p.id, loan_index, row_index];
                let mut witnesses = vec![];
                line(&mut gf, &json!({"graph_id":graph_id,"graph":g}));
                for index in corners {
                    let args = [
                        Val::f64(rate),
                        Val::f64(pv),
                        latent[index],
                        Val::f64(0.0),
                        Val::f64(0.0),
                    ];
                    let values = graph_values_typed::<F>(&g, &args);
                    let got = values[g.output].to_f64::<F>();
                    assert!(got.is_finite());
                    let trace_id = trace_count;
                    trace_count += 1;
                    line(
                        &mut tf,
                        &json!({"trace_id":trace_id,"graph_id":graph_id,"latent_index":index,
                        "inputs":args.map(Val::exact_encoding),"values":values.into_iter().map(Val::exact_encoding).collect::<Vec<_>>(),
                        "actual":format!("0x{:016x}",got.to_bits())}),
                    );
                    witnesses.push(json!({"latent_index":index,"trace_id":trace_id,"actual":format!("0x{:016x}",got.to_bits())}));
                }
                line(
                    &mut cf,
                    &json!({"profile_id":p.id,"loan_index":loan_index,"row_index":row_index,
                    "before":[before.0,before.1],"after":[low,high],"boundaries":witnesses,
                    "expected":rows[row_index]["expected"],"graph_id":graph_id}),
                );
                if low == high {
                    first_kill = Some((loan_index, row_index));
                    break;
                }
            }
            if first_kill.is_some() {
                break;
            }
            feasible_loans.push(json!({"loan_index":loan_index,"latent_interval":[low,high]}));
        }
        if let Some((loan, row)) = first_kill {
            killed += 1;
            *histogram.entry(row).or_default() += 1;
            assert!(loan < 60);
        } else {
            surviving += 1;
            assert_eq!(feasible_loans.len(), 60);
            line(
                &mut sf,
                &json!({"profile_id":p.id,"loans":feasible_loans,
                "meaning":"necessary conditional consumer; payment-producing program remains unknown"}),
            );
        }
        completed += 1;
        if completed % 250 == 0 {
            println!(
                "balance profiles {completed}/7500; survivors {surviving}; evaluations {evaluations}"
            );
        }
    }
    for f in [&mut pf, &mut cf, &mut gf, &mut tf, &mut sf] {
        f.flush().unwrap();
    }
    let summary = json!({"declared_profiles":7500,"completed_profiles":completed,"unsearched_profiles":7500-completed,
        "complete":completed==7500,"killed":killed,"survivors":surviving,"boundary_traces":trace_count,
        "bisection_evaluations":evaluations,"original_available_rows":540,"loan_groups":60,
        "first_kill_histogram":histogram,"seconds":clock.elapsed().as_secs_f64(),
        "initial_domain":if free_domains.is_some(){"complete first-singleton inverse; no external PMT constraint"}else{"external PMT publication"},
        "interest_sum":interest_sum,
        "claim":"Exact lattice constraints for positive-coefficient monotone balance/accumulation consumers; no payment fitted or installed. Surviving profiles are partial programs conditional on a common latent payment per loan."});
    crate::json_write(root.join("summary.json"), &summary);
    println!("{summary}");
}
