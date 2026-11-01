//! R11 bounded cumulative-principal programs, from ordinary annuity algebra.
//! No observed payment or principal is installed as a component or constant.
use crate::{dsl::*, eval::*, ported::excel_numeric::Fpu};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::{
    collections::{BTreeMap, HashMap},
    fs,
    io::{BufWriter, Write},
    path::Path,
};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Program {
    pub id: usize,
    pub power: usize,
    pub growth_model: EvalModel,
    pub payment_form: usize,
    pub payment_model: EvalModel,
    pub payment_store: bool,
    pub principal_form: usize,
    pub principal_model: EvalModel,
    pub sum_order: usize,
    pub sum_model: EvalModel,
}

pub fn programs() -> Vec<Program> {
    let mut result = Vec::new();
    for power in 0..2 {
        for growth_model in WIDE_MODELS {
            for payment_form in 0..3 {
                for payment_model in WIDE_MODELS {
                    for payment_store in [false, true] {
                        for principal_form in 0..4 {
                            for principal_model in WIDE_MODELS {
                                for sum_order in 0..3 {
                                    for sum_model in WIDE_MODELS {
                                        result.push(Program {
                                            id: result.len(),
                                            power,
                                            growth_model,
                                            payment_form,
                                            payment_model,
                                            payment_store,
                                            principal_form,
                                            principal_model,
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
        }
    }
    assert_eq!(result.len(), 90000);
    result
}

struct Builder<'a> {
    g: GraphBuilder,
    p: &'a Program,
    one: usize,
    minus_one: usize,
    zero: usize,
    q: usize,
    powers: HashMap<usize, usize>,
}
impl Builder<'_> {
    fn add(&mut self, a: usize, b: usize, m: EvalModel) -> usize {
        self.g.push(Op::Add(a, b), m)
    }
    fn mul(&mut self, a: usize, b: usize, m: EvalModel) -> usize {
        self.g.push(Op::Mul(a, b), m)
    }
    fn div(&mut self, a: usize, b: usize, m: EvalModel) -> usize {
        self.g.push(Op::Div(a, b), m)
    }
    fn power(&mut self, n: usize) -> usize {
        if let Some(&value) = self.powers.get(&n) {
            return value;
        }
        let m = self.p.growth_model;
        let value = if self.p.power == 0 {
            // A sequential product, reusing its preceding powers.
            let mut value = self.one;
            for k in 1..=n {
                value = if let Some(&old) = self.powers.get(&k) {
                    old
                } else {
                    let next = self.mul(value, self.q, m);
                    self.powers.insert(k, next);
                    next
                };
            }
            value
        } else {
            let mut k = n;
            let mut b = self.q;
            let mut value = self.one;
            while k > 0 {
                if k & 1 != 0 {
                    value = self.mul(value, b, m);
                }
                k >>= 1;
                if k > 0 {
                    b = self.mul(b, b, m);
                }
            }
            value
        };
        self.powers.insert(n, value);
        value
    }
}

pub fn graph(p: &Program, periods: usize, start: usize, end: usize, beginning: bool) -> Graph {
    assert!(1 <= start && start <= end && end <= periods && periods <= 512);
    let mut g = GraphBuilder::default();
    let rate = g.strict(Op::Arg(0));
    let pv = g.strict(Op::Arg(1));
    let negative_pv = g.strict(Op::Arg(2));
    let one = g.strict(Op::Const {
        bits: 1f64.to_bits(),
    });
    let minus_one = g.strict(Op::Const {
        bits: (-1f64).to_bits(),
    });
    let zero = g.strict(Op::Const { bits: 0 });
    let q = g.push(Op::Add(one, rate), p.growth_model);
    let mut b = Builder {
        g,
        p,
        one,
        minus_one,
        zero,
        q,
        powers: HashMap::from([(0, one)]),
    };
    let growth = b.power(periods);
    let pm = p.payment_model;
    let mut payment = match p.payment_form {
        0 => {
            let interest = b.mul(negative_pv, rate, pm);
            let numerator = b.mul(interest, growth, pm);
            let denominator = b.add(growth, minus_one, pm);
            b.div(numerator, denominator, pm)
        }
        1 | 2 => {
            let negative_discount = b.div(minus_one, growth, pm);
            let denominator = b.add(one, negative_discount, pm);
            if p.payment_form == 1 {
                let interest = b.mul(negative_pv, rate, pm);
                b.div(interest, denominator, pm)
            } else {
                let factor = b.div(rate, denominator, pm);
                b.mul(negative_pv, factor, pm)
            }
        }
        _ => unreachable!(),
    };
    if beginning {
        payment = b.div(payment, q, pm);
    }
    // Explicit binary64 store, represented by exact multiplication by one in
    // strict mode. No extra precision is introduced by a hidden helper.
    if p.payment_store {
        payment = b.mul(payment, one, EvalModel::Strict);
    }
    let m = p.principal_model;
    let mut terms = Vec::with_capacity(end - start + 1);
    let mut balance = pv;
    let mut previous = zero;
    for period in 1..=end {
        let principal = match p.principal_form {
            0 => {
                if beginning && period == 1 {
                    payment
                } else {
                    let discount = b.power(periods - period + 1);
                    b.div(payment, discount, m)
                }
            }
            1 => {
                if period == 1 {
                    if beginning {
                        payment
                    } else {
                        let interest = b.mul(pv, rate, m);
                        b.add(payment, interest, m)
                    }
                } else if beginning && period == 2 {
                    let adjusted = b.mul(payment, q, m);
                    let interest = b.mul(pv, rate, m);
                    b.add(adjusted, interest, m)
                } else {
                    b.mul(previous, q, m)
                }
            }
            2 | 3 => {
                if beginning && period == 1 {
                    payment
                } else {
                    let interest = b.mul(balance, rate, m);
                    b.add(payment, interest, m)
                }
            }
            _ => unreachable!(),
        };
        if period >= start {
            terms.push(principal);
        }
        if p.principal_form >= 2 && period < end {
            balance = if beginning && period == 1 || p.principal_form == 2 {
                b.add(balance, principal, m)
            } else {
                let grown = b.mul(balance, q, m);
                b.add(grown, payment, m)
            };
        }
        previous = principal;
    }
    let output = sum_graph(&mut b.g, &terms, p.sum_order, p.sum_model, zero);
    b.g.finish(output)
}

fn arguments(row: &Value) -> ([f64; 5], usize, usize, usize, bool) {
    assert_eq!(row["function"], "CUMPRINC");
    assert!(row["result_index"].is_null());
    let args = row["args"].as_array().unwrap();
    assert_eq!(args.len(), 6);
    let a: Vec<f64> = args.iter().map(crate::annuity::number).collect();
    assert!(
        a.iter()
            .all(|v| v.is_finite() && (*v == 0.0 || v.is_normal()))
    );
    assert!(a[0] > 0.0 && a[2] > 0.0 && matches!(a[5], 0.0 | 1.0));
    assert!(a[1..=1].iter().chain(&a[3..=4]).all(|v| v.fract() == 0.0));
    let (n, s, e) = (a[1] as usize, a[3] as usize, a[4] as usize);
    assert!(1 <= s && s <= e && e <= n && n <= 512);
    ([a[0], a[2], -a[2], 0.0, 0.0], n, s, e, a[5] != 0.0)
}

fn line(writer: &mut impl Write, value: &impl Serialize) {
    serde_json::to_writer(&mut *writer, value).unwrap();
    writeln!(writer).unwrap();
}

pub fn run<F: Fpu>(input: &str, out: &str) {
    let rows = crate::lines::<Value>(input);
    assert_eq!(rows.len(), 540);
    let cases: Vec<_> = rows.iter().map(arguments).collect();
    fs::create_dir_all(out).unwrap();
    let dir = Path::new(out);
    let mut pf = BufWriter::new(fs::File::create(dir.join("programs.jsonl")).unwrap());
    let mut kf = BufWriter::new(fs::File::create(dir.join("first-kills.jsonl")).unwrap());
    let mut tf = BufWriter::new(fs::File::create(dir.join("traces.jsonl")).unwrap());
    let mut survivors = Vec::new();
    let mut histogram = BTreeMap::<usize, usize>::new();
    // Reuse exactly equal instantiated graphs, preserving every logical program.
    let mut cache = HashMap::<Vec<u8>, (usize, String)>::new();
    let mut next_trace = 0;
    let mut evaluations = 0usize;
    for p in programs() {
        line(&mut pf, &p);
        let mut killed = false;
        for (i, (args, n, s, e, t)) in cases.iter().enumerate() {
            let g = graph(&p, *n, *s, *e, *t);
            let key = serde_json::to_vec(&json!([i, &g])).unwrap();
            let (trace_id, actual) = if let Some(found) = cache.get(&key) {
                found.clone()
            } else {
                let value = eval_graph::<F>(&g, args);
                let actual = outcome::<F>(value);
                let trace_id = next_trace;
                next_trace += 1;
                evaluations += 1;
                line(
                    &mut tf,
                    &json!({"trace_id":trace_id,"row_index":i,"graph":g,
                    "input_bits":args.map(f64::to_bits),"values":trace_graph::<F>(&g,args),"actual":actual}),
                );
                cache.insert(key, (trace_id, actual.clone()));
                (trace_id, actual)
            };
            if actual != rows[i]["expected"].as_str().unwrap() {
                line(
                    &mut kf,
                    &json!({"program_id":p.id,"row_index":i,"trace_id":trace_id,
                    "expected":rows[i]["expected"],"actual":actual,"raw":rows[i]}),
                );
                *histogram.entry(i).or_default() += 1;
                killed = true;
                break;
            }
        }
        if !killed {
            survivors.push(p.id);
        }
    }
    pf.flush().unwrap();
    kf.flush().unwrap();
    tf.flush().unwrap();
    crate::json_write(dir.join("survivors.json"), &survivors);
    let summary = json!({"programs":90000,"original_rows":540,"unique_instantiated_graphs":evaluations,
        "survivors":survivors.len(),"first_kill_histogram":histogram,
        "constants":[0,1,-1],"growth":"integer sequential or binary-squaring; all arithmetic explicit",
        "claim":"Exact finite-family elimination; no ranking, fitted constants, answer components or identified Excel program. Logical alternatives with inactive choices are preserved."});
    crate::json_write(dir.join("summary.json"), &summary);
    println!("{summary}");
}

/// Growth-only audit after the concrete complete-program family has failed.
/// These are explicitly calculated components, not observations of Excel internals.
pub fn growth_tables<F: Fpu>(input: &str, output: &str) {
    let rows = crate::lines::<Value>(input);
    let mut groups = BTreeMap::<String, usize>::new();
    for (i, row) in rows.iter().enumerate() {
        let _ = arguments(row);
        let key = serde_json::to_string(&json!([
            row["args"][0],
            row["args"][1],
            row["source"],
            row["reading"]
        ]))
        .unwrap();
        groups.entry(key).or_insert(i);
    }
    let mut writer = BufWriter::new(fs::File::create(output).unwrap());
    for (group, row_index) in groups {
        let (args, n, _, _, _) = arguments(&rows[row_index]);
        for power in 0..2 {
            for growth_model in WIDE_MODELS {
                let p = Program {
                    id: 0,
                    power,
                    growth_model,
                    payment_form: 0,
                    payment_model: EvalModel::Strict,
                    payment_store: false,
                    principal_form: 0,
                    principal_model: EvalModel::Strict,
                    sum_order: 0,
                    sum_model: EvalModel::Strict,
                };
                let mut g = GraphBuilder::default();
                let rate = g.strict(Op::Arg(0));
                let one = g.strict(Op::Const {
                    bits: 1f64.to_bits(),
                });
                let minus_one = g.strict(Op::Const {
                    bits: (-1f64).to_bits(),
                });
                let zero = g.strict(Op::Const { bits: 0 });
                let q = g.push(Op::Add(one, rate), growth_model);
                let mut b = Builder {
                    g,
                    p: &p,
                    one,
                    minus_one,
                    zero,
                    q,
                    powers: HashMap::from([(0, one)]),
                };
                let outputs: Vec<usize> = (0..=n).map(|k| b.power(k)).collect();
                let graph = b.g.finish(outputs[n]);
                let values = trace_graph::<F>(&graph, &args);
                let factors: Vec<&str> = outputs.iter().map(|&i| values[i].as_str()).collect();
                line(
                    &mut writer,
                    &json!({"group":group,"row_index":row_index,
                "power":power,"growth_model":growth_model,"periods":n,
                "input_bits":args.map(f64::to_bits),"outputs":outputs,
                "factors":factors,"graph":graph,"values":values,
                "meaning":"calculated conditional components, not Excel helper observations"}),
                );
            }
        }
    }
    writer.flush().unwrap();
}
