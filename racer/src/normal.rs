//! Bounded normal-CDF sketches from Cody's published rational tables.
//! CALERF: netlib SPECFUN/erf; direct tables: SciPy1.14.1 cdflib.c:cumnor.
//! See notices/SciPy-cdflib.txt. No application normal kernel/history is read.
use crate::{
    boundary::{neg, store},
    dsl::*,
    eval::*,
    json_write,
    ported::excel_numeric::*,
};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::{
    collections::{BTreeMap, HashMap},
    fs,
    io::{BufWriter, Write},
    path::Path,
};

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct Program {
    pub id: usize,
    pub sketch: usize, // 0: direct normal; 1: scaled CALERF; 2: published Hart body
    pub extended_constants: bool,
    pub normalization: usize, // divide sqrt2, multiply reciprocal, divide sqrt(2)
    pub reduction: EvalModel,
    pub polynomial: EvalModel,
    pub quotient: EvalModel,
    pub exponent: EvalModel,
    pub reconstruction: EvalModel,
    pub exponential: usize, // resident cleared, stored cleared, CRT scalar, libm scalar
    pub split: usize,       // direct square, split at1/16, split at1/1.6
    pub association: usize,
    #[serde(default)]
    pub reflection: usize, // four algebraic publication forms, in the published middle region
}

fn family(wrappers: bool, hart: bool) -> Vec<Program> {
    let mut out = vec![];
    for sketch in if hart { 2..3 } else { 0..2 } {
        for extended_constants in [false, true] {
            for normalization in 0..if sketch != 1 { 1 } else { 3 } {
                for reduction in if sketch != 1 {
                    vec![EvalModel::Strict]
                } else {
                    MODELS.to_vec()
                } {
                    for polynomial in WIDE_MODELS {
                        for quotient in [polynomial] {
                            for reconstruction in WIDE_MODELS {
                                for exponent in [reconstruction] {
                                    for exponential in 0..4 {
                                        for split in 0..3 {
                                            for association in 0..3 {
                                                for reflection in 0..if wrappers { 4 } else { 1 } {
                                                    out.push(Program {
                                                        id: out.len(),
                                                        sketch,
                                                        extended_constants,
                                                        normalization,
                                                        reduction,
                                                        polynomial,
                                                        quotient,
                                                        exponent,
                                                        reconstruction,
                                                        exponential,
                                                        split,
                                                        association,
                                                        reflection,
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
            }
        }
    }
    assert_eq!(
        out.len(),
        (if hart { 1800 } else { 18_000 }) * (if wrappers { 4 } else { 1 })
    );
    out
}

#[derive(Clone, Copy)]
struct V {
    value: Val,
    id: usize,
}
struct Machine<'a, F: Fpu> {
    nodes: Vec<Value>,
    next: usize,
    record: bool,
    constants: &'a Value,
    leaves: &'a mut BTreeMap<String, Value>,
    marker: std::marker::PhantomData<F>,
}
impl<'a, F: Fpu> Machine<'a, F> {
    fn push(&mut self, mut node: Value, value: Val) -> V {
        let id = self.next;
        self.next += 1;
        if self.record {
            node["value"] = json!(value.exact_encoding());
            self.nodes.push(node);
        }
        V { value, id }
    }
    fn literal(&mut self, x: f64) -> V {
        self.push(
            json!({"op":"literal", "bits":format!("{:016x}",x.to_bits())}),
            Val::f64(x),
        )
    }
    fn constant(&mut self, table: &str, index: usize, extended: bool) -> V {
        let c = &self.constants["tables"][table][index];
        let value = if extended {
            let s = c["ext80_le"].as_str().unwrap();
            let raw: Vec<_> = (0..20)
                .step_by(2)
                .map(|i| u8::from_str_radix(&s[i..i + 2], 16).unwrap())
                .collect();
            Val::Ext(Extended(raw.try_into().unwrap()))
        } else {
            Val::f64(f64::from_bits(
                u64::from_str_radix(c["binary64"].as_str().unwrap(), 16).unwrap(),
            ))
        };
        self.push(
            json!({"op":"constant","table":table,"index":index,"extended":extended}),
            value,
        )
    }
    fn op(&mut self, kind: &str, model: EvalModel, a: V, b: V) -> V {
        let operation = match kind {
            "add" => Op::Add(0, 1),
            "mul" => Op::Mul(0, 1),
            "div" => Op::Div(0, 1),
            _ => unreachable!(),
        };
        let value = binop::<F>(operation, model, a.value, b.value);
        self.push(json!({"op":kind,"model":model,"args":[a.id,b.id]}), value)
    }
    fn minus(&mut self, a: V) -> V {
        self.push(json!({"op":"neg","args":[a.id]}), neg(a.value))
    }
    fn sub(&mut self, model: EvalModel, a: V, b: V) -> V {
        let b = self.minus(b);
        self.op("add", model, a, b)
    }
    fn sqrt(&mut self, a: V, model: EvalModel) -> V {
        let cw = if matches!(model, EvalModel::StoredPc53 | EvalModel::ResidentPc53) {
            (CW & !0x0300) | 0x0200
        } else {
            CW
        };
        let value = if model == EvalModel::Strict {
            Val::f64(a.value.to_f64::<F>().sqrt())
        } else {
            store::<F>(
                Val::Ext(F::sqrt(&a.value.to_ext(), cw)),
                matches!(model, EvalModel::StoredPc64 | EvalModel::StoredPc53),
            )
        };
        self.push(json!({"op":"sqrt","model":model,"args":[a.id]}), value)
    }
    fn floor(&mut self, a: V) -> V {
        let value = Val::Ext(F::rndint(&a.value.to_ext(), (CW & !0x0c00) | 0x0400));
        self.push(json!({"op":"floor","args":[a.id]}), value)
    }
    fn exp(&mut self, a: V, provider: usize) -> V {
        let key = format!("{provider}:{}", a.value.exact_encoding());
        let value = match provider {
            0 => Val::Ext(research_exp::<F>(&a.value.to_ext())),
            1 => Val::f64(F::store(&research_exp::<F>(&a.value.to_ext()), CW)),
            2 => Val::f64(crt_exp(a.value.to_f64::<F>())),
            3 => Val::f64(libm::exp(a.value.to_f64::<F>())),
            _ => unreachable!(),
        };
        let leaf = json!({"provider":provider,"input":a.value.exact_encoding(),"output":value.exact_encoding()});
        if let Some(old) = self.leaves.insert(key.clone(), leaf.clone()) {
            assert_eq!(old, leaf, "helper is not deterministic");
        }
        self.push(
            json!({"op":"exp","args":[a.id],"provider":provider,"concrete_leaf":key}),
            value,
        )
    }
    fn horner(&mut self, table: &str, z: V, ext: bool, model: EvalModel) -> V {
        let n = self.constants["tables"][table].as_array().unwrap().len();
        let mut v = self.constant(table, n - 1, ext);
        for i in 0..n - 1 {
            v = self.op("mul", model, v, z);
            let c = self.constant(table, i, ext);
            v = self.op("add", model, v, c);
        }
        v
    }
    fn denominator(&mut self, table: &str, z: V, ext: bool, model: EvalModel) -> V {
        let n = self.constants["tables"][table].as_array().unwrap().len();
        let mut v = z;
        for i in 0..n {
            let c = self.constant(table, i, ext);
            v = self.op("add", model, v, c);
            if i + 1 < n {
                v = self.op("mul", model, v, z);
            }
        }
        v
    }
    fn ordinary_polynomial(&mut self, table: &str, z: V, ext: bool, model: EvalModel) -> V {
        let n = self.constants["tables"][table].as_array().unwrap().len();
        let mut v = self.constant(table, n - 1, ext);
        for i in (0..n - 1).rev() {
            v = self.op("mul", model, v, z);
            let c = self.constant(table, i, ext);
            v = self.op("add", model, v, c);
        }
        v
    }
    fn triple(&mut self, a: V, b: V, c: V, model: EvalModel, association: usize) -> V {
        match association {
            0 => {
                let ab = self.op("mul", model, a, b);
                self.op("mul", model, ab, c)
            }
            1 => {
                let bc = self.op("mul", model, b, c);
                self.op("mul", model, a, bc)
            }
            _ => {
                let ab = self.op("mul", model, a, b);
                self.op("mul", model, c, ab)
            }
        }
    }
}

fn evaluate<F: Fpu>(
    p: &Program,
    x: f64,
    constants: &Value,
    leaves: &mut BTreeMap<String, Value>,
    record: bool,
) -> (String, Vec<Value>, String) {
    let mut m = Machine::<F> {
        nodes: vec![],
        next: 0,
        record,
        constants,
        leaves,
        marker: std::marker::PhantomData,
    };
    let _xnode = m.push(
        json!({"op":"input","bits":format!("{:016x}",x.to_bits())}),
        Val::f64(x),
    );
    let half = m.literal(0.5);
    let one = m.literal(1.0);
    let positive = m.literal(x.abs());
    let y = if p.sketch != 1 {
        positive
    } else {
        let c = if p.normalization == 2 {
            let two = m.literal(2.0);
            m.sqrt(two, p.reduction)
        } else {
            m.constant(
                if p.normalization == 0 {
                    "sqrt2"
                } else {
                    "invsqrt2"
                },
                0,
                p.extended_constants,
            )
        };
        m.op(
            if p.normalization == 1 { "mul" } else { "div" },
            p.reduction,
            positive,
            c,
        )
    };
    // Comparisons use full extended operands, not their binary64 publication.
    let yf = y.value.to_ext();
    let le = |bound: f64| {
        let b = research_load(bound);
        let d = F::sub(&yf, &b, CW);
        d.0[9] & 0x80 != 0 || d.0[..8].iter().all(|b| *b == 0)
    };
    let central = p.sketch != 2 && le(if p.sketch == 0 { 0.66291 } else { 0.46875 });
    let middle = p.sketch == 2
        || le(if p.sketch == 0 {
            std::f64::consts::SQRT_2 * 4.0
        } else {
            4.0
        });
    let region = if central {
        "central"
    } else if middle {
        "middle"
    } else {
        "tail"
    };
    let prefix = if p.sketch == 0 {
        "direct"
    } else if p.sketch == 1 {
        "erfc"
    } else {
        "hart"
    };
    let table = |suffix: &str| format!("{prefix}_{suffix}");
    let mut answer = if central {
        let square = m.op("mul", p.polynomial, y, y);
        let n = m.horner(&table("a"), square, p.extended_constants, p.polynomial);
        let d = m.denominator(&table("b"), square, p.extended_constants, p.polynomial);
        let delta = match p.association {
            0 => {
                let n = m.op("mul", p.quotient, y, n);
                m.op("div", p.quotient, n, d)
            }
            1 => {
                let r = m.op("div", p.quotient, n, d);
                m.op("mul", p.quotient, y, r)
            }
            _ => {
                let r = m.op("div", p.quotient, y, d);
                m.op("mul", p.quotient, r, n)
            }
        };
        if p.sketch == 0 {
            if x < 0.0 {
                m.sub(p.reconstruction, half, delta)
            } else {
                m.op("add", p.reconstruction, half, delta)
            }
        } else {
            // Keep direct half-plus-erf and the erfc wrapper as explicit sketches.
            let erfc = if x < 0.0 {
                m.sub(p.reconstruction, one, delta)
            } else {
                m.op("add", p.reconstruction, one, delta)
            };
            m.op("mul", p.reconstruction, half, erfc)
        }
    } else if x.abs() > 40.0 {
        m.literal(if x < 0.0 { 0.0 } else { 1.0 })
    } else {
        let ratio = if p.sketch == 2 {
            let n = m.ordinary_polynomial("hart_p", y, p.extended_constants, p.polynomial);
            let d = m.ordinary_polynomial("hart_q", y, p.extended_constants, p.polynomial);
            m.op("div", p.quotient, n, d)
        } else if middle {
            let n = m.horner(&table("c"), y, p.extended_constants, p.polynomial);
            let d = m.denominator(&table("d"), y, p.extended_constants, p.polynomial);
            m.op("div", p.quotient, n, d)
        } else {
            let square = m.op("mul", p.polynomial, y, y);
            let z = m.op("div", p.polynomial, one, square);
            let n = m.horner(&table("p"), z, p.extended_constants, p.polynomial);
            let d = m.denominator(&table("q"), z, p.extended_constants, p.polynomial);
            let n = m.op("mul", p.quotient, z, n);
            let r = m.op("div", p.quotient, n, d);
            let c = m.constant(&table("sqrpi"), 0, p.extended_constants);
            let r = m.sub(p.quotient, c, r);
            m.op("div", p.quotient, r, y)
        };
        let scale = if p.sketch != 1 { half } else { one };
        let (a, b) = if p.split == 0 {
            let sq = m.op("mul", p.exponent, y, y);
            let sq = m.op("mul", p.exponent, sq, scale);
            (m.minus(sq), m.literal(0.0))
        } else {
            let grid = m.literal(if p.split == 1 { 16.0 } else { 1.6 });
            let t = m.op("mul", p.exponent, y, grid);
            let t = m.floor(t);
            let t = m.op("div", p.exponent, t, grid);
            let lo = m.sub(p.exponent, y, t);
            let hi = m.op("add", p.exponent, y, t);
            let delta = m.op("mul", p.exponent, lo, hi);
            let delta = m.op("mul", p.exponent, delta, scale);
            let sq = m.op("mul", p.exponent, t, t);
            let sq = m.op("mul", p.exponent, sq, scale);
            (m.minus(sq), m.minus(delta))
        };
        let ea = m.exp(a, p.exponential);
        let eb = m.exp(b, p.exponential);
        let mut small = m.triple(ea, eb, ratio, p.reconstruction, p.association);
        if p.sketch == 1 {
            small = m.op("mul", p.reconstruction, small, half);
        }
        if x > 0.0 {
            m.sub(p.reconstruction, one, small)
        } else {
            small
        }
    };
    if region == "middle" {
        answer = match p.reflection {
            0 => answer,
            1 => {
                let other = m.sub(p.reconstruction, one, answer);
                m.sub(p.reconstruction, one, other)
            }
            2 => {
                let other = m.sub(p.reconstruction, half, answer);
                m.sub(p.reconstruction, half, other)
            }
            3 => {
                let two = m.literal(2.0);
                let doubled = m.op("mul", p.reconstruction, answer, two);
                let erf = m.sub(p.reconstruction, one, doubled);
                let erfc = m.sub(p.reconstruction, one, erf);
                m.op("mul", p.reconstruction, half, erfc)
            }
            _ => unreachable!(),
        };
    }
    let v = answer.value.to_f64::<F>();
    // CDFLIB's published underflow rule: values below the least normal are zero.
    let encoded = if v == 0.0 || v.is_subnormal() {
        "0x0000000000000000".into()
    } else if !v.is_finite() {
        "error:Num".into()
    } else {
        format!("0x{:016x}", v.to_bits())
    };
    (encoded, m.nodes, region.into())
}

pub fn run<F: Fpu>(
    trace: &str,
    constant_path: &str,
    output: &str,
    build: &str,
    wrappers: bool,
    hart: bool,
    central_only: bool,
) {
    let constants: Value = serde_json::from_slice(&fs::read(constant_path).unwrap()).unwrap();
    let out = Path::new(output);
    fs::create_dir_all(out).unwrap();
    let programs = family(wrappers, hart);
    let mut file = BufWriter::new(fs::File::create(out.join("programs.jsonl")).unwrap());
    for p in &programs {
        writeln!(file, "{}", serde_json::to_string(p).unwrap()).unwrap();
    }
    file.flush().unwrap();
    let mut scoped = vec![];
    let mut counts = BTreeMap::<String, usize>::new();
    let mut excluded = BufWriter::new(fs::File::create(out.join("excluded.jsonl")).unwrap());
    for (i, raw) in crate::lines::<Value>(trace).into_iter().enumerate() {
        let args = raw["args"].as_array().unwrap();
        let a: Result<Vec<_>, _> = args
            .iter()
            .map(|a| calc_verdict::decode_argument(a, "arg"))
            .collect();
        let mut x = None;
        let reason = if raw["reading"]["excel_build"].as_str() != Some(build) {
            "other_environment"
        } else if !raw["result_index"].is_null() || args.len() != 2 {
            "unsupported_shape"
        } else if let Ok(a) = a {
            let cdf = match a[1].core() {
                crate::value::CoreValue::Logical(b) => Some(*b),
                crate::value::CoreValue::Number(n) => Some(*n != 0.0),
                _ => None,
            };
            if cdf != Some(true) {
                "density_or_flag_control"
            } else if let Some(z) = a[0]
                .as_number()
                .filter(|n| n.is_finite() && !n.is_subnormal())
            {
                if central_only && z.abs() > 0.66291 {
                    "outside_published_central_region"
                } else if hart && z.abs() >= 7.07106781186547 {
                    "outside_hart_published_body"
                } else {
                    x = Some(z);
                    "scoped"
                }
            } else {
                "argument_admission_control"
            }
        } else {
            "failed_decode"
        };
        *counts.entry(reason.into()).or_default() += 1;
        if let Some(x) = x {
            calc_verdict::decode_result(raw["expected"].as_str().unwrap()).unwrap();
            scoped.push((i, raw, x));
        } else {
            writeln!(
                excluded,
                "{}",
                json!({"trace_index":i,"reason":reason,"raw":raw})
            )
            .unwrap();
        }
    }
    // All rows are discovery. Scheduling is independent of residual magnitude.
    scoped.sort_by_key(|(i, _, x)| {
        (
            if x.abs() > 0.66291 && x.abs() < 5.0 {
                0
            } else if x.abs() < 0.66291 {
                1
            } else {
                2
            },
            *i,
        )
    });
    json_write(
        out.join("scope.json"),
        &json!({"counts":counts,"rows":scoped.len(),"build":build,"discovery_only":true,"central_only":central_only,"hart_body_only":hart}),
    );
    let mut alive: Vec<_> = (0..programs.len()).collect();
    let mut leaves = BTreeMap::new();
    let mut kills = BufWriter::new(fs::File::create(out.join("first-kills.jsonl")).unwrap());
    let mut total = 0;
    let mut evaluated = 0;
    let mut killed = 0;
    let mut stats = vec![];
    for (order, (trace_index, raw, x)) in scoped.iter().enumerate() {
        if alive.is_empty() {
            break;
        }
        let mut next = vec![];
        let mut groups: HashMap<String, Vec<usize>> = HashMap::new();
        // Retain physical programs even if some blocks are unused on this row.
        let mut cache = HashMap::<String, (String, String)>::new();
        for id in alive {
            let mut key = programs[id].clone();
            key.id = 0;
            // For central direct rational, exponential axes have no executable use.
            if key.sketch == 0 && x.abs() <= 0.66291 {
                key.exponent = EvalModel::Strict;
                key.exponential = 0;
                key.split = 0;
            }
            let signature = serde_json::to_string(&key).unwrap();
            let result = cache.entry(signature.clone()).or_insert_with(|| {
                evaluated += 1;
                let (got, _, region) =
                    evaluate::<F>(&programs[id], *x, &constants, &mut leaves, false);
                (got, region)
            });
            total += 1;
            if result.0 == raw["expected"].as_str().unwrap() { // ubs:ignore — public numerical observation, no secret comparison
                next.push(id);
            } else {
                groups.entry(signature).or_default().push(id);
                killed += 1;
            }
        }
        let mut group_keys: Vec<_> = groups.keys().cloned().collect();
        group_keys.sort();
        for key in group_keys {
            let ids = &groups[&key];
            let (got, nodes, region) =
                evaluate::<F>(&programs[ids[0]], *x, &constants, &mut leaves, true);
            assert_eq!((&got, &region), (&cache[&key].0, &cache[&key].1));
            writeln!(kills,"{}",json!({"candidate_ids":ids,"program":programs[ids[0]],"trace_index":trace_index,"order":order,
                "witness":raw,"got":got,"expected":raw["expected"],"nodes":nodes,"region":region})).unwrap();
        }
        alive = next;
        stats.push(json!({"order":order,"trace_index":trace_index,"id":raw["id"],"remaining":alive.len(),"physical_evaluations":total,"unique_evaluations":evaluated}));
        kills.flush().unwrap();
        json_write(
            out.join("checkpoint.json"),
            &json!({"order":order+1,"alive":alive,"complete":alive.is_empty()||order+1==scoped.len()}),
        );
        println!("{}", stats.last().unwrap());
    }
    json_write(
        out.join("survivors.json"),
        &alive.iter().map(|id| &programs[*id]).collect::<Vec<_>>(),
    );
    json_write(out.join("concrete-helpers.json"), &leaves);
    json_write(
        out.join("summary.json"),
        &json!({"programs":programs.len(),"survivors":alive.len(),"killed":killed,
        "scoped_rows":scoped.len(),"rows_executed":stats.len(),"evaluations":total,"unique_evaluations":evaluated,
        "build":build,"complete":true,"statistics":stats,"constant_file":constant_path,
        "conditional_helpers":true,"new_excel_answers":0,"product_changed":false}),
    );
}
