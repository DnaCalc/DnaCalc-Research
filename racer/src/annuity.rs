//! R8: named growth diagnostics before joint program search.
//! Own annuity algebra; cleared Fpu operations; no fitted constants or branches.
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
    fs,
    io::{BufWriter, Write},
    path::Path,
};

fn add<F: Fpu>(a: Val, b: Val, m: EvalModel) -> Val {
    binop::<F>(Op::Add(0, 1), m, a, b)
}
fn mul<F: Fpu>(a: Val, b: Val, m: EvalModel) -> Val {
    binop::<F>(Op::Mul(0, 1), m, a, b)
}
fn div<F: Fpu>(a: Val, b: Val, m: EvalModel) -> Val {
    binop::<F>(Op::Div(0, 1), m, a, b)
}
fn one() -> Val {
    Val::f64(1.0)
}
pub(crate) fn number(v: &Value) -> f64 {
    let d = calc_verdict::decode_argument(v, "annuity diagnostic arg").unwrap();
    match d.core() {
        crate::value::CoreValue::Number(x) => *x,
        _ => panic!("numeric diagnostic input"),
    }
}

pub fn growth_diagnostics<F: Fpu>(input: &str, output: &str) {
    let rows = crate::lines::<Value>(input);
    let mut f = BufWriter::new(fs::File::create(output).unwrap());
    for (idx, case) in rows.iter().enumerate() {
        let raw = &case["witness"];
        let r = number(&raw["args"][0]);
        let n = number(&raw["args"][1]);
        let mut candidates = Vec::new();
        let mut push = |name: String, v: Val| {
            candidates
                .push(json!({"name":name,"value":v.exact_encoding(),"publication":outcome::<F>(v)}))
        };
        if r > -1.0 && r.is_finite() && n.is_finite() {
            push(
                "libm_pow_stored_base".into(),
                Val::f64(libm::pow(1.0 + r, n)),
            );
            push(
                "cleared_crt_pow_stored_base".into(),
                Val::f64(crt_pow(1.0 + r, n)),
            );
            for m in MODELS {
                let base = add::<F>(one(), Val::f64(r), m);
                let power = Val::Ext(research_pow2::<F>(&F::fyl2x(
                    &research_load(n),
                    &base.to_ext(),
                    CW,
                )));
                push(format!("x87_pow_base_{m:?}"), power);
                push(
                    format!("x87_pow_base_{m:?}_stored"),
                    store::<F>(power, true),
                );
                let log = Val::f64(excel_log1p(r));
                let t = mul::<F>(Val::f64(n), log, m);
                let exp = Val::Ext(research_exp::<F>(&t.to_ext()));
                push(format!("exp_n_log1p_product_{m:?}"), exp);
                push(
                    format!("exp_n_log1p_product_{m:?}_stored"),
                    store::<F>(exp, true),
                );
                if n.fract() == 0.0 && n.abs() <= 8192.0 {
                    let mut repeated = one();
                    for _ in 0..n.abs() as usize {
                        repeated = mul::<F>(repeated, base, m);
                    }
                    if n < 0.0 {
                        repeated = div::<F>(one(), repeated, m);
                    }
                    push(format!("repeated_product_{m:?}"), repeated);
                    let mut k = n.abs() as u32;
                    let mut b = base;
                    let mut acc = one();
                    while k > 0 {
                        if k & 1 != 0 {
                            acc = mul::<F>(acc, b, m)
                        }
                        k >>= 1;
                        if k > 0 {
                            b = mul::<F>(b, b, m)
                        }
                    }
                    if n < 0.0 {
                        acc = div::<F>(one(), acc, m)
                    }
                    push(format!("binary_squaring_{m:?}"), acc);
                }
                // Mathematical convergence region fixed before observation.
                // These are explicit finite Taylor programs, not exact sums.
                if r.abs() <= 1.0 / 256.0 && (n * r).abs() <= 0.5 {
                    let mut term = one();
                    let mut sum = one();
                    for k in 1..=16 {
                        term = div::<F>(
                            mul::<F>(
                                mul::<F>(term, Val::f64(n - (k - 1) as f64), m),
                                Val::f64(r),
                                m,
                            ),
                            Val::f64(k as f64),
                            m,
                        );
                        sum = add::<F>(sum, term, m);
                        if [2, 4, 8, 16].contains(&k) {
                            push(format!("binomial_degree_{k}_{m:?}"), sum)
                        }
                    }
                }
            }
            let t = n * excel_log1p(r);
            push("libm_exp_n_log1p".into(), Val::f64(libm::exp(t)));
            push("cleared_exp_n_log1p".into(), Val::f64(exp_portable(t)));
        }
        serde_json::to_writer(&mut f,&json!({"case":idx,"rate":format!("0x{:016x}",r.to_bits()),"nper":format!("0x{:016x}",n.to_bits()),"candidates":candidates})).unwrap();
        writeln!(f).unwrap();
    }
    println!(
        "growth diagnostics: {} original isolated observations; named providers with explicit applicability",
        rows.len()
    );
}

#[derive(Clone, Serialize, Deserialize)]
pub struct AnnuityCase {
    pub trace_index: usize,
    pub raw: Value,
    pub args: [u64; 5],
}
impl AnnuityCase {
    pub fn numbers(&self) -> [f64; 5] {
        self.args.map(f64::from_bits)
    }
}

pub fn audit_decimal_inputs(input: &str, out: &str) {
    let mut records = vec![];
    let mut changed = 0;
    for (idx, row) in crate::lines::<Value>(input).iter().enumerate() {
        let mut differences = vec![];
        for j in 0..5 {
            let actual = row["args"][j].as_f64().unwrap();
            let requested = number(&row["raw"]["args"][j]);
            if actual.to_bits() != requested.to_bits() {
                changed += 1;
                differences.push(json!({"arg":j,"requested":format!("0x{:016x}",requested.to_bits()),"reloaded":format!("0x{:016x}",actual.to_bits())}));
            }
        }
        if !differences.is_empty() {
            records.push(json!({"case":idx,"function":row["raw"]["function"],"id":row["raw"]["id"],"differences":differences}));
        }
    }
    let result = json!({"rows":crate::lines::<Value>(input).len(),"changed_rows":records.len(),"changed_arguments":changed,"records":records,"claim":"These decimal-normalized research inputs fail exact original-bit identity; prior searches over them cannot exclude programmes on original Excel inputs."});
    json_write(out, &result);
    println!(
        "{}",
        json!({"changed_rows":records.len(),"changed_arguments":changed})
    );
}

pub fn transport_collisions(input: &str, old_decimal: &str, out: &str) {
    use std::io::BufRead;
    let transform =
        |x: f64| serde_json::from_str::<f64>(&serde_json::to_string(&x).unwrap()).unwrap();
    let mut controls = 0usize;
    let mut changed_controls = 0usize;
    for row in crate::lines::<Value>(old_decimal) {
        for j in 0..5 {
            let original = number(&row["raw"]["args"][j]);
            let t = transform(original);
            let observed = row["args"][j].as_f64().unwrap();
            assert_eq!(
                t.to_bits(),
                observed.to_bits(),
                "transport control matches the actual earlier normalization"
            );
            controls += 1;
            changed_controls += (t.to_bits() != original.to_bits()) as usize;
        }
    }
    let mut seen = [
        std::collections::HashMap::<String, (String, usize)>::new(),
        std::collections::HashMap::<String, (String, usize)>::new(),
    ];
    let mut stats =
        std::collections::BTreeMap::<String, std::collections::BTreeMap<String, usize>>::new();
    let mut conflicts = vec![];
    for (idx, line) in std::io::BufReader::new(fs::File::open(input).unwrap())
        .lines()
        .enumerate()
    {
        let row: Value = serde_json::from_str(&line.unwrap()).unwrap();
        let source = row["source"].as_str().unwrap();
        let c = stats.entry(source.into()).or_default();
        *c.entry("rows".into()).or_default() += 1;
        let mut args = row["args"].clone();
        let mut changed = 0usize;
        let mut unsupported = false;
        for arg in args.as_array_mut().unwrap() {
            if let Some(text) = arg
                .as_str()
                .filter(|s| s.starts_with("0x") && s.len() == 18)
            {
                let value = f64::from_bits(u64::from_str_radix(&text[2..], 16).unwrap());
                if !value.is_finite() {
                    unsupported = true;
                    break;
                }
                let result = transform(value);
                changed += (result.to_bits() != value.to_bits()) as usize;
                *arg = json!(format!("0x{:016x}", result.to_bits()));
            }
        }
        if unsupported {
            *c.entry("nonfinite_transport_control".into()).or_default() += 1;
            continue;
        }
        *c.entry("changed_arguments".into()).or_default() += changed;
        *c.entry("changed_rows".into()).or_default() += (changed != 0) as usize;
        let expected = row["expected"].as_str().unwrap();
        for mode in 0..2 {
            let encoded = if mode == 0 { &row["args"] } else { &args };
            let key = format!(
                "{}|{}|{}|{}|{}",
                row["function"], source, row["reading"], encoded, row["result_index"]
            );
            if let Some((old, first)) = seen[mode].get(&key) {
                *c.entry(format!("mode{mode}_merged_inputs")).or_default() += 1;
                if old != expected {
                    *c.entry(format!("mode{mode}_conflicting_observations"))
                        .or_default() += 1;
                    if conflicts.len() < 40 {
                        conflicts.push(json!({"mode":mode,"first_trace_index":first,"new_trace_index":idx,"first_expected":old,"new_expected":expected,"new_original_witness":row,"candidate_actual_input_args":encoded}));
                    }
                }
            } else {
                seen[mode].insert(key, (expected.into(), idx));
            }
        }
    }
    let result = json!({"transport_control_fields":controls,"changed_control_fields":changed_controls,"models":["identity","serde_json shortest-decimal serialization then default f64deserialization"],"source_statistics":stats,"conflict_examples":conflicts,
        "policy":"These are measurement-pipeline hypotheses, not modified Excel answers or a proposed worksheet kernel. Input transformation is uniform across every numeric argument and never chosen by disagreement. A same-source/context input collision with different outcomes contradicts a deterministic composition under that transport model.","limitation":"The legacy capture implementation is not inferred from this diagnostic; no source or build correction is applied to the evidence."});
    json_write(out, &result);
    println!(
        "{}",
        json!({"controls":controls,"changed_controls":changed_controls,"source_statistics":stats,"conflict_examples":conflicts.len()})
    );
}

pub fn prepare(pv: &str, fv: &str, out: &str, all_periods: bool) {
    fs::create_dir_all(out).unwrap();
    let mut rows = vec![];
    let mut excluded =
        BufWriter::new(fs::File::create(Path::new(out).join("excluded.jsonl")).unwrap());
    let mut counts = std::collections::BTreeMap::<String, usize>::new();
    for path in [fv, pv] {
        for (idx, raw) in crate::lines::<Value>(path).into_iter().enumerate() {
            let decoded: Option<Vec<f64>> = raw["args"].as_array().and_then(|a| {
                a.iter()
                    .map(|v| {
                        let d = calc_verdict::decode_argument(v, "annuity input").ok()?;
                        match d.core() {
                            crate::value::CoreValue::Number(x) => Some(*x),
                            _ => None,
                        }
                    })
                    .collect()
            });
            let args = decoded.and_then(|a| <[f64; 5]>::try_from(a).ok());
            let reason = match args {
                None => "typed_or_shape_control",
                Some(_) if !raw["result_index"].is_null() => "selected_element_control",
                Some(a) if a.iter().any(|x| !x.is_finite() || x.is_subnormal()) => {
                    "input_admission_control"
                }
                Some(a) if a[0] == 0.0 => "zero_rate_control",
                Some(a) if a[0] <= -1.0 => "nonpositive_base_control",
                Some(a)
                    if !all_periods
                        && (a[1] < 0.0 || a[1].fract() != 0.0 || a[1] >= u32::MAX as f64) =>
                {
                    "fractional_or_negative_period_control"
                }
                Some(_) => {
                    if all_periods {
                        "all_period_scoped"
                    } else {
                        "integer_scoped"
                    }
                }
            };
            *counts
                .entry(format!("{}:{reason}", raw["function"].as_str().unwrap()))
                .or_default() += 1;
            if reason == "integer_scoped" || reason == "all_period_scoped" {
                calc_verdict::decode_result(raw["expected"].as_str().unwrap()).unwrap();
                rows.push(AnnuityCase {
                    trace_index: idx,
                    raw,
                    args: args.unwrap().map(f64::to_bits),
                });
            } else {
                serde_json::to_writer(
                    &mut excluded,
                    &json!({"trace_index":idx,"reason":reason,"raw":raw}),
                )
                .unwrap();
                writeln!(excluded).unwrap();
            }
        }
    }
    // Start with both functions on shared rate/period/environment groups.
    // Selection uses input identities alone; all are inspected discovery data.
    let mut functions = std::collections::HashMap::<String, u8>::new();
    let group = |r: &AnnuityCase| {
        format!(
            "{}|{}|{}",
            r.raw["args"][0], r.raw["args"][1], r.raw["reading"]
        )
    };
    for r in &rows {
        *functions.entry(group(r)).or_default() |= if r.raw["function"] == "FV" { 1 } else { 2 };
    }
    rows.sort_by_key(|r| {
        (
            functions[&group(r)] != 3,
            r.args[0],
            r.args[1],
            r.raw["function"] == "FV",
            r.trace_index,
        )
    });
    let mut f = BufWriter::new(fs::File::create(Path::new(out).join("scoped.jsonl")).unwrap());
    for r in &rows {
        serde_json::to_writer(&mut f, r).unwrap();
        writeln!(f).unwrap()
    }
    json_write(
        Path::new(out).join("coverage.json"),
        &json!({"counts":counts,"scoped":rows.len(),"reader":"Original xtask witness traces; calc-verdict decode_argument/decode_result; no answer synthesis"}),
    );
    println!("{}", json!({"counts":counts,"scoped":rows.len()}));
}

fn fractional_growth<F: Fpu>(h: &Growth, r: f64, n: f64, provider: usize) -> Val {
    if n >= 0.0 && n.fract() == 0.0 && n < u32::MAX as f64 {
        return integer_growth::<F>(h, r, n);
    }
    if n < 0.0 && n.fract() == 0.0 && n.abs() < u32::MAX as f64 {
        return div::<F>(
            one(),
            store::<F>(integer_growth::<F>(h, r, -n), true),
            EvalModel::StoredPc64,
        );
    }
    let b = add::<F>(one(), Val::f64(r), h.base);
    let v = match provider {
        0..=7 => {
            let log_stored = provider & 1 == 0;
            let reciprocal_stage = provider & 2 == 0;
            let sqrt_half = provider & 4 == 0;
            let exponent = if reciprocal_stage { n.abs() } else { n };
            let positive = if sqrt_half && exponent == 0.5 {
                Val::f64(b.to_f64::<F>().sqrt())
            } else {
                let log = store::<F>(Val::Ext(research_log::<F>(&b.to_ext())), log_stored);
                let t = mul::<F>(Val::f64(exponent), log, h.product);
                Val::Ext(research_exp::<F>(&t.to_ext()))
            };
            if reciprocal_stage && n < 0.0 {
                div::<F>(one(), store::<F>(positive, true), EvalModel::StoredPc64)
            } else {
                positive
            }
        }
        8 => Val::f64(crt_pow(b.to_f64::<F>(), n)),
        9 => Val::f64(libm::pow(b.to_f64::<F>(), n)),
        10 => Val::f64(exp_portable(n * excel_log1p(r))),
        11 => {
            let t = mul::<F>(Val::f64(n), Val::f64(excel_log1p(r)), h.product);
            Val::Ext(research_exp::<F>(&t.to_ext()))
        }
        12 => Val::Ext(research_pow2::<F>(&F::fyl2x(
            &research_load(n),
            &b.to_ext(),
            CW,
        ))),
        _ => unreachable!(),
    };
    store::<F>(v, h.store)
}

pub(crate) fn loan_growth<F: Fpu>(growth_id: usize, rate: f64, periods: f64) -> Val {
    assert!(periods < 0.0 && periods.fract() == 0.0 && periods.abs() < u32::MAX as f64);
    fractional_growth::<F>(&growths()[growth_id], rate, periods, 0)
}

pub fn refine_fractional<F: Fpu>(input: &str, prior: &str, out: &str) {
    fs::create_dir_all(out).unwrap();
    let rows = crate::lines::<AnnuityCase>(input);
    let gs = growths();
    let prior: Vec<usize> = serde_json::from_slice(&fs::read(prior).unwrap()).unwrap();
    let mut alive: Vec<(usize, usize)> = prior
        .iter()
        .flat_map(|&p| (0..13).map(move |h| (p, h)))
        .collect();
    let count = alive.len();
    let mut evaluations = 0u64;
    let providers:Vec<_>=(0..13).map(|p|if p<8{format!("rounded_base_x87_log_exp: log_store={}, negative_reciprocal_stage={}, sqrt_half={}",p&1==0,p&2==0,p&4==0)}else{["cleared_crt_pow","libm_pow","portable_exp_n_log1p","x87_exp_n_log1p","x87_direct_pow"][p-8].into()}).collect();
    json_write(
        Path::new(out).join("grammar.json"),
        &json!({"integer_programs":prior,"fractional_providers":providers,
        "layout":"id=integer_program_id*13+provider","programs":count,
        "integer_growth":"Original nonnegative integer binary squaring. Negative integers use a stored-PC64 reciprocal of the stored positive power.",
        "fractional_growth":"Same base/product/store models as the integer helper; direct or abs-exponent reciprocal log/exp alternatives, sqrt-half predicate, or cleared scalar providers. No learned threshold.",
        "roots":"Original joint PV/FV integer graph and all its fixed site models; no per-row root selection",
        "input_policy":"All raw argument bits verified on load; positive-base nonzero-rate scope, zero/nonpositive-base controls retained separately"}),
    );
    let mut kills =
        BufWriter::new(fs::File::create(Path::new(out).join("first-kills.jsonl")).unwrap());
    let mut stats =
        BufWriter::new(fs::File::create(Path::new(out).join("witness-kills.jsonl")).unwrap());
    let mut killing = vec![];
    for (idx, row) in rows.iter().enumerate() {
        for j in 0..5 {
            assert_eq!(row.args[j], number(&row.raw["args"][j]).to_bits())
        }
        let before = alive.len();
        let mut rejected = std::collections::BTreeMap::<String, Vec<usize>>::new();
        if before > 0 {
            let a = row.numbers();
            let growth: Vec<Vec<_>> = gs
                .iter()
                .map(|g| {
                    (0..13)
                        .map(|p| fractional_growth::<F>(g, a[0], a[1], p))
                        .collect()
                })
                .collect();
            alive.retain(|&(id, p)| {
                let gi = id / (3 * PROFILES);
                let g = growth[gi][p];
                let got = if row.raw["function"] == "PV" && g.to_ext().0 == research_load(0.0).0 {
                    "error:Div0".into()
                } else {
                    outcome::<F>(evaluate::<F>(
                        row,
                        g,
                        (id / PROFILES) % 3,
                        models(id % PROFILES),
                    ))
                };
                evaluations += 1;
                if got == row.raw["expected"].as_str().unwrap() {
                    true
                } else {
                    rejected.entry(got).or_default().push(id * 13 + p);
                    false
                }
            });
            for (got, ids) in rejected {
                serde_json::to_writer(&mut kills,&json!({"candidate_ids":ids,"case":idx,"witness":row,"expected":row.raw["expected"],"got":got})).unwrap();
                writeln!(kills).unwrap();
            }
        }
        let rec = json!({"case":idx,"function":row.raw["function"],"id":row.raw["id"],"killed":before-alive.len(),"remaining":alive.len()});
        serde_json::to_writer(&mut stats, &rec).unwrap();
        writeln!(stats).unwrap();
        if before != alive.len() {
            println!("{rec}");
            killing.push(rec);
            kills.flush().unwrap();
            stats.flush().unwrap()
        }
    }
    let ids: Vec<_> = alive.iter().map(|&(id, p)| id * 13 + p).collect();
    json_write(Path::new(out).join("survivors.json"), &ids);
    let summary = json!({"programs":count,"scoped_rows":rows.len(),"survivors":ids.len(),"evaluations":evaluations,"first_kill_rows":killing,"complete":true,"claim":"Exact refinement on all recorded positive-base nonzero-rate PV/FV observations, including fractional and negative periods; no algorithm identification or application promotion."});
    json_write(Path::new(out).join("summary.json"), &summary);
    println!("{summary}");
}

fn refined_prediction<F: Fpu>(row: &AnnuityCase, id: usize, gs: &[Growth]) -> u64 {
    let integer = id / 13;
    let provider = id % 13;
    let gidx = integer / (3 * PROFILES);
    let a = row.numbers();
    let g = fractional_growth::<F>(&gs[gidx], a[0], a[1], provider);
    refined_prediction_with_growth::<F>(row, id, g)
}
fn refined_prediction_with_growth<F: Fpu>(row: &AnnuityCase, id: usize, g: Val) -> u64 {
    let integer = id / 13;
    if row.raw["function"] == "PV" && g.to_ext().0 == research_load(0.0).0 {
        return f64::NAN.to_bits() + 2;
    }
    publication_bits::<F>(evaluate::<F>(
        row,
        g,
        (integer / PROFILES) % 3,
        models(integer % PROFILES),
    ))
}
fn question_outcome(bits: u64) -> String {
    if bits == f64::NAN.to_bits() + 2 {
        "error:Div0".into()
    } else {
        bits_outcome(bits)
    }
}
fn question(fn_name: &str, args: [f64; 5]) -> AnnuityCase {
    let bits = args.map(f64::to_bits);
    let encoded: Vec<_> = bits.iter().map(|b| format!("0x{b:016x}")).collect();
    AnnuityCase {
        trace_index: 0,
        args: bits,
        raw: json!({"function":fn_name,"args":encoded,"result_index":Value::Null,"kind":"unanswered research request, not an Excel observation"}),
    }
}
pub fn questions(prior: &str, observations: &str, out: &str) {
    fs::create_dir(out).unwrap();
    let root = Path::new(out);
    let gs = growths();
    let ids: Vec<usize> = serde_json::from_slice(&fs::read(prior).unwrap()).unwrap();
    assert!(!ids.is_empty());
    let old = crate::lines::<AnnuityCase>(observations);
    let identity = |r: &AnnuityCase| format!("{}|{}", r.raw["function"], r.raw["args"]);
    let known: std::collections::HashSet<_> = old.iter().map(&identity).collect();
    let control_keys = [
        ("PV", "w111extended-2026092903-PV-00023"),
        ("FV", "w111extended-2026092903-FV-00023"),
        ("PV", "oc27"),
        ("FV", "rate-fv-companion-v1-0129-x-plus-h-00"),
    ];
    let mut pool = vec![];
    let mut seen = std::collections::HashSet::new();
    for (f, id) in control_keys {
        let r = old
            .iter()
            .find(|r| r.raw["function"] == f && r.raw["id"] == id)
            .unwrap();
        let q = question(f, r.numbers());
        assert!(seen.insert(identity(&q)));
        pool.push(q);
    }
    let rechecks = pool.len();
    let power = |k: i32| f64::from_bits(((1023 + k) as u64) << 52);
    let mut control_inputs = vec![];
    for r in [
        libm::nextafter(power(-53), 0.0),
        power(-53),
        libm::nextafter(power(-53), f64::INFINITY),
        libm::nextafter(-power(-54), f64::NEG_INFINITY),
        -power(-54),
        libm::nextafter(-power(-54), 0.0),
    ] {
        control_inputs.push([r, 2.0, 1.0, 1.0, 1.0]);
    }
    for n in [
        libm::nextafter(0.5, 0.0),
        0.5,
        libm::nextafter(0.5, 1.0),
        libm::nextafter(2.0, 1.0),
        2.0,
        libm::nextafter(2.0, 3.0),
    ] {
        control_inputs.push([0.05, n, 0.0, 1.0, 0.0]);
    }
    for a in control_inputs {
        for f in ["PV", "FV"] {
            let q = question(f, a);
            if seen.insert(identity(&q)) {
                pool.push(q)
            }
        }
    }
    let controls = pool.len();
    let mut rates = vec![
        power(-53),
        libm::nextafter(power(-53), 0.0),
        libm::nextafter(power(-53), f64::INFINITY),
        -power(-54),
        libm::nextafter(-power(-54), 0.0),
        libm::nextafter(-power(-54), f64::NEG_INFINITY),
        power(-28),
        libm::nextafter(power(-28), 0.0),
        libm::nextafter(power(-28), f64::INFINITY),
        power(-20),
        0.05,
        0.1,
        -0.1,
        -0.5,
        1.0,
    ];
    let mut state = 0x524f554e44385056u64;
    for _ in 0..16 {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        rates.push(f64::from_bits(
            (1022u64 << 52) | (state & 0x000fffffffffffff),
        ));
    }
    let periods = [
        1.0,
        2.0,
        3.0,
        7.0,
        12.0,
        360.0,
        0.5,
        libm::nextafter(0.5, 0.0),
        libm::nextafter(0.5, 1.0),
        libm::nextafter(2.0, 1.0),
        libm::nextafter(2.0, 3.0),
        -0.5,
        libm::nextafter(-0.5, 0.0),
        libm::nextafter(-0.5, -1.0),
        -3.0,
        -360.0,
    ];
    let balances = [
        (-100.0, 1000.0),
        (1.0, 0.0),
        (0.0, 1.0),
        (
            f64::from_bits(0x3fa3d69fba8f4fb0),
            f64::from_bits(0x3ff043028d916c0b),
        ),
    ];
    for rate in rates {
        for n in periods {
            for (p, b) in balances {
                for timing in [0.0, 1.0] {
                    for f in ["PV", "FV"] {
                        let q = question(f, [rate, n, p, b, timing]);
                        if seen.insert(identity(&q)) {
                            pool.push(q)
                        }
                    }
                }
            }
        }
    }
    // Both backends are explicit hypotheses; finite vector equality is only
    // evaluation sharing and never a claim of global programme equivalence.
    let n = ids.len();
    let mut matrix = Vec::<Vec<u64>>::new();
    let mut uncached_checks = 0usize;
    let active: std::collections::BTreeSet<_> = ids
        .iter()
        .map(|&id| (id / 13 / (3 * PROFILES), id % 13))
        .collect();
    for (idx, q) in pool.iter().enumerate() {
        let a = q.numbers();
        let mut hard = [Val::F64(0); 18 * 13];
        let mut soft = hard;
        for &(g, p) in &active {
            hard[g * 13 + p] = fractional_growth::<Hard>(&gs[g], a[0], a[1], p);
            soft[g * 13 + p] = fractional_growth::<Soft>(&gs[g], a[0], a[1], p);
        }
        let mut row: Vec<_> = ids
            .iter()
            .map(|&id| {
                refined_prediction_with_growth::<Hard>(
                    q,
                    id,
                    hard[(id / 13 / (3 * PROFILES)) * 13 + id % 13],
                )
            })
            .collect();
        row.extend(ids.iter().map(|&id| {
            refined_prediction_with_growth::<Soft>(
                q,
                id,
                soft[(id / 13 / (3 * PROFILES)) * 13 + id % 13],
            )
        }));
        if idx % 701 == 0 {
            for (j, &id) in ids.iter().enumerate() {
                assert_eq!(row[j], refined_prediction::<Hard>(q, id, &gs));
                assert_eq!(row[n + j], refined_prediction::<Soft>(q, id, &gs));
                uncached_checks += 2;
            }
        }
        matrix.push(row);
        if idx % 1000 == 0 {
            println!("question pool {idx}/{}", pool.len())
        }
    }
    let mut groups = Vec::<Vec<usize>>::new();
    let mut hashes = std::collections::HashMap::<u64, Vec<usize>>::new();
    for cid in 0..n * 2 {
        let mut hash = 0xcbf29ce484222325u64;
        for row in &matrix {
            hash = (hash ^ row[cid]).wrapping_mul(0x100000001b3);
        }
        let existing = hashes.get(&hash).and_then(|bucket| {
            bucket
                .iter()
                .copied()
                .find(|&g| matrix.iter().all(|r| r[cid] == r[groups[g][0]]))
        });
        if let Some(g) = existing {
            groups[g].push(cid)
        } else {
            let g = groups.len();
            groups.push(vec![cid]);
            hashes.entry(hash).or_default().push(g);
        }
    }
    let reps: Vec<_> = groups.iter().map(|g| g[0]).collect();
    let weights: Vec<u64> = groups.iter().map(|g| g.len() as u64).collect();
    let mut classes = vec![0usize; groups.len()];
    let mut chosen = vec![];
    let mut remaining: std::collections::BTreeSet<_> = (0..pool.len()).collect();
    for step in 0..64 {
        let pick = if step < controls {
            Some(step)
        } else {
            let mut best = None;
            let mut score_best = 0u64;
            let mut base = std::collections::HashMap::<usize, u64>::new();
            for (i, &c) in classes.iter().enumerate() {
                *base.entry(c).or_default() += weights[i]
            }
            let pairs: u64 = base.values().map(|&x| x * x).sum();
            for &qi in &remaining {
                let mut partitions = std::collections::HashMap::<(usize, u64), u64>::new();
                for (i, &cid) in reps.iter().enumerate() {
                    *partitions.entry((classes[i], matrix[qi][cid])).or_default() += weights[i]
                }
                let score = pairs - partitions.values().map(|&x| x * x).sum::<u64>();
                if score > score_best {
                    score_best = score;
                    best = Some(qi)
                }
            }
            best
        };
        let Some(qi) = pick else { break };
        remaining.remove(&qi);
        chosen.push(qi);
        let mut map = std::collections::HashMap::<(usize, u64), usize>::new();
        for (i, &cid) in reps.iter().enumerate() {
            let next = map.len();
            classes[i] = *map.entry((classes[i], matrix[qi][cid])).or_insert(next)
        }
    }
    let mut requests = std::collections::BTreeMap::<String, Vec<Value>>::new();
    let mut native =
        BufWriter::new(fs::File::create(root.join("predictions-native.jsonl")).unwrap());
    let mut software =
        BufWriter::new(fs::File::create(root.join("predictions-software.jsonl")).unwrap());
    let mut notes = vec![];
    let mut differing = 0usize;
    let mut new_inputs = 0usize;
    for (idx, &qi) in chosen.iter().enumerate() {
        let q = &pool[qi];
        let function = q.raw["function"].as_str().unwrap();
        let qid = format!("r8-annuity-{function}-{idx:02}");
        requests
            .entry(function.into())
            .or_default()
            .push(json!({"id":qid,"args":q.raw["args"]}));
        let fresh = !known.contains(&identity(q));
        new_inputs += fresh as usize;
        let mut diff = 0usize;
        for backend in 0..2 {
            let mut parts = std::collections::BTreeMap::<u64, Vec<usize>>::new();
            for (j, &id) in ids.iter().enumerate() {
                let value = matrix[qi][backend * n + j];
                parts.entry(value).or_default().push(id);
                if backend == 0 && value != matrix[qi][n + j] {
                    diff += 1
                }
            }
            let parts: Vec<_> = parts
                .into_iter()
                .map(|(v, ids)| json!({"outcome":question_outcome(v),"candidate_ids":ids}))
                .collect();
            let f = if backend == 0 {
                &mut native
            } else {
                &mut software
            };
            serde_json::to_writer(&mut *f,&json!({"id":qid,"kind":"frozen programme predictions, never Excel answers","partitions":parts})).unwrap();
            writeln!(f).unwrap();
        }
        differing += diff;
        notes.push(json!({"id":qid,"fresh_typed_identity":fresh,"intentional_recheck":idx<rechecks,"independent_substrate_control":idx>=rechecks&&idx<controls,"backend_prediction_differences":diff}));
    }
    for (function, witnesses) in requests {
        json_write(
            root.join(format!("requests-{function}.json")),
            &json!({"function":function,"witnesses":witnesses}),
        );
    }
    json_write(
        root.join("hypotheses.json"),
        &json!({"candidate_ids":ids,"layout":"id=integer_program_id*13+fractional_provider; integer_program=(growth*3+term_association)*6561+base3_site_code","backends":["native x87","software x87"],"scope":"Exact on3785positive-base nonzero-rate recorded PV/FV observations; zero and nonpositive base remain outside this unit","reader":"cargo xtask witness read --func PV --inputs requests-PV.json; separately FV. Omit --new-only because controls repeat old inputs."}),
    );
    json_write(
        root.join("selection.json"),
        &json!({"pool":pool.len(),"questions":chosen.len(),"answers":0,"fresh_typed_inputs":new_inputs,"intentional_rechecks":rechecks,"independent_substrate_controls":controls-rechecks,"uncached_prediction_checks":uncached_checks,"physical_programs":n,"backend_hypotheses":n*2,"finite_pool_prediction_signatures":groups.len(),"selected_prediction_classes":classes.iter().copied().collect::<std::collections::HashSet<_>>().len(),"backend_prediction_differences":differing,"method":"four labelled rechecks,24fixed independently suspected branch controls, then exact weighted pair disagreement over both explicit backend hypotheses; no answer labels for new requests","notes":notes}),
    );
    println!(
        "{}",
        json!({"questions":chosen.len(),"fresh":new_inputs,"programs":n,"finite_signatures":groups.len(),"backend_differences":differing})
    );
}

#[derive(Clone, Serialize, Deserialize)]
pub struct Growth {
    id: usize,
    base: EvalModel,
    product: EvalModel,
    store: bool,
}
fn growths() -> Vec<Growth> {
    let mut gs = vec![];
    for base in MODELS {
        for product in MODELS {
            for stored in [true, false] {
                gs.push(Growth {
                    id: gs.len(),
                    base,
                    product,
                    store: stored,
                });
            }
        }
    }
    gs
}
fn integer_growth<F: Fpu>(h: &Growth, r: f64, n: f64) -> Val {
    assert!(n >= 0.0 && n.fract() == 0.0 && n < u32::MAX as f64);
    let mut b = add::<F>(one(), Val::f64(r), h.base);
    let mut acc = one();
    let mut k = n as u32;
    while k > 0 {
        if k & 1 != 0 {
            acc = mul::<F>(acc, b, h.product)
        }
        k >>= 1;
        if k > 0 {
            b = mul::<F>(b, b, h.product)
        }
    }
    store::<F>(acc, h.store)
}
const PROFILES: usize = 6561;
fn models(mut code: usize) -> [EvalModel; 8] {
    std::array::from_fn(|_| {
        let m = MODELS[code % 3];
        code /= 3;
        m
    })
}
fn evaluate<F: Fpu>(row: &AnnuityCase, g: Val, association: usize, m: [EvalModel; 8]) -> Val {
    let [rate, _, payment, balance, timing] = row.numbers();
    let r = Val::f64(rate);
    let p = Val::f64(payment);
    let b = Val::f64(balance);
    let tf = if timing.trunc() == 0.0 {
        one()
    } else {
        add::<F>(one(), r, m[0])
    };
    let d = add::<F>(g, Val::f64(-1.0), m[1]);
    let term = match association {
        0 => div::<F>(mul::<F>(tf, d, m[2]), r, m[3]),
        1 => mul::<F>(tf, div::<F>(d, r, m[2]), m[3]),
        2 => mul::<F>(div::<F>(tf, r, m[2]), d, m[3]),
        _ => unreachable!(),
    };
    let paid = mul::<F>(p, term, m[4]);
    if row.raw["function"] == "FV" {
        neg(add::<F>(mul::<F>(b, g, m[5]), paid, m[6]))
    } else {
        neg(div::<F>(add::<F>(b, paid, m[6]), g, m[7]))
    }
}

pub fn joint_integer<F: Fpu>(input: &str, out: &str) {
    fs::create_dir_all(out).unwrap();
    let rows = crate::lines::<AnnuityCase>(input);
    let gs = growths();
    let count = gs.len() * 3 * PROFILES;
    let grammar = json!({"unit":"annuity-r8-unit.json","programs":count,"growths":gs,
        "constants":"0 and 1 only; exact rate==0 excluded as a separate control",
        "growth":"right-to-left binary squaring with the declared base/product models and optional final binary64 store",
        "associations":["(tf*(g-1))/rate","tf*((g-1)/rate)","(tf/rate)*(g-1)"],
        "models":MODELS,"sites":["tf addition","g-1","term first op","term second op","payment*term","FV balance*g","balance sum","PV final /g"],
        "FV":"-(balance*g + payment*term)","PV":"-(balance + payment*term)/g",
        "sharing":"One g value is shared by both occurrences and both function programs. The growth and annuity site models are fixed across every row/function. No per-row winner.",
        "layout":"id=(growth_id*3+association)*6561+code; code is little-endian base3 over the eight site models"});
    json_write(Path::new(out).join("grammar.json"), &grammar);
    let mut alive: Vec<usize> = (0..count).collect();
    let mut kills =
        BufWriter::new(fs::File::create(Path::new(out).join("first-kills.jsonl")).unwrap());
    let mut stats =
        BufWriter::new(fs::File::create(Path::new(out).join("witness-kills.jsonl")).unwrap());
    let mut evaluations = 0u64;
    let mut kill_rows = vec![];
    for (idx, row) in rows.iter().enumerate() {
        let before = alive.len();
        let mut by_outcome = std::collections::BTreeMap::<String, Vec<usize>>::new();
        if before > 0 {
            let a = row.numbers();
            let growth: Vec<_> = gs
                .iter()
                .map(|g| integer_growth::<F>(g, a[0], a[1]))
                .collect();
            alive.retain(|&id| {
                let g = id / (3 * PROFILES);
                let association = (id / PROFILES) % 3;
                let ms = models(id % PROFILES);
                let prediction = outcome::<F>(evaluate::<F>(row, growth[g], association, ms));
                evaluations += 1;
                if prediction == row.raw["expected"].as_str().unwrap() {
                    true
                } else {
                    by_outcome.entry(prediction).or_default().push(id);
                    false
                }
            });
            for (got, ids) in by_outcome {
                serde_json::to_writer(&mut kills,&json!({"candidate_ids":ids,"case":idx,"witness":row,"expected":row.raw["expected"],"got":got})).unwrap();
                writeln!(kills).unwrap();
            }
        }
        let rec = json!({"case":idx,"function":row.raw["function"],"id":row.raw["id"],"killed":before-alive.len(),"remaining":alive.len()});
        serde_json::to_writer(&mut stats, &rec).unwrap();
        writeln!(stats).unwrap();
        if before != alive.len() {
            println!("{rec}");
            kill_rows.push(rec);
            kills.flush().unwrap();
            stats.flush().unwrap()
        }
    }
    json_write(Path::new(out).join("survivors.json"), &alive);
    json_write(
        Path::new(out).join("summary.json"),
        &json!({"programs":count,"scoped_rows":rows.len(),"survivors":alive.len(),"evaluations":evaluations,"first_kill_rows":kill_rows,"complete":true,"limitation":"Integer nonnegative period positive-base PV/FV family only. A surviving finite program is a scoped hypothesis, not identification."}),
    );
    println!(
        "{}",
        json!({"programs":count,"scoped_rows":rows.len(),"survivors":alive.len(),"evaluations":evaluations})
    );
}

#[derive(Clone, Serialize, Deserialize)]
pub struct CoreHelper {
    pub id: usize,
    growth: Growth,
    difference: EvalModel,
    division: EvalModel,
    store: bool,
}
fn cores() -> Vec<CoreHelper> {
    let mut out = vec![];
    for growth in growths() {
        for difference in MODELS {
            for division in MODELS {
                for stored in [true, false] {
                    out.push(CoreHelper {
                        id: out.len(),
                        growth: growth.clone(),
                        difference,
                        division,
                        store: stored,
                    });
                }
            }
        }
    }
    out
}
fn core_values<F: Fpu>(c: &CoreHelper, r: f64, g: Val) -> (Val, Val) {
    let d = add::<F>(g, Val::f64(-1.0), c.difference);
    (g, store::<F>(div::<F>(d, Val::f64(r), c.division), c.store))
}
pub const ROOTS: usize = 729;
fn root_models(mut code: usize) -> [EvalModel; 5] {
    std::array::from_fn(|_| {
        let m = MODELS[code % 3];
        code /= 3;
        m
    })
}
fn root_value<F: Fpu>(row: &AnnuityCase, g: Val, a: Val, id: usize) -> Val {
    let m = root_models(id % 243);
    let association = id / 243;
    let [rate, _, payment, balance, timing] = row.numbers();
    let tf = if timing.trunc() == 0.0 {
        one()
    } else {
        add::<F>(one(), Val::f64(rate), m[0])
    };
    let p = Val::f64(payment);
    let b = Val::f64(balance);
    let paid = match association {
        0 => mul::<F>(mul::<F>(p, tf, m[1]), a, m[2]),
        1 => mul::<F>(mul::<F>(p, a, m[1]), tf, m[2]),
        2 => mul::<F>(p, mul::<F>(tf, a, m[1]), m[2]),
        _ => unreachable!(),
    };
    if row.raw["function"] == "FV" {
        neg(add::<F>(mul::<F>(b, g, m[3]), paid, m[4]))
    } else {
        neg(div::<F>(add::<F>(b, paid, m[3]), g, m[4]))
    }
}

pub fn factorized<F: Fpu>(input: &str, out: &str, family: u8) {
    let only_growth = family != 0;
    let wide = family == 2;
    fs::create_dir_all(out).unwrap();
    let rows = crate::lines::<AnnuityCase>(input);
    let gs = growths();
    let cs = if only_growth {
        gs.iter()
            .map(|g| CoreHelper {
                id: g.id,
                growth: g.clone(),
                difference: EvalModel::Strict,
                division: EvalModel::Strict,
                store: true,
            })
            .collect()
    } else {
        cores()
    };
    let expressions = crate::boundary::expressions_with_powers([1, 1, 1, -1]);
    let roots = if wide {
        expressions.len() * 2187
    } else if only_growth {
        6561
    } else {
        ROOTS
    };
    let count = cs.len() * roots;
    let joint_count = cs.len() as u64 * roots as u64 * roots as u64;
    json_write(
        Path::new(out).join("grammar.json"),
        &json!({"cores":cs,"roots_per_function":roots,"only_growth_shared":only_growth,
        "wide_scalar":wide,"wide_scalar_expressions":if wide{serde_json::to_value(&expressions).unwrap()}else{Value::Null},
        "wide_scalar_layout":if wide{json!({"leaves":["payment","tf","g-1","rate"],"powers":[1,1,1,-1],"root":"expression*2187+little_endian_base3_code","FV_sites":["tf","g-1","scalar op1","scalar op2","scalar op3","balance*g","sum"],"PV_sites":["tf","g-1","scalar op1","scalar op2","scalar op3","sum","/g"],"all_arithmetic_sites":7,"expression_count":expressions.len(),"claim":"All ordered four-leaf mul/div trees having exponents [1,1,1,-1]; payment can precede the division."})}else{Value::Null},
        "core":if only_growth{"shared g from binary squaring; independent annuity terms per function"}else{"shared g from binary squaring, a=store?(RN_model((g-1)/rate)); difference and division models explicit"},
        "root_models":MODELS,"root_layout":if wide{"root=scalar_expression*2187+little_endian_base3_code"}else if only_growth{"root=association*2187+little_endian_base3_code"}else{"root=association*243+little_endian_base3_code"},
        "root_sites_FV":if wide{json!(["tf","g-1","scalar op1","scalar op2","scalar op3","balance*g","sum"])}else if only_growth{json!(["tf","g-1","term op1","term op2","payment*term","balance*g","sum"])}else{json!(["tf addition","term op1","term op2","balance*g","sum"])},
        "root_sites_PV":if wide{json!(["tf","g-1","scalar op1","scalar op2","scalar op3","sum","/g"])}else if only_growth{json!(["tf","g-1","term op1","term op2","payment*term","sum","/g"])}else{json!(["tf addition","term op1","term op2","sum","/g"])},
        "term_associations":if wide{json!({"catalogue":"wide_scalar_expressions","leaves":["payment","tf","g-1","rate"]})}else if only_growth{json!(["(tf*(g-1))/rate","tf*((g-1)/rate)","(tf/rate)*(g-1)"])}else{json!(["(payment*tf)*a","(payment*a)*tf","payment*(tf*a)"])},
        "FV":"-(balance*g + term)","PV":"-(balance + term)/g",
        "projected_profile_layout":format!("profile=core*{roots}+root"),
        "joint_layout":format!("joint=(core*{roots}+FV_root)*{roots}+PV_root"),
        "growth_only_layout":if only_growth&&!wide{json!({"root":"association*2187+little_endian_base3_code","FV_sites":["tf","g-1","term op1","term op2","payment*term","balance*g","sum"],"PV_sites":["tf","g-1","term op1","term op2","payment*term","sum","/g"],"term_associations":["(tf*(g-1))/rate","tf*((g-1)/rate)","(tf/rate)*(g-1)"],"profile":"core*6561+root","joint":"(core*6561+FV_root)*6561+PV_root","sharing":"Only g is shared; each function has its own untimed/timed term and root operations. The unused CoreHelper difference/division/store fields do not participate."})}else{Value::Null},
        "programs":joint_count,"projected_profiles_per_function":count,
        "factorization_proof":"For fixed core parameters PV and FV roots have disjoint independent parameters; a global program survives iff both projected roots survive all their own function's rows. The joint survivor set is exactly their Cartesian product per core. The declared shared core is a fixed program for every row/function; no per-row helper values are chosen."}),
    );
    let mut alive = [
        (0..count).collect::<Vec<usize>>(),
        (0..count).collect::<Vec<usize>>(),
    ];
    let mut kills =
        BufWriter::new(fs::File::create(Path::new(out).join("first-kills.jsonl")).unwrap());
    let mut stats =
        BufWriter::new(fs::File::create(Path::new(out).join("witness-kills.jsonl")).unwrap());
    let mut evaluations = 0u64;
    let mut killing = vec![];
    for (idx, row) in rows.iter().enumerate() {
        let which = if row.raw["function"] == "FV" { 0 } else { 1 };
        let before = alive[which].len();
        let mut rejected = std::collections::BTreeMap::<String, Vec<usize>>::new();
        if before > 0 {
            let numbers = row.numbers();
            let growth: Vec<_> = gs
                .iter()
                .map(|g| integer_growth::<F>(g, numbers[0], numbers[1]))
                .collect();
            let vals: Vec<_> = cs
                .iter()
                .map(|c| core_values::<F>(c, numbers[0], growth[c.growth.id]))
                .collect();
            alive[which].retain(|&id| {
                let (g, a) = vals[id / roots];
                let v = if wide {
                    let root = id % roots;
                    let mut code = root % 2187;
                    let ms: [EvalModel; 7] = std::array::from_fn(|_| {
                        let m = MODELS[code % 3];
                        code /= 3;
                        m
                    });
                    let [r, _, p, b, t] = row.numbers();
                    let tf = if t.trunc() == 0.0 {
                        one()
                    } else {
                        add::<F>(one(), Val::f64(r), ms[0])
                    };
                    let d = add::<F>(g, Val::f64(-1.0), ms[1]);
                    let mut site = 0;
                    let paid = crate::boundary::eval_expr::<F>(
                        &expressions[root / 2187],
                        &[Val::f64(p), tf, d, Val::f64(r)],
                        &[ms[2], ms[3], ms[4], EvalModel::Strict],
                        &mut site,
                    );
                    assert_eq!(site, 3);
                    if which == 0 {
                        neg(add::<F>(mul::<F>(Val::f64(b), g, ms[5]), paid, ms[6]))
                    } else {
                        neg(div::<F>(add::<F>(Val::f64(b), paid, ms[5]), g, ms[6]))
                    }
                } else if only_growth {
                    let root = id % roots;
                    let mut code = root % 2187;
                    let digits: [EvalModel; 7] = std::array::from_fn(|_| {
                        let m = MODELS[code % 3];
                        code /= 3;
                        m
                    });
                    let ms = if which == 0 {
                        [
                            digits[0],
                            digits[1],
                            digits[2],
                            digits[3],
                            digits[4],
                            digits[5],
                            digits[6],
                            EvalModel::Strict,
                        ]
                    } else {
                        [
                            digits[0],
                            digits[1],
                            digits[2],
                            digits[3],
                            digits[4],
                            EvalModel::Strict,
                            digits[5],
                            digits[6],
                        ]
                    };
                    evaluate::<F>(row, g, root / 2187, ms)
                } else {
                    root_value::<F>(row, g, a, id % roots)
                };
                let got = outcome::<F>(v);
                evaluations += 1;
                if got == row.raw["expected"].as_str().unwrap() {
                    true
                } else {
                    rejected.entry(got).or_default().push(id);
                    false
                }
            });
            for (got, ids) in rejected {
                serde_json::to_writer(&mut kills,&json!({"function":row.raw["function"],"projected_profile_ids":ids,"case":idx,"witness":row,"expected":row.raw["expected"],"got":got})).unwrap();
                writeln!(kills).unwrap();
            }
        }
        let rec = json!({"case":idx,"function":row.raw["function"],"id":row.raw["id"],"killed":before-alive[which].len(),"remaining_projected":alive[which].len()});
        serde_json::to_writer(&mut stats, &rec).unwrap();
        writeln!(stats).unwrap();
        if before != alive[which].len() {
            println!("{rec}");
            killing.push(rec);
            kills.flush().unwrap();
            stats.flush().unwrap()
        }
    }
    let mut joint = vec![];
    let mut total = 0u64;
    for c in &cs {
        let selected: Vec<Vec<usize>> = alive
            .iter()
            .map(|ids| {
                ids.iter()
                    .filter_map(|&id| (id / roots == c.id).then_some(id % roots))
                    .collect()
            })
            .collect();
        let n = selected[0].len() as u64 * selected[1].len() as u64;
        total += n;
        joint.push(
            json!({"core":c.id,"FV_roots":selected[0],"PV_roots":selected[1],"joint_survivors":n}),
        );
    }
    json_write(Path::new(out).join("joint-survivors.json"), &joint);
    for (i, name) in ["FV", "PV"].iter().enumerate() {
        json_write(
            Path::new(out).join(format!("{name}-survivors.json")),
            &alive[i],
        );
    }
    let summary = json!({"programs":joint_count,"projected_profiles_per_function":count,"surviving_FV_profiles":alive[0].len(),"surviving_PV_profiles":alive[1].len(),"joint_survivors":total,"scoped_rows":rows.len(),"evaluations":evaluations,"first_kill_rows":killing,"complete":true,"claim":"Exact factorized joint constraints, not count fitting. Shared growth and annuity core, independent function-specific timing/payment/root operation graphs."});
    json_write(Path::new(out).join("summary.json"), &summary);
    println!("{summary}");
}

pub(crate) fn pmt_growth_value<F: Fpu>(
    r: &crate::Row,
    g: Val,
    em: Val,
    association: usize,
    mut code: usize,
) -> Val {
    let m: [EvalModel; 5] = std::array::from_fn(|_| {
        let v = MODELS[code % 3];
        code /= 3;
        v
    });
    let [rate, pv, _, timing, _] = r.args();
    let numerator = mul::<F>(Val::f64(pv), g, m[0]);
    let tf = if timing == 0.0 {
        one()
    } else {
        add::<F>(one(), Val::f64(rate), m[1])
    };
    let v = match association {
        0 => mul::<F>(
            div::<F>(div::<F>(numerator, em, m[2]), tf, m[3]),
            Val::f64(rate),
            m[4],
        ),
        1 => div::<F>(
            mul::<F>(numerator, Val::f64(rate), m[2]),
            mul::<F>(em, tf, m[3]),
            m[4],
        ),
        2 => mul::<F>(
            div::<F>(numerator, tf, m[2]),
            div::<F>(Val::f64(rate), em, m[3]),
            m[4],
        ),
        _ => unreachable!(),
    };
    neg(v)
}

pub fn pmt_growth_search<F: Fpu>(input: &str, joint_pv_fv: &str, out: &str) {
    fs::create_dir_all(out).unwrap();
    let root = Path::new(out);
    let original = crate::lines::<crate::Row>(input);
    let gs = growths();
    let hs = crate::boundary::helpers();
    let pv_fv: Vec<usize> = serde_json::from_slice(&fs::read(joint_pv_fv).unwrap()).unwrap();
    let mut by_growth = vec![0u64; gs.len()];
    for id in &pv_fv {
        by_growth[id / (3 * PROFILES)] += 1;
    }
    let selected: Vec<_> = gs.iter().filter(|g| by_growth[g.id] > 0).cloned().collect();
    let mut rows = vec![];
    let mut exclusions = std::collections::BTreeMap::<String, usize>::new();
    let mut excluded = BufWriter::new(fs::File::create(root.join("excluded.jsonl")).unwrap());
    for r in original {
        let n = number(&r.raw["args"][1]);
        let a = r.args();
        let reason = if a[2] != 0.0 {
            "FV_nonzero_control"
        } else if n <= 0.0 || n.fract() != 0.0 || n >= u32::MAX as f64 {
            "nonpositive_or_fractional_period_control"
        } else if selected
            .iter()
            .any(|g| !integer_growth::<F>(g, a[0], n).to_f64::<F>().is_finite())
        {
            "nonfinite_growth_control"
        } else {
            "scoped"
        };
        *exclusions.entry(reason.into()).or_default() += 1;
        if reason == "scoped" {
            rows.push(r)
        } else {
            serde_json::to_writer(&mut excluded, &json!({"reason":reason,"row":r})).unwrap();
            writeln!(excluded).unwrap();
        }
    }
    let per_growth = hs.len() * 3 * 243;
    let count = selected.len() * per_growth;
    let joint_count = pv_fv.len() as u64 * per_growth as u64;
    json_write(
        root.join("grammar.json"),
        &json!({"growths":selected,"positive_expm1_helpers":hs,
        "PV_FV_precondition_survivors":pv_fv.len(),"PV_FV_profiles_per_growth":by_growth,
        "shared_growth":"Same integer binary-squaring program in PV,FV and PMT",
        "helper_contract":"boundary::helper(h, rate, -nper, frozen_log1p(rate)); positive exponent nper*log1p(rate), no refit",
        "numerator":"PV*growth; FV=0 scope, with strict/stored-PC64/resident-PC64 product",
        "models":MODELS,"sites":["PV*growth","tf addition","combine op1","combine op2","combine op3"],
        "combines":["((N/em)/tf)*rate","(N*rate)/(em*tf)","(N/tf)*(rate/em)"],
        "publication":"negate then same PMT outcome policy; full Ext80 until declared store/final binary64 publication",
        "component_programs":count,"joint_programs_with_PV_FV":joint_count,
        "layout":"component=((selected_growth_index*1980+helper_id)*3+combine)*243+little_endian_base3_code",
        "scope_counts":exclusions,"all_inputs":"Original u64 PMT leaves and raw bit-encoded nper, no decimal normalization"}),
    );
    let mut alive: Vec<usize> = (0..count).collect();
    let mut evaluations = 0u64;
    let mut kills = BufWriter::new(fs::File::create(root.join("first-kills.jsonl")).unwrap());
    let mut stats = BufWriter::new(fs::File::create(root.join("witness-kills.jsonl")).unwrap());
    let mut killers = vec![];
    for (idx, row) in rows.iter().enumerate() {
        let before = alive.len();
        let mut rejected = std::collections::BTreeMap::<String, Vec<usize>>::new();
        if before > 0 {
            let n = number(&row.raw["args"][1]);
            let a = row.args();
            let growth: Vec<_> = selected
                .iter()
                .map(|g| integer_growth::<F>(g, a[0], n))
                .collect();
            let em: Vec<_> = hs
                .iter()
                .map(|h| crate::boundary::helper::<F>(h, a[0], -n, f64::from_bits(row.log1p)))
                .collect();
            alive.retain(|&id| {
                let gi = id / per_growth;
                let hi = (id / (3 * 243)) % hs.len();
                let association = (id / 243) % 3;
                let got = outcome::<F>(pmt_growth_value::<F>(
                    row,
                    growth[gi],
                    em[hi],
                    association,
                    id % 243,
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
            killers.push(rec);
            kills.flush().unwrap();
            stats.flush().unwrap()
        }
    }
    let mut joint_survivors = 0u64;
    for id in &alive {
        joint_survivors += by_growth[selected[id / per_growth].id];
    }
    json_write(root.join("survivors.json"), &alive);
    let summary = json!({"component_programs":count,"joint_programs_with_PV_FV":joint_count,"component_survivors":alive.len(),"joint_survivors":joint_survivors,"scope_counts":exclusions,"scoped_rows":rows.len(),"evaluations":evaluations,"first_kill_rows":killers,"complete":true,"claim":"Exact early-numerator structure test jointly constrained by PV/FV integer growth; no fitted parameters, no application promotion."});
    json_write(root.join("summary.json"), &summary);
    println!("{summary}");
}
