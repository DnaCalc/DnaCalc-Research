//! Finite variance programs from defining identities and published updating formulas.
//! Original implementation, MIT. See variance-r14-unit.json for scope and sources.
//! This module is research only; no candidate is selected by residual or match count.
use crate::dsl::{EvalModel, Op, WIDE_MODELS};
use crate::eval::{self, Val};
use crate::ported::excel_numeric::Fpu;
use calc_func::value::CoreRef;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::collections::BTreeMap;
use std::fs::{self, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Profile {
    id: usize,
    algorithm: usize,
    linear_order: usize,
    square_order: usize,
    // Linear reduction; mean/residual; products; squared reduction; finish.
    models: [EvalModel; 5],
}

pub fn profiles() -> Vec<Profile> {
    let mut all = Vec::new();
    for algorithm in 0..14 {
        let orders: Vec<_> = if algorithm < 5 || algorithm >= 7 {
            (0..3).flat_map(|a| (0..3).map(move |b| (a, b))).collect()
        } else {
            vec![(0, 0), (1, 1)]
        };
        for (linear_order, square_order) in orders {
            for index in 0..5usize.pow(5) {
                let mut label = index;
                let models = std::array::from_fn(|_| {
                    let model = WIDE_MODELS[label % 5];
                    label /= 5;
                    model
                });
                all.push(Profile {
                    id: all.len(),
                    algorithm,
                    linear_order,
                    square_order,
                    models,
                });
            }
        }
    }
    assert_eq!(all.len(), 350000);
    all
}

/// AVERAGE has just a reduction and division. Other regions do not enter its grammar.
fn function_profiles(function: &str) -> Vec<Profile> {
    if function != "AVERAGE" { return profiles(); }
    let mut all = Vec::new();
    for linear_order in 0..3 {
        for &division in &WIDE_MODELS {
            for &sum in &WIDE_MODELS {
                all.push(Profile { id: all.len(), algorithm: 0, linear_order,
                    square_order: 0, models: [sum, division, EvalModel::Strict, EvalModel::Strict, EvalModel::Strict] });
            }
        }
    }
    all
}

// A separate, explicitly named family preserves the R14/R16 grammar and its
// receipts. Structural choices are mathematical algorithms, not answer rules.
fn structural_profiles(function: &str) -> Vec<Profile> {
    let mut all = Vec::new();
    let algorithms = if function == "AVERAGE" { 1 } else { 22 };
    for algorithm in 0..algorithms {
        let orders: Vec<_> = if [5, 6, 19, 20, 21].contains(&algorithm) {
            vec![(0, 0), (1, 1)]
        } else if function == "AVERAGE" {
            (0..13).map(|order| (order, 0)).collect()
        } else {
            (0..3).flat_map(|x| (0..3).map(move |y| (x, y)))
                .chain((3..13).map(|order| (order, order))).collect()
        };
        for (linear_order, square_order) in orders {
            let policies = if function == "AVERAGE" { 25 } else { 5usize.pow(5) };
            for index in 0..policies {
                let mut label = index;
                let models = std::array::from_fn(|region| {
                    if function == "AVERAGE" && region >= 2 { return EvalModel::Strict; }
                    let model = WIDE_MODELS[label % 5];
                    label /= 5;
                    model
                });
                all.push(Profile { id: all.len(), algorithm, linear_order, square_order, models });
            }
        }
    }
    all
}

fn neg(x: Val) -> Val {
    match x {
        Val::F64(b) => Val::F64(b ^ (1 << 63)),
        Val::Ext(mut e) => {
            e.0[9] ^= 128;
            Val::Ext(e)
        }
    }
}

struct Arithmetic<'a, F> {
    p: &'a Profile,
    trace: Option<Vec<Value>>,
    marker: std::marker::PhantomData<F>,
}
impl<F: Fpu> Arithmetic<'_, F> {
    fn op(&mut self, op: Op, region: usize, a: Val, b: Val) -> Val {
        let z = eval::binop::<F>(op, self.p.models[region], a, b);
        if let Some(t) = &mut self.trace {
            t.push(json!({"op":match op {Op::Add(..)=>"add",Op::Mul(..)=>"mul",Op::Div(..)=>"div",_=>unreachable!()},
                "region":region,"a":a.exact_encoding(),"b":b.exact_encoding(),"value":z.exact_encoding()}));
        }
        z
    }
    fn add(&mut self, region: usize, a: Val, b: Val) -> Val {
        self.op(Op::Add(0, 1), region, a, b)
    }
    fn sub(&mut self, region: usize, a: Val, b: Val) -> Val {
        self.add(region, a, neg(b))
    }
    fn mul(&mut self, region: usize, a: Val, b: Val) -> Val {
        self.op(Op::Mul(0, 1), region, a, b)
    }
    fn div(&mut self, region: usize, a: Val, b: Val) -> Val {
        self.op(Op::Div(0, 1), region, a, b)
    }
    fn sum(&mut self, xs: &[Val], order: usize, region: usize) -> Val {
        match order {
            0 => xs
                .iter()
                .fold(Val::f64(0.0), |sum, &x| self.add(region, sum, x)),
            1 => xs
                .iter()
                .rev()
                .fold(Val::f64(0.0), |sum, &x| self.add(region, sum, x)),
            2 => match xs.len() {
                0 => Val::f64(0.0),
                1 => xs[0],
                n => {
                    let a = self.sum(&xs[..n / 2], 2, region);
                    let b = self.sum(&xs[n / 2..], 2, region);
                    self.add(region, a, b)
                }
            },
            3 | 4 => {
                let mut sorted = xs.to_vec();
                sorted.sort_by(|a, b| a.to_f64::<F>().abs().total_cmp(&b.to_f64::<F>().abs()));
                self.sum(&sorted, usize::from(order == 4), region)
            }
            5 | 6 => {
                let mut total = Val::f64(0.0);
                let mut correction = Val::f64(0.0);
                let mut data = xs.to_vec();
                if order == 6 { data.reverse(); }
                for x in data {
                    let y = self.sub(region, x, correction);
                    let next = self.add(region, total, y);
                    let added = self.sub(region, next, total);
                    correction = self.sub(region, added, y);
                    total = next;
                }
                total
            }
            7 | 8 => {
                let mut total = Val::f64(0.0);
                let mut correction = Val::f64(0.0);
                let mut data = xs.to_vec();
                if order == 8 { data.reverse(); }
                for x in data {
                    let next = self.add(region, total, x);
                    let error = if total.to_f64::<F>().abs() >= x.to_f64::<F>().abs() {
                        let retained = self.sub(region, total, next);
                        self.add(region, retained, x)
                    } else {
                        let retained = self.sub(region, x, next);
                        self.add(region, retained, total)
                    };
                    correction = self.add(region, correction, error);
                    total = next;
                }
                self.add(region, total, correction)
            }
            9..=11 => {
                let lanes = 1 << (order - 8);
                let mut totals = vec![Val::f64(0.0); lanes];
                for (i, &x) in xs.iter().enumerate() {
                    totals[i % lanes] = self.add(region, totals[i % lanes], x);
                }
                self.sum(&totals, 0, region)
            }
            12 => {
                let mut values = xs.to_vec();
                while values.len() > 1 {
                    values = values.chunks(2).map(|pair| if pair.len() == 2 {
                        self.add(region, pair[0], pair[1])
                    } else { pair[0] }).collect();
                }
                values.first().copied().unwrap_or(Val::f64(0.0))
            }
            _ => unreachable!(),
        }
    }
    fn sqrt(&mut self, mut x: Val) -> Val {
        let model = self.p.models[4];
        let z = if model == EvalModel::Strict {
            x = Val::f64(x.to_f64::<F>());
            Val::f64(x.to_f64::<F>().sqrt())
        } else {
            let cw = if matches!(model, EvalModel::StoredPc53 | EvalModel::ResidentPc53) {
                0x123f
            } else {
                0x133f
            };
            let z = F::sqrt(&x.to_ext(), cw);
            if matches!(model, EvalModel::StoredPc53 | EvalModel::StoredPc64) {
                Val::f64(F::store(&z, 0x133f))
            } else {
                Val::Ext(z)
            }
        };
        if let Some(t) = &mut self.trace {
            t.push(
                json!({"op":"sqrt","region":4,"a":x.exact_encoding(),"value":z.exact_encoding()}),
            );
        }
        z
    }
}

fn population(function: &str) -> bool {
    matches!(
        function,
        "VAR.P" | "VARP" | "VARPA" | "STDEV.P" | "STDEVP" | "STDEVPA"
    )
}

/// Compare finite tagged values without inserting a binary64 store. Stored
/// policies already contain binary64; a resident guard must keep its low bits.
fn finite_less(a:Val,b:Val)->bool {
    let parts=|v:Val| {
        let x=v.to_ext().0;
        let exponent=u16::from_le_bytes([x[8],x[9]])&0x7fff;
        assert_ne!(exponent,0x7fff,"variance counterfactual outside finite comparison scope");
        let significand=u64::from_le_bytes(x[..8].try_into().unwrap());
        (x[9]&0x80!=0,(exponent,significand))
    };
    let (sa,ma)=parts(a);let (sb,mb)=parts(b);
    if ma==(0,0)&&mb==(0,0){return false;}
    if sa!=sb{return sa;}
    if sa{ma>mb}else{ma<mb}
}

#[cfg(test)]
mod guard_tests {
    use super::*;
    #[test]
    fn resident_guard_keeps_bits_that_a_double_store_would_lose() {
        let one=Val::f64(1.0);
        let above=eval::binop::<crate::ported::excel_numeric::Soft>(Op::Add(0,1),
            EvalModel::ResidentPc64,one,Val::f64(2f64.powi(-60)));
        assert_eq!(one.to_f64::<crate::ported::excel_numeric::Soft>(),
                   above.to_f64::<crate::ported::excel_numeric::Soft>());
        assert!(finite_less(one,above));assert!(!finite_less(above,one));
        assert!(finite_less(neg(above),neg(one)));
        assert!(!finite_less(Val::f64(-0.0),Val::f64(0.0)));
    }
}

/// Raw moments with the previously measured cancellation fallback. This is a
/// structural counterfactual, not a per-row choice between expected answers.
pub fn stored_steps<F: Fpu>(input: &str, output: &str) {
    use EvalModel::{Strict, StoredPc64, ResidentPc64};
    let rows = crate::lines::<Value>(input);
    let mut counts = [[0usize; 6]; 9];
    let mut changes = Vec::new();
    let mut predictions = Vec::new();
    let mut excluded = BTreeMap::<String,usize>::new();
    for (index, raw) in rows.iter().enumerate() {
        let (raw, bits) = if raw.get("values_bits").is_some() {
            (&raw["raw"], serde_json::from_value::<Vec<u64>>(raw["values_bits"].clone()).unwrap())
        } else {
            match numeric_values(raw) {
                Ok(bits) => (raw, bits),
                Err(reason) => { *excluded.entry(reason.into()).or_default() += 1; continue; }
            }
        };
        let function = raw["function"].as_str().unwrap();
        let pop = population(function);
        let mut outputs = Vec::new();
        for (id, (sum, step)) in [Strict,StoredPc64,ResidentPc64].into_iter()
            .flat_map(|sum| [Strict,StoredPc64,ResidentPc64].map(|step|(sum,step))).enumerate() {
            let p = Profile { id, algorithm: if pop {8} else {7}, linear_order:0,
                square_order:0, models:[sum,step,step,sum,step] };
            let mut a = Arithmetic::<F> {p:&p,trace:None,marker:std::marker::PhantomData};
            let xs:Vec<_> = bits.iter().copied().map(Val::F64).collect();
            let total = a.sum(&xs,0,0);
            let terms:Vec<_> = xs.iter().map(|&x|a.mul(2,x,x)).collect();
            let q = a.sum(&terms,0,3);
            let n = Val::f64(bits.len() as f64);
            let square = a.mul(2,total,total);
            let left = if pop {a.mul(4,n,q)} else {q};
            let right = if pop {square} else {a.div(1,square,n)};
            let numerator = a.sub(4,left,right);
            let cutoff = a.mul(4,left,Val::f64(0.01));
            let fallback = finite_less(numerator,cutoff);
            let mut p = p;
            if fallback {p.algorithm=0;}
            let prediction = evaluate::<F>(&p,&bits,function,false).0.to_f64::<F>();
            let actual = if prediction.is_finite() && !prediction.is_subnormal() {
                format!("0x{:016x}", if prediction==0.0 {0} else {prediction.to_bits()})
            } else {"error:Num".into()};
            let old = raw["actual"] == raw["expected"];
            let now = raw["expected"] == actual;
            for (i, flag) in [true,now,old,now&&!old,old&&!now,fallback].into_iter().enumerate() {
                counts[id][i] += usize::from(flag);
            }
            if actual != raw["actual"] {
                changes.push(json!({"program":id,"trace_index":index,"id":raw["id"],
                    "function":function,"args":raw["args"],"expected":raw["expected"],
                    "baseline":raw["actual"],"prediction":actual,"fallback":fallback,
                    "ratio":numerator.to_f64::<F>()/left.to_f64::<F>()}));
            }
            outputs.push(actual);
        }
        let center = Profile {id:0,algorithm:0,linear_order:0,square_order:0,
            models:[StoredPc64;5]};
        let centered = evaluate::<F>(&center,&bits,function,false).0.to_f64::<F>();
        predictions.push(json!({"trace_index":index,"id":raw["id"],"predictions":outputs,
            "centered_pc64":format!("0x{:016x}",centered.to_bits())}));
    }
    let destination = Path::new(output); fs::create_dir(destination).unwrap();
    crate::json_write(destination.join("changes.json"),&changes);
    crate::json_write(destination.join("predictions.json"),&predictions);
    let summary=json!({"complete":true,"counts_all_exact_baseline_gains_losses_fallback":counts,
        "excluded":excluded,"programs":"Cartesian (sum,step) of strict, stored-PC64, resident-PC64; forward folds; sample (Q-S*S/n)/(n-1), population (n*Q-S*S)/(n*n); centered two-pass when the numerator <0.01*left.",
        "limits":"The existing corpus brackets the guard near1%; its exact transition needs new Excel questions. No answer-dependent dispatch; resident models are counterfactuals. Original answers remain unchanged."});
    crate::json_write(destination.join("summary.json"),&summary); println!("{summary}");
}

/// Every region applies the same declared precision/store policy at every loop iteration.
pub fn evaluate<F: Fpu>(
    p: &Profile,
    bits: &[u64],
    function: &str,
    traced: bool,
) -> (Val, Option<Vec<Value>>) {
    let mut a = Arithmetic::<F> {
        p,
        trace: traced.then(Vec::new),
        marker: std::marker::PhantomData,
    };
    let mut xs: Vec<_> = bits.iter().copied().map(Val::F64).collect();
    let n = Val::f64(xs.len() as f64);
    if function == "AVERAGE" {
        let total = a.sum(&xs, p.linear_order, 0);
        let mean = a.div(1, total, n);
        return (mean, a.trace);
    }
    let d = Val::f64((xs.len() - usize::from(!population(function))) as f64);
    if p.algorithm >= 14 {
        let ss = if p.algorithm == 21 {
            // Chan's merge identity: M2(A+B)=M2(A)+M2(B)+delta^2*nA*nB/(nA+nB).
            fn merge<F: Fpu>(a: &mut Arithmetic<'_, F>, xs: &[Val]) -> (Val, Val) {
                if xs.len() == 1 { return (xs[0], Val::f64(0.0)); }
                let k = xs.len()/2;
                let (left, left_ss) = merge(a, &xs[..k]);
                let (right, right_ss) = merge(a, &xs[k..]);
                let delta = a.sub(1, right, left);
                let fraction = a.div(1, Val::f64((xs.len()-k) as f64), Val::f64(xs.len() as f64));
                let step = a.mul(1, delta, fraction);
                let mean = a.add(1, left, step);
                let square = a.mul(2, delta, delta);
                let weighted = a.mul(2, square, Val::f64((k*(xs.len()-k)) as f64));
                let correction = a.div(4, weighted, Val::f64(xs.len() as f64));
                let subtotal = a.add(3, left_ss, right_ss);
                (mean, a.add(3, subtotal, correction))
            }
            if p.linear_order == 1 { xs.reverse(); }
            merge(&mut a, &xs).1
        } else if p.algorithm >= 19 {
            if p.linear_order == 1 { xs.reverse(); }
            let mut mean = xs[0];
            let mut ss = Val::f64(0.0);
            for (i, &x) in xs.iter().enumerate().skip(1) {
                let k = Val::f64((i+1) as f64);
                let delta = a.sub(1, x, mean);
                if p.algorithm == 19 {
                    let weighted = a.mul(1, mean, Val::f64(i as f64));
                    let numerator = a.add(1, weighted, x);
                    mean = a.div(1, numerator, k);
                    let remainder = a.sub(1, x, mean);
                    let term = a.mul(2, delta, remainder);
                    ss = a.add(3, ss, term);
                } else {
                    let increment = a.div(1, delta, k);
                    mean = a.add(1, mean, increment);
                    let square = a.mul(2, delta, delta);
                    let weighted = a.mul(2, square, Val::f64(i as f64));
                    let term = a.div(4, weighted, k);
                    ss = a.add(3, ss, term);
                }
            }
            ss
        } else {
            if p.algorithm == 15 {
                let first = xs[0];
                xs = xs.into_iter().map(|x| a.sub(1, x, first)).collect();
            }
            let total = a.sum(&xs, p.linear_order, 0);
            let mut mean = a.div(1, total, n);
            if [14, 18].contains(&p.algorithm) {
                // Correct the computed center once or twice, then perform a
                // new pass. No convergence threshold is chosen from answers.
                for _ in 0..if p.algorithm == 14 { 1 } else { 2 } {
                    let deltas: Vec<_> = xs.iter().map(|&x| a.sub(1, x, mean)).collect();
                    let residual = a.sum(&deltas, p.linear_order, 0);
                    let correction = a.div(1, residual, n);
                    mean = a.add(1, mean, correction);
                }
            }
            let deltas: Vec<_> = xs.iter().map(|&x| a.sub(1, x, mean)).collect();
            let terms: Vec<_> = deltas.iter().zip(&xs).map(|(&delta, &x)|
                a.mul(2, delta, if [16, 17].contains(&p.algorithm) { x } else { delta })).collect();
            let mut ss = a.sum(&terms, p.square_order, 3);
            if p.algorithm == 15 || p.algorithm == 17 {
                let residual = a.sum(&deltas, p.linear_order, 0);
                let correction = if p.algorithm == 15 {
                    let square = a.mul(2, residual, residual);
                    a.div(4, square, n)
                } else { a.mul(2, mean, residual) };
                ss = a.sub(4, ss, correction);
            }
            ss
        };
        let mut result = if function == "DEVSQ" { a.add(4, ss, Val::f64(0.0)) } else { a.div(4, ss, d) };
        if function.starts_with("STDEV") { result = a.sqrt(result); }
        return (result, a.trace);
    }
    // Algebraically equal raw-moment forms have distinct store boundaries.
    // Constants here are only the count and its sample/population divisor.
    if p.algorithm >= 7 {
        let total = a.sum(&xs, p.linear_order, 0);
        let terms: Vec<_> = xs.iter().map(|&x| a.mul(2, x, x)).collect();
        let squares = a.sum(&terms, p.square_order, 3);
        let square_total = a.mul(2, total, total);
        let mut value = match p.algorithm {
            7 => {
                let correction = a.div(1, square_total, n);
                let numerator = a.sub(4, squares, correction);
                if function == "DEVSQ" { numerator } else { a.div(4, numerator, d) }
            }
            8 => {
                let left = a.mul(4, n, squares);
                let numerator = a.sub(4, left, square_total);
                let divisor = if function == "DEVSQ" { n } else { a.mul(4, n, d) };
                a.div(4, numerator, divisor)
            }
            9 => {
                let divisor = if function == "DEVSQ" { Val::f64(1.0) } else { d };
                let left = a.div(4, squares, divisor);
                let denominator = a.mul(4, n, divisor);
                let right = a.div(4, square_total, denominator);
                a.sub(4, left, right)
            }
            10 => {
                let divisor = if function == "DEVSQ" { Val::f64(1.0) } else { d };
                let left = a.div(4, squares, divisor);
                let mean = a.div(1, total, n);
                let scaled = a.div(1, total, divisor);
                let right = a.mul(2, mean, scaled);
                a.sub(4, left, right)
            }
            11 => {
                let mean = a.div(1, total, n);
                let left = a.div(4, squares, n);
                let right = a.mul(2, mean, mean);
                let difference = a.sub(4, left, right);
                let scale = if function == "DEVSQ" { n } else { a.div(4, n, d) };
                a.mul(4, difference, scale)
            }
            12 => {
                let mean = a.div(1, total, n);
                let terms: Vec<_> = xs.iter().map(|&x| {
                    let delta = a.sub(1, x, mean);
                    let product = a.mul(2, delta, delta);
                    if function == "DEVSQ" { product } else { a.div(4, product, d) }
                }).collect();
                a.sum(&terms, p.square_order, 3)
            }
            13 => {
                let mut terms = Vec::new();
                for i in 0..xs.len() {
                    for j in i+1..xs.len() {
                        let difference = a.sub(1, xs[i], xs[j]);
                        terms.push(a.mul(2, difference, difference));
                    }
                }
                let numerator = a.sum(&terms, p.square_order, 3);
                let divisor = if function == "DEVSQ" { n } else { a.mul(4, n, d) };
                a.div(4, numerator, divisor)
            }
            _ => unreachable!(),
        };
        if function.starts_with("STDEV") { value = a.sqrt(value); }
        return (value, a.trace);
    }
    let ss = if p.algorithm >= 5 {
        if p.linear_order == 1 {
            xs.reverse();
        }
        let mut center = xs[0];
        let mut m2 = Val::f64(0.0);
        let mut total = xs[0];
        for (i, &x) in xs.iter().enumerate().skip(1) {
            let count = Val::f64((i + 1) as f64);
            if p.algorithm == 5 {
                let delta = a.sub(1, x, center);
                let step = a.div(1, delta, count);
                center = a.add(1, center, step);
                let after = a.sub(1, x, center);
                let term = a.mul(2, delta, after);
                m2 = a.add(3, m2, term);
            } else {
                total = a.add(0, total, x);
                let kx = a.mul(1, count, x);
                let delta = a.sub(1, kx, total);
                let square = a.mul(2, delta, delta);
                let den = a.mul(4, count, Val::f64(i as f64));
                let term = a.div(4, square, den);
                m2 = a.add(3, m2, term);
            }
        }
        m2
    } else {
        if p.algorithm == 4 {
            let anchor = xs[0];
            xs = xs.into_iter().map(|x| a.sub(1, x, anchor)).collect();
        }
        let total = a.sum(&xs, p.linear_order, 0);
        let mean = a.div(1, total, n);
        if p.algorithm < 2 {
            let deltas: Vec<_> = xs.iter().map(|&x| a.sub(1, x, mean)).collect();
            let terms: Vec<_> = deltas.iter().map(|&x| a.mul(2, x, x)).collect();
            let squares = a.sum(&terms, p.square_order, 3);
            if p.algorithm == 0 {
                squares
            } else {
                let residual = a.sum(&deltas, p.square_order, 0);
                let product = a.mul(2, residual, residual);
                let correction = a.div(4, product, n);
                a.sub(4, squares, correction)
            }
        } else {
            let terms: Vec<_> = xs.iter().map(|&x| a.mul(2, x, x)).collect();
            let squares = a.sum(&terms, p.square_order, 3);
            if p.algorithm == 3 {
                let left = a.mul(4, n, squares);
                let right = a.mul(2, total, total);
                let numerator = a.sub(4, left, right);
                // Divide by n here, then the common sample/population divisor below.
                a.div(4, numerator, n)
            } else {
                let correction = a.mul(2, total, mean);
                a.sub(4, squares, correction)
            }
        }
    };
    let mut value = if function == "DEVSQ" { ss } else { a.div(4, ss, d) };
    if function.starts_with("STDEV") {
        value = a.sqrt(value);
    }
    (value, a.trace)
}

#[derive(Clone, Serialize, Deserialize)]
struct Observation {
    trace_index: usize,
    raw: Value,
    values_bits: Vec<u64>,
}

fn numeric_values(raw: &Value) -> Result<Vec<u64>, &'static str> {
    let mut bits = Vec::new();
    for (i, arg) in raw["args"]
        .as_array()
        .ok_or("argument_shape")?
        .iter()
        .enumerate()
    {
        let v = calc_verdict::decode_argument(arg, &format!("args[{i}]"))
            .map_err(|_| "unsupported_argument")?;
        if let Some(array) = v.as_array() {
            for cell in array.iter_row_major() {
                let CoreRef::Number(x) = cell.core_ref() else {
                    return Err("non_numeric_argument");
                };
                bits.push(x.to_bits());
            }
        } else {
            let CoreRef::Number(x) = v.core_ref() else {
                return Err("non_numeric_argument");
            };
            bits.push(x.to_bits());
        }
    }
    if !raw["result_index"].is_null() {
        return Err("selected_result");
    }
    if !(2..=32).contains(&bits.len()) {
        return Err("length_outside_2_to_32");
    }
    for &b in &bits {
        let x = f64::from_bits(b).abs();
        if !x.is_finite() || (x != 0.0 && !(2f64.powi(-400)..=2f64.powi(400)).contains(&x)) {
            return Err("input_outside_bounded_normal_scope");
        }
    }
    let e = raw["expected"].as_str().ok_or("unsupported_expected")?;
    if e.len() != 18 || !e.starts_with("0x") {
        return Err("non_numeric_expected");
    }
    let x = f64::from_bits(u64::from_str_radix(&e[2..], 16).map_err(|_| "unsupported_expected")?);
    if !x.is_finite() || x.is_subnormal() || (x < 0.0 && raw["function"] != "AVERAGE") {
        return Err("expected_outside_bounded_scope");
    }
    Ok(bits)
}

fn file(path: impl AsRef<Path>) -> BufWriter<fs::File> {
    BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path)
            .unwrap(),
    )
}
fn line(out: &mut impl Write, value: &impl Serialize) {
    serde_json::to_writer(&mut *out, value).unwrap();
    writeln!(out).unwrap();
}

pub fn prepare(input: &str, out: &str) {
    fs::create_dir(out).unwrap();
    let rows = crate::lines::<Value>(input);
    let mut scoped = file(Path::new(out).join("scoped.jsonl"));
    let mut excluded = file(Path::new(out).join("excluded.jsonl"));
    let mut counts = BTreeMap::<String, usize>::new();
    for (trace_index, raw) in rows.iter().enumerate() {
        assert!(matches!(
            raw["function"].as_str().unwrap(),
            "AVERAGE"
                | "DEVSQ"
                | "VAR"
                | "VAR.S"
                | "VAR.P"
                | "VARP"
                | "VARA"
                | "VARPA"
                | "STDEV"
                | "STDEV.S"
                | "STDEV.P"
                | "STDEVP"
                | "STDEVA"
                | "STDEVPA"
        ));
        match numeric_values(raw) {
            Ok(values_bits) => {
                line(
                    &mut scoped,
                    &Observation {
                        trace_index,
                        raw: raw.clone(),
                        values_bits,
                    },
                );
                *counts.entry("scoped".into()).or_default() += 1;
            }
            Err(reason) => {
                line(
                    &mut excluded,
                    &json!({"trace_index":trace_index,"raw":raw,"reason":reason}),
                );
                *counts.entry(reason.into()).or_default() += 1;
            }
        }
    }
    let mut summary = file(Path::new(out).join("summary.json"));
    line(
        &mut summary,
        &json!({"recorded":rows.len(),"counts":counts,"unread":"not represented in decoded traces; retain original judge summary","new_excel_readings":0}),
    );
}

pub fn run<F: Fpu>(input: &str, function: &str, out: &str) {
    run_family::<F>(input, function, out, false)
}

pub fn run_family<F: Fpu>(input: &str, function: &str, out: &str, structural: bool) {
    let start = std::time::Instant::now();
    fs::create_dir(out).unwrap();
    let observations = crate::lines::<Observation>(input);
    let rows: Vec<_> = observations
        .iter()
        .filter(|r| r.raw["function"] == function)
        .collect();
    assert!(!rows.is_empty());
    for r in &rows {
        assert_eq!(numeric_values(&r.raw).unwrap(), r.values_bits);
    }
    let ps = if structural { structural_profiles(function) } else { function_profiles(function) };
    let mut pf = file(Path::new(out).join("profiles.jsonl"));
    for p in &ps {
        line(&mut pf, p);
    }
    pf.flush().unwrap();
    let mut alive: Vec<_> = (0..ps.len()).collect();
    let mut kills = file(Path::new(out).join("first-kills.jsonl"));
    let mut witnesses = file(Path::new(out).join("witnesses.jsonl"));
    let mut executions = 0;
    for (order, row) in rows.iter().enumerate() {
        let before = alive.len();
        alive.retain(|&id| {
            let (value, _) = evaluate::<F>(&ps[id], &row.values_bits, function, false);
            executions += 1;
            let got = eval::outcome::<F>(value);
            if got == row.raw["expected"] {
                true
            } else {
                line(
                    &mut kills,
                    &json!({"candidate_id":id,"order":order,"trace_index":row.trace_index,
                    "root":value.exact_encoding(),"got":got,"expected":row.raw["expected"]}),
                );
                false
            }
        });
        line(
            &mut witnesses,
            &json!({"order":order,"trace_index":row.trace_index,"killed":before-alive.len(),"remaining":alive.len()}),
        );
    }
    let mut surviving = file(Path::new(out).join("survivors.json"));
    line(&mut surviving, &alive);
    let summary = json!({"complete":true,"function":function,"programs":ps.len(),"observations":rows.len(),"executions":executions,
        "survivors":alive.len(),"seconds":start.elapsed().as_secs_f64(),
        "family":if structural {"structural-r17"} else {"r16"},
        "claim":"Exhaustive exact first-contradiction elimination of the declared finite family, not all variance programs or Excel kernels."});
    line(&mut file(Path::new(out).join("summary.json")), &summary);
    println!("{summary}");
}

/// Mathematical controls, plus actual registered-route replay of the unmodified product.
pub fn controls(input: &str, out: &str) {
    use crate::ported::excel_numeric::{Hard, Soft};
    let mut output = file(out);
    let ps = profiles();
    for p in ps.iter().step_by(625) {
        for (function, xs, expected) in [
            ("VAR.S", vec![-3.0, 0.0, 3.0], 9.0_f64),
            ("VAR.P", vec![-3.0, 0.0, 3.0], 6.0),
            ("STDEV.S", vec![-3.0, 0.0, 3.0], 3.0),
            ("DEVSQ", vec![-3.0, 0.0, 3.0], 18.0),
            ("AVERAGE", vec![-3.0, 0.0, 3.0], 0.0),
        ] {
            let bits: Vec<_> = xs.into_iter().map(f64::to_bits).collect();
            let (hard, trace) = evaluate::<Hard>(p, &bits, function, true);
            let (soft, _) = evaluate::<Soft>(p, &bits, function, false);
            assert_eq!(hard.exact_encoding(), soft.exact_encoding());
            assert_eq!(hard.to_f64::<Hard>().to_bits(), expected.to_bits());
            line(
                &mut output,
                &json!({"kind":"mathematical_control","profile":p,"function":function,"values_bits":bits,"root":hard.exact_encoding(),"trace":trace}),
            );
        }
    }
    for raw in crate::lines::<Value>(input) {
        let args = raw["args"]
            .as_array()
            .unwrap()
            .iter()
            .enumerate()
            .map(|(i, a)| calc_verdict::decode_argument(a, &format!("args[{i}]")).unwrap())
            .collect();
        let (got, _) = crate::invoke_control(raw["function"].as_str().unwrap(), args);
        assert_eq!(got, raw["actual"]);
        line(
            &mut output,
            &json!({"kind":"registered_product_replay","raw":raw,"got":got}),
        );
    }
}
