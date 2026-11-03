//! R20 fixed rational Bessel seeds, from SheetJS/bessel (Apache-2.0).
//! Reference revision 2e2331329891dc5c16bda3461db0b9fe590ff2e9, bessel.js.
//! Copyright (C) 2013-present SheetJS. See the pinned source and licence in
//! /home/ubuntu/codex-runs/dna-firehorse-20261103-round20/permitted-bessel-source.
//! Implementer: approved facts and this Tier2 source only; no replaced body/history.
use crate::{dsl::*, eval::*, json_write, lines, ported::excel_numeric::*};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::{
    collections::BTreeMap,
    fs,
    io::{BufWriter, Write},
    path::Path,
};

// Tables retain their published decimal values, in ascending power order.
const J0_NUM: &[f64] = &[
    57568490574.0,
    -13362590354.0,
    651619640.7,
    -11214424.18,
    77392.33017,
    -184.9052456,
];
const J0_DEN: &[f64] = &[
    57568490411.0,
    1029532985.0,
    9494680.718,
    59272.64853,
    267.8532712,
    1.0,
];
const J1_NUM: &[f64] = &[
    72362614232.0,
    -7895059235.0,
    242396853.1,
    -2972611.439,
    15704.48260,
    -30.16036606,
];
const J1_DEN: &[f64] = &[
    144725228442.0,
    2300535178.0,
    18583304.74,
    99447.43394,
    376.9991397,
    1.0,
];
const Y0_NUM: &[f64] = &[
    -2957821389.0,
    7062834065.0,
    -512359803.6,
    10879881.29,
    -86327.92757,
    228.4622733,
];
const Y0_DEN: &[f64] = &[
    40076544269.0,
    745249964.8,
    7189466.438,
    47447.26470,
    226.1030244,
    1.0,
];
const Y1_NUM: &[f64] = &[
    -0.4900604943e13,
    0.1275274390e13,
    -0.5153438139e11,
    0.7349264551e9,
    -0.4237922726e7,
    0.8511937935e4,
];
const Y1_DEN: &[f64] = &[
    0.2499580570e14,
    0.4244419664e12,
    0.3733650367e10,
    0.2245904002e8,
    0.1020426050e6,
    0.3549632885e3,
    1.0,
];
const AMP0: &[f64] = &[
    1.0,
    -0.1098628627e-2,
    0.2734510407e-4,
    -0.2073370639e-5,
    0.2093887211e-6,
];
const PHASE_J0: &[f64] = &[
    -0.1562499995e-1,
    0.1430488765e-3,
    -0.6911147651e-5,
    0.7621095161e-6,
    -0.934935152e-7,
];
// The published Y0 table has a different last decimal from J0; retain it.
const PHASE_Y0: &[f64] = &[
    -0.1562499995e-1,
    0.1430488765e-3,
    -0.6911147651e-5,
    0.7621095161e-6,
    -0.934945152e-7,
];
const AMP1: &[f64] = &[
    1.0,
    0.183105e-2,
    -0.3516396496e-4,
    0.2457520174e-5,
    -0.240337019e-6,
];
const PHASE1: &[f64] = &[
    0.04687499995,
    -0.2002690873e-3,
    0.8449199096e-5,
    -0.88228987e-6,
    0.105787412e-6,
];
const W: f64 = 0.636619772;

/// The pinned reference's complete adaptive backward program, without tuning.
/// The application admits integer orders first; the selected route has n >= 2.
pub fn published_j_backward(order: i32, x: f64) -> f64 {
    if x == 0.0 {
        return 0.0;
    }
    let magnitude = x.abs();
    let n = f64::from(order);
    let start = 2.0 * ((n + (40.0 * n).sqrt().floor()) / 2.0).floor();
    let twice_over_x = 2.0 / magnitude;
    let mut result = 0.0;
    let mut current = 1.0;
    let mut next = 0.0;
    let mut total = 0.0;
    let mut add_to_total = false;
    for index in (1..=start as i32).rev() {
        let previous = f64::from(index) * twice_over_x * current - next;
        next = current;
        current = previous;
        if current.abs() > 1E10 {
            current *= 1E-10;
            next *= 1E-10;
            result *= 1E-10;
            total *= 1E-10;
        }
        if add_to_total {
            total += current;
        }
        add_to_total = !add_to_total;
        if index == order {
            result = next;
        }
    }
    total = 2.0 * total - current;
    result /= total;
    if x < 0.0 && order % 2 != 0 {
        -result
    } else {
        result
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Program {
    pub id: usize,
    pub polynomial: EvalModel,
    pub reconstruction: EvalModel,
    pub assignment_stores: bool,
    pub helpers: usize,     // 0 worksheet, 1 cleared CRT, 2 libm
    pub association: usize, // three fixed algebraic forms, applied uniformly
    pub reduction: usize,   // 64/(x*x) or (8/x)*(8/x)
    #[serde(default)]
    pub group_y0_log: bool, // one fixed Y0 expression: J0*(W*log(x))
}

pub fn family() -> Vec<Program> {
    let mut out = Vec::new();
    for polynomial in WIDE_MODELS {
        for reconstruction in WIDE_MODELS {
            for assignment_stores in [false, true] {
                for helpers in 0..3 {
                    for association in 0..3 {
                        for reduction in 0..2 {
                            out.push(Program {
                                id: out.len(),
                                polynomial,
                                reconstruction,
                                assignment_stores,
                                helpers,
                                association,
                                reduction,
                                group_y0_log: false,
                            });
                        }
                    }
                }
            }
        }
    }
    assert_eq!(out.len(), 900);
    out
}

struct Machine<'a, F: Fpu> {
    p: &'a Program,
    marker: std::marker::PhantomData<F>,
}
impl<'a, F: Fpu> Machine<'a, F> {
    fn new(p: &'a Program) -> Self {
        Self {
            p,
            marker: std::marker::PhantomData,
        }
    }
    fn add(&self, a: Val, b: Val, m: EvalModel) -> Val {
        binop::<F>(Op::Add(0, 1), m, a, b)
    }
    fn mul(&self, a: Val, b: Val, m: EvalModel) -> Val {
        binop::<F>(Op::Mul(0, 1), m, a, b)
    }
    fn div(&self, a: Val, b: Val, m: EvalModel) -> Val {
        binop::<F>(Op::Div(0, 1), m, a, b)
    }
    fn neg(&self, a: Val) -> Val {
        match a {
            Val::F64(b) => Val::F64(b ^ (1 << 63)),
            Val::Ext(mut e) => {
                e.0[9] ^= 0x80;
                Val::Ext(e)
            }
        }
    }
    fn sub(&self, a: Val, b: Val, m: EvalModel) -> Val {
        self.add(a, self.neg(b), m)
    }
    fn assign(&self, a: Val) -> Val {
        if self.p.assignment_stores {
            Val::f64(a.to_f64::<F>())
        } else {
            a
        }
    }
    fn poly(&self, y: Val, c: &[f64]) -> Val {
        let m = self.p.polynomial;
        let mut v = Val::f64(*c.last().unwrap());
        for &a in c[..c.len() - 1].iter().rev() {
            v = self.add(self.mul(y, v, m), Val::f64(a), m);
        }
        self.assign(v)
    }
    fn sqrt(&self, a: Val) -> Val {
        let m = self.p.reconstruction;
        if m == EvalModel::Strict {
            return Val::f64(a.to_f64::<F>().sqrt());
        }
        let cw = if matches!(m, EvalModel::StoredPc53 | EvalModel::ResidentPc53) {
            (CW & !0x0300) | 0x0200
        } else {
            CW
        };
        let v = Val::Ext(F::sqrt(&a.to_ext(), cw));
        if matches!(m, EvalModel::StoredPc64 | EvalModel::StoredPc53) {
            Val::f64(v.to_f64::<F>())
        } else {
            v
        }
    }
    fn elementary(&self, kind: &str, a: Val) -> Val {
        let x = a.to_f64::<F>(); // scalar elementary call boundary
        Val::f64(match (kind, self.p.helpers) {
            ("sin", 0) => excel_sin(x),
            ("cos", 0) => excel_cos(x),
            ("log", 0) => excel_log(x),
            ("sin", 1) => crt_sin(x),
            ("cos", 1) => crt_cos(x),
            ("log", 1) => crt_log(x),
            ("sin", 2) => libm::sin(x),
            ("cos", 2) => libm::cos(x),
            ("log", 2) => libm::log(x),
            _ => unreachable!(),
        })
    }
    fn seed(&self, second: bool, order: usize, x: f64) -> Val {
        let x = Val::f64(x);
        let pm = self.p.polynomial;
        let m = self.p.reconstruction;
        let y = self.assign(self.mul(x, x, pm));
        if x.to_f64::<F>() < 8.0 {
            let (num, den) = match (second, order) {
                (false, 0) => (J0_NUM, J0_DEN),
                (false, 1) => (J1_NUM, J1_DEN),
                (true, 0) => (Y0_NUM, Y0_DEN),
                (true, 1) => (Y1_NUM, Y1_DEN),
                _ => unreachable!(),
            };
            let a = self.poly(y, num);
            let b = self.poly(y, den);
            let rational = if order == 0 {
                self.div(a, b, m)
            } else {
                match self.p.association {
                    0 => self.div(self.assign(self.mul(x, a, m)), b, m),
                    1 => self.mul(x, self.div(a, b, m), m),
                    2 => self.mul(self.div(x, b, m), a, m),
                    _ => unreachable!(),
                }
            };
            if !second {
                return rational;
            }
            let j = Val::f64(self.seed(false, order, x.to_f64::<F>()).to_f64::<F>());
            let log = self.elementary("log", x);
            let w = Val::f64(W);
            let term = if order == 0 {
                match if self.p.group_y0_log {
                    1
                } else {
                    self.p.association
                } {
                    0 => self.mul(self.mul(w, j, m), log, m),
                    1 => self.mul(j, self.mul(w, log, m), m),
                    2 => self.mul(self.mul(j, log, m), w, m),
                    _ => unreachable!(),
                }
            } else {
                let inv = self.div(Val::f64(1.0), x, m);
                match self.p.association {
                    0 => self.mul(w, self.sub(self.mul(j, log, m), inv, m), m),
                    1 => self.sub(self.mul(w, self.mul(j, log, m), m), self.mul(w, inv, m), m),
                    2 => self.sub(self.mul(self.mul(w, j, m), log, m), self.mul(w, inv, m), m),
                    _ => unreachable!(),
                }
            };
            self.add(rational, term, m)
        } else {
            let phase = self.assign(self.sub(
                x,
                Val::f64(if order == 0 { 0.785398164 } else { 2.356194491 }),
                m,
            ));
            let y = if self.p.reduction == 0 {
                self.div(Val::f64(64.0), y, pm)
            } else {
                let z = self.div(Val::f64(8.0), x, pm);
                self.mul(z, z, pm)
            };
            let y = self.assign(y);
            let a = self.poly(y, if order == 0 { AMP0 } else { AMP1 });
            let b = self.poly(
                y,
                if order == 0 {
                    if second { PHASE_Y0 } else { PHASE_J0 }
                } else {
                    PHASE1
                },
            );
            let sine = self.elementary("sin", phase);
            let cosine = self.elementary("cos", phase);
            let (leading, trailing) = if second {
                (sine, cosine)
            } else {
                (cosine, sine)
            };
            let first = self.mul(leading, a, m);
            let tail = match self.p.association {
                0 => self.div(self.mul(self.mul(trailing, b, m), Val::f64(8.0), m), x, m),
                1 => self.mul(self.mul(trailing, b, m), self.div(Val::f64(8.0), x, m), m),
                2 => self.mul(trailing, self.mul(b, self.div(Val::f64(8.0), x, m), m), m),
                _ => unreachable!(),
            };
            let combined = if second {
                self.add(first, tail, m)
            } else {
                self.sub(first, tail, m)
            };
            self.mul(self.sqrt(self.div(Val::f64(W), x, m)), combined, m)
        }
    }
}

fn decode(s: &Value) -> Option<f64> {
    let s = s.as_str()?;
    let h = s.strip_prefix("0x")?;
    Some(f64::from_bits(u64::from_str_radix(h, 16).ok()?))
}
fn prediction<F: Fpu>(p: &Program, r: &Value) -> (String, String) {
    let f = r["function"].as_str().unwrap();
    let Some(a) = r["args"].as_array() else {
        return ("refusal:input-shape".into(), "unsupported".into());
    };
    if a.len() != 2 {
        return ("refusal:input-shape".into(), "unsupported".into());
    }
    let (Some(x), Some(n)) = (decode(&a[0]), decode(&a[1])) else {
        return ("refusal:input-type".into(), "unsupported".into());
    };
    if !x.is_finite()
        || !n.is_finite()
        || n < 0.0
        || n.trunc() > 1.0
        || (f != "BESSELJ" && f != "BESSELY")
    {
        return ("refusal:outside-seed-scope".into(), "unsupported".into());
    }
    let second = f == "BESSELY";
    let region = format!(
        "{f}:{}:{}",
        n.trunc(),
        if x.abs() < 8.0 { "low" } else { "high" }
    );
    if second && x <= 0.0 {
        return ("error:Num".into(), region);
    }
    let mut v = Machine::<F>::new(p)
        .seed(second, n.trunc() as usize, x.abs())
        .to_f64::<F>();
    if !second && x < 0.0 && n.trunc() == 1.0 {
        v = -v;
    }
    let result = if v.is_finite() {
        format!("0x{:016x}", v.to_bits())
    } else {
        "error:Num".into()
    };
    (result, region)
}

pub fn run<F: Fpu>(input: &str, out: &str) {
    run_programs::<F>(input, out, family(), false);
}

pub fn replay_fixed<F: Fpu>(input: &str, out: &str) {
    let mut p = family().remove(0);
    p.id = 900;
    p.group_y0_log = true;
    run_programs::<F>(input, out, vec![p], false);
}

pub fn replay_followup<F: Fpu>(input: &str, out: &str) {
    let mut p = family().remove(0);
    p.id = 1000;
    run_programs::<F>(input, out, vec![p], true);
}

fn prediction_followup<F: Fpu>(p: &Program, r: &Value) -> (String, String) {
    let Some(args) = r["args"].as_array() else {
        return ("refusal:input-shape".into(), "unsupported".into());
    };
    if args.len() != 2 {
        return ("refusal:input-shape".into(), "unsupported".into());
    }
    let (Some(x), Some(raw_order)) = (decode(&args[0]), decode(&args[1])) else {
        return ("refusal:input-type".into(), "unsupported".into());
    };
    let n = raw_order.trunc();
    if !x.is_finite() || !raw_order.is_finite() || raw_order < 0.0 || n > 10_000_000.0 {
        return ("refusal:outside-admission".into(), "unsupported".into());
    }
    let function = r["function"].as_str().unwrap();
    let (value, region) = if function == "BESSELJ" && n >= 2.0 && x != 0.0 && x.abs() <= n {
        (
            published_j_backward(n as i32, x),
            format!(
                "BESSELJ:backward:{}",
                if x.abs() < 8.0 { "low" } else { "high" }
            ),
        )
    } else if function == "BESSELY" && n == 1.0 && x >= 8.0 {
        (
            Machine::<F>::new(p).seed(true, 1, x).to_f64::<F>(),
            "BESSELY:1:high".into(),
        )
    } else {
        return (
            "refusal:outside-fixed-followup-scope".into(),
            "unsupported".into(),
        );
    };
    let prediction = if value.is_finite() {
        format!("0x{:016x}", value.to_bits())
    } else {
        "error:Num".into()
    };
    (prediction, region)
}

fn run_programs<F: Fpu>(input: &str, out: &str, programs: Vec<Program>, followup: bool) {
    fs::create_dir_all(out).unwrap();
    let rows = lines::<Value>(input);
    json_write(Path::new(out).join("programs.json"), &programs);
    let mut ranking = Vec::new();
    let mut region_rankings: BTreeMap<String, Vec<Value>> = BTreeMap::new();
    for p in &programs {
        let mut regions: BTreeMap<String, [usize; 5]> = BTreeMap::new();
        let mut total = [0usize; 5];
        for r in &rows {
            let (actual, region) = if followup {
                prediction_followup::<F>(p, r)
            } else {
                prediction::<F>(p, r)
            };
            let expected = r["expected"].as_str().unwrap();
            let baseline = r["actual"].as_str().unwrap();
            let matched = actual == expected;
            let old = baseline == expected;
            for t in [&mut total, regions.entry(region).or_default()] {
                t[0] += 1;
                t[1] += usize::from(matched);
                t[2] += usize::from(old);
                t[3] += usize::from(matched && !old);
                t[4] += usize::from(!matched && old);
            }
        }
        for (region, c) in &regions {
            region_rankings
                .entry(region.clone())
                .or_default()
                .push(json!({"id":p.id,"counts":c}));
        }
        ranking.push(json!({"program":p,"counts":total,"regions":regions}));
    }
    let order = |a: &Value, b: &Value| {
        b["counts"][1]
            .as_u64()
            .cmp(&a["counts"][1].as_u64())
            .then_with(|| a["counts"][4].as_u64().cmp(&b["counts"][4].as_u64()))
    };
    ranking.sort_by(order);
    for ranks in region_rankings.values_mut() {
        ranks.sort_by(order);
        ranks.truncate(10);
    }
    json_write(Path::new(out).join("ranking.json"), &ranking);
    json_write(Path::new(out).join("region-top.json"), &region_rankings);
    let top = &ranking[0];
    let p: Program = serde_json::from_value(top["program"].clone()).unwrap();
    let mut f =
        BufWriter::new(fs::File::create(Path::new(out).join("best-predictions.jsonl")).unwrap());
    for r in &rows {
        let (actual, region) = if followup {
            prediction_followup::<F>(&p, r)
        } else {
            prediction::<F>(&p, r)
        };
        serde_json::to_writer(&mut f,&json!({"raw":r,"program":p.id,"prediction":actual,"region":region,"matched":actual==r["expected"].as_str().unwrap()})).unwrap();
        writeln!(f).unwrap();
    }
    json_write(
        Path::new(out).join("summary.json"),
        &json!({
            "input":input,"rows":rows.len(),"programs":programs.len(),"count_order":["all","exact","baseline_exact","new_exact","lost_exact"],
            "scope":if followup {"Exact published backward J for admitted n>=2 and0<abs(x)<=n; complete Y1-high survivor. All other rows explicitly refused and counted fail in research totals."} else {"J0/J1/Y0/Y1; every unsupported and numeric/type mismatch remains a failure in overall counts"},
            "best":top,"region_top":region_rankings,
            "source":"SheetJS/bessel 2e2331329891dc5c16bda3461db0b9fe590ff2e9 (Apache-2.0), approved R20 facts; IEEE754/Intel x87; original unchanged Excel readings",
            "discovery":"All captured readings are discovery. No held-out or product parity claim."
        }),
    );
    println!(
        "bessel-rational {} rows {} programs best {} exact, {} lost baseline matches",
        rows.len(),
        programs.len(),
        top["counts"][1],
        top["counts"][4]
    );
}
