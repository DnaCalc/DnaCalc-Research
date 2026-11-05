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

/// One fixed high-J seed tree; the two literal edits are declared in the unit.
fn high_j_literal_value<F:Fpu>(x:f64,order:usize,program:usize)->f64 {
    high_seed_grouped_value::<F>(x,order,false,program,false,0)
}

fn high_seed_grouped_value<F:Fpu>(x:f64,order:usize,second:bool,program:usize,z_square:bool,association:usize)->f64 {
    let mut p=family().remove(0);p.group_y0_log=true;
    let machine=Machine::<F>::new(&p);let model=EvalModel::Strict;
    let edited_amp0=[1.0,-0.1098628267e-2,0.2734510407e-4,-0.2073370639e-5,0.2093887211e-6];
    let repeated_amp1=[AMP1[0],AMP1[1],AMP1[2],AMP1[2],AMP1[3],AMP1[4]];
    let argument=Val::f64(x);
    let square=machine.mul(argument,argument,model);
    let z=machine.div(Val::f64(8.0),argument,model);
    let reduction=if z_square{machine.mul(z,z,model)}else{machine.div(Val::f64(64.0),square,model)};
    let phase=machine.sub(argument,Val::f64(if order==0{0.785398164}else{2.356194491}),model);
    let sine=machine.elementary("sin",phase);let cosine=machine.elementary("cos",phase);
    let scale=machine.sqrt(machine.div(Val::f64(W),argument,model));
    let amplitude=machine.poly(reduction,if order==0{if program==0||second{AMP0}else{&edited_amp0}}
        else if program==0||second{AMP1}else{&repeated_amp1});
    let correction=machine.poly(reduction,if order==0{if program==2||second{PHASE_Y0}else{PHASE_J0}}else{PHASE1});
    let (leading,trailing)=if second{(sine,cosine)}else{(cosine,sine)};
    let leading=machine.mul(leading,amplitude,model);
    let trailing=match association{
        0=>machine.div(machine.mul(machine.mul(trailing,correction,model),Val::f64(8.0),model),argument,model),
        1=>machine.mul(machine.mul(z,trailing,model),correction,model),
        2=>machine.mul(trailing,machine.mul(z,correction,model),model),
        3=>machine.mul(machine.mul(trailing,correction,model),z,model),
        _=>unreachable!(),
    };
    let combined=if second{machine.add(leading,trailing,model)}else{machine.sub(leading,trailing,model)};
    machine.mul(scale,combined,model).to_f64::<F>()
}

/// Simple left-to-right z*trig*polynomial comes first. The existing
/// quotient-last reconstruction is a control, never an accuracy-ranked choice.
pub fn replay_high_seed_groups<F:Fpu>(input:&str,out:&str){
    let destination=Path::new(out);fs::create_dir(destination).unwrap();let rows=lines::<Value>(input);
    let profiles=[(true,1),(false,1),(true,2),(false,2),(true,3),(false,3),(true,0),(false,0)];
    let mut counts=vec![[0_usize;5];profiles.len()];let mut excluded=0;let mut predictions=Vec::new();
    for row in &rows{
        let name=row["function"].as_str().unwrap();let second=name=="BESSELY";
        if !second&&name!="BESSELJ"{excluded+=1;continue;}
        let args=row["args"].as_array().unwrap();let(Some(x),Some(n))=(decode(&args[0]),decode(&args[1]))else{excluded+=1;continue;};
        let order=n.trunc();
        if !x.is_finite()||x<8.0||n<0.0||order>1.0||!row["actual"].as_str().unwrap().starts_with("0x"){excluded+=1;continue;}
        let baseline=crate::invoke_control(name,vec![crate::value::CalcValue::number(x),crate::value::CalcValue::number(n)]).0;
        assert_eq!(baseline,row["actual"],"frozen high-seed baseline {}",row["id"]);
        let values:Vec<_>=profiles.iter().enumerate().map(|(i,&(reduction,association))|{
            let value=high_seed_grouped_value::<F>(x,order as usize,second,1,reduction,association);
            let prediction=format!("0x{:016x}",value.to_bits());let old=baseline==row["expected"];let now=prediction==row["expected"];
            for(j,b)in [true,now,old,now&&!old,old&&!now].into_iter().enumerate(){counts[i][j]+=usize::from(b);}
            prediction
        }).collect();predictions.push(json!({"original":row,"predictions":values}));
    }
    json_write(destination.join("predictions.json"),&predictions);
    let result=json!({"complete":true,"input_rows":rows.len(),"excluded_rows":excluded,"profiles":profiles,
        "count_order":["all","exact","baseline_exact","gains","losses"],"counts":counts,
        "candidate_order":"Named z=8/x, then z*trig*correction first; z-square then64/(x*x). Other fixed associations follow. The pinned SheetJS quotient-last tree is the last control. No match-count ranking or fitted coefficient. The first draft incorrectly called the z-left alternative published; it is a simple algebraic reconstruction hypothesis, not SheetJS's grouping.",
        "limits":"High positive admitted seed orders0/1 only, original recorded J and Y observations. J amplitude edits are held fixed, Y literals unchanged. Strict arithmetic and compiled worksheet SIN/COS leaves. Every contradiction retained; complete recurrence/application effect remains separate."});
    json_write(destination.join("summary.json"),&result);println!("{result}");
}

/// The existing caller owns admission; this helper adds only the original parity.
pub(crate) fn inferred_high_j_seed(order:usize,x:f64)->f64 {
    let value=high_j_literal_value::<Native>(x.abs(),order,1);
    if x<0.0&&order==1{-value}else{value}
}

/// Two handwritten transcription hypotheses, in fixed code-shape order.
pub fn replay_high_j_literals<F:Fpu>(input:&str,out:&str) {
    let destination=Path::new(out);fs::create_dir(destination).unwrap();let rows=lines::<Value>(input);
    let mut counts=[[0_usize;5];3];let mut predictions=Vec::new();let mut excluded=0;
    for row in &rows {
        let Some(args)=row["args"].as_array()else{excluded+=1;continue;};
        if row["function"]!="BESSELJ"||args.len()!=2||!row["result_index"].is_null(){excluded+=1;continue;}
        let(Some(x),Some(n))=(decode(&args[0]),decode(&args[1]))else{excluded+=1;continue;};
        if !(x>=8.0&&x<=134217727.0&&(n==0.0||n==1.0))||!row["expected"].as_str().unwrap().starts_with("0x"){excluded+=1;continue;}
        let mut values=Vec::new();
        for program in 0..3 {
            let value=high_j_literal_value::<F>(x,n as usize,program);
            let value=format!("0x{:016x}",value.to_bits());
            let now=value==row["expected"];let old=row["actual"]==row["expected"];
            for(j,b)in [true,now,old,now&&!old,old&&!now].into_iter().enumerate(){counts[program][j]+=usize::from(b);}
            values.push(value);
        }
        predictions.push(json!({"original":row,"predictions":values}));
    }
    json_write(destination.join("predictions.json"),&predictions);
    let summary=json!({"complete":true,"input_rows":rows.len(),"excluded_rows":excluded,
        "program_order":["published-literal control in fixed study graph","J0 digit transposition and J1 repeated term","same edits with published shared Y0 phase table"],
        "count_order":["all","exact","baseline_exact","gains","losses"],"counts":counts,
        "sources":"Pinned Apache2 SheetJS/bessel and original Excel phase-window observations. J0 decimal transposition and J1 repeated published coefficient are fixed code-shape hypotheses. No numeric coefficient fit or match-count ranking.",
        "limits":"Positive admitted x>=8, seed orders0/1 only. Strict arithmetic and current worksheet SIN/COS are conditional leaves selected by the compiled feature, not the CLI backend selector. Other recurrence/admission behavior is outside this seed study. Excluded rows retain their original outcomes and are not scored here. Every contradiction is retained. Discovery, not held-out."});
    json_write(destination.join("summary.json"),&summary);println!("{summary}");
}

#[cfg(bessel_high_j_audit)]
pub fn replay_high_j_route(input:&str,out:&str) {
    replay_bessel_component(input,out,false,|_,x,n|crate::ported::bessel_high_j::besselj_kernel(x,n));
}

#[cfg(bessel_high_groups_audit)]
pub fn replay_high_seed_route(input:&str,out:&str) {
    replay_bessel_component(input,out,true,|name,x,n|if name=="BESSELJ"{
        crate::ported::bessel_grouped::besselj_kernel(x,n)
    }else{crate::ported::bessel_grouped::bessely_kernel(x,n)});
}

#[cfg(any(bessel_high_j_audit,bessel_high_groups_audit))]
fn replay_bessel_component(input:&str,out:&str,include_y:bool,
    calculate:impl Fn(&str,f64,f64)->Result<f64,calc_func::value::WorksheetErrorCode>) {
    let destination=Path::new(out);fs::create_dir(destination).unwrap();let rows=lines::<Value>(input);
    let mut counts=[0_usize;5];let mut changes=Vec::new();let mut excluded=0;
    let encode=|v:Result<f64,calc_func::value::WorksheetErrorCode>|v.map(|x|format!("0x{:016x}",x.to_bits())).unwrap_or_else(|e|format!("error:{e:?}"));
    for row in &rows {
        let Some(args)=row["args"].as_array()else{excluded+=1;continue;};
        let name=row["function"].as_str().unwrap();
        if !(name=="BESSELJ"||(include_y&&name=="BESSELY"))||args.len()!=2||!row["result_index"].is_null(){excluded+=1;continue;}
        let(Some(x),Some(n))=(decode(&args[0]),decode(&args[1]))else{excluded+=1;continue;};
        let baseline=encode(if name=="BESSELJ"{crate::ported::bessel::besselj_kernel(x,n)}else{crate::ported::bessel::bessely_kernel(x,n)});
        assert_eq!(baseline,row["actual"],"current component and registered original {}",row["id"]);
        let prediction=encode(calculate(name,x,n));
        let old=baseline==row["expected"];let now=prediction==row["expected"];
        for(j,b)in [true,now,old,now&&!old,old&&!now].into_iter().enumerate(){counts[j]+=usize::from(b);}
        if prediction!=baseline{changes.push(json!({"original":row,"prediction":prediction}));}
    }
    assert_eq!(counts[0]+excluded,rows.len());
    json_write(destination.join("changes.json"),&changes);
    let result=json!({"complete":true,"input_rows":rows.len(),"excluded_rows":excluded,"counts":counts,
        "count_order":["all","exact","baseline_exact","gains","losses"],
        "variant":if include_y{"high J/Y correction uses ((8/x)*trig)*polynomial"}else{"two high-J seed calls replaced"},
        "limits":"Whole current cleared Bessel component with only the declared high-region seed change. Original admission, parity, integer-order handling, forward/backward recurrence and publication retained. Every original registered baseline is checked. This is a component study; the application route still requires before/after replay. No coefficient, cutoff, arithmetic mask or row selector fit to match counts."});
    json_write(destination.join("summary.json"),&result);println!("{result}");
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

/// Three fixed algebraic realizations of DLMF10.6.1, now across five coherent
/// arithmetic models and optional binary64 iteration stores. Seeds and coercion use the current invoke
/// route; this diagnostic never substitutes its result into that route.
pub fn replay_y_forward<F:Fpu>(input: &str, out: &str) {
    use calc_func::value::{CalcValue, CoreValue};
    let destination = Path::new(out);
    fs::create_dir(destination).unwrap();
    let rows = lines::<Value>(input);
    let mut profiles=Vec::new();for shape in 0..3 {for model in WIDE_MODELS {for iteration_store in [false,true] {profiles.push((shape,model,iteration_store));}}}
    let mut counts = vec![[0_usize;5];profiles.len()];
    let mut excluded = BTreeMap::<String, usize>::new();
    let mut file =
        BufWriter::new(fs::File::create_new(destination.join("predictions.jsonl")).unwrap());
    for (index, row) in rows.iter().enumerate() {
        if row["function"] != "BESSELY" {
            *excluded.entry("other_function".into()).or_default() += 1;
            continue;
        }
        let a = row["args"].as_array().unwrap();
        let (Some(x), Some(raw_order)) = (decode(&a[0]), decode(&a[1])) else {
            *excluded.entry("input_type".into()).or_default() += 1;
            continue;
        };
        if !(x.is_finite()
            && raw_order.is_finite()
            && raw_order >= 0.0
            && raw_order.trunc() <= 10_000_000.0)
        {
            *excluded
                .entry("outside_order_admission".into())
                .or_default() += 1;
            continue;
        }
        let args = vec![CalcValue::number(x), CalcValue::number(raw_order)];
        let (baseline, _) = crate::invoke_control("BESSELY", args);
        let expected = row["expected"].as_str().unwrap();
        let seeds: Vec<_> = [0.0, 1.0]
            .map(|n| {
                crate::invoke_control("BESSELY", vec![CalcValue::number(x), CalcValue::number(n)])
            })
            .into();
        for (program,&(shape,model,iteration_store)) in profiles.iter().enumerate() {
            let order = raw_order.trunc() as usize;
            let prediction = if order <= 1 {
                seeds[order].0.clone()
            } else if let (Some(previous), Some(current)) = (seeds[0].1, seeds[1].1) {
                let mut previous=Val::f64(previous);let mut current=Val::f64(current);
                let mul=|a,b|binop::<F>(Op::Mul(0,1),model,a,b);
                let div=|a,b|binop::<F>(Op::Div(0,1),model,a,b);
                let sub=|a,b|binop::<F>(Op::Add(0,1),model,a,crate::boundary::neg(b));
                let twice_over_x = div(Val::f64(2.0),Val::f64(x));
                for j in 1..order {
                    let j = Val::f64(j as f64);
                    let next = match shape {
                        0 => sub(mul(mul(current,j),twice_over_x),previous),
                        1 => sub(mul(mul(j,twice_over_x),current),previous),
                        2 => sub(mul(div(mul(Val::f64(2.0),j),Val::f64(x)),current),previous),
                        _ => unreachable!(),
                    };
                    previous = current;
                    current = if iteration_store {Val::f64(next.to_f64::<F>())}else{next};
                }
                let current=current.to_f64::<F>();
                if !current.is_finite() {
                    "error:Num".into()
                } else {
                    calc_verdict::encode_result(&CalcValue::number(if current.is_subnormal() {
                        0.0
                    } else {
                        current
                    }))
                    .unwrap()
                }
            } else {
                seeds.iter().find(|s| s.1.is_none()).unwrap().0.clone()
            };
            if shape==1&&model==EvalModel::Strict {assert_eq!(prediction,baseline,"strict Y forward route {}",row["id"]);}
            let matched = prediction == expected;
            let was = baseline == expected;
            for (slot, increment) in counts[program].iter_mut().zip([
                1,
                usize::from(matched),
                usize::from(was),
                usize::from(matched && !was),
                usize::from(!matched && was),
            ]) {
                *slot += increment;
            }
            writeln!(file, "{}", json!({"trace_index":index,"id":row["id"],"program":program,"prediction":prediction,"expected":expected,"baseline":baseline,"matched":matched})).unwrap();
        }
    }
    let result = json!({"complete":true,"input_rows":rows.len(),"profiles":profiles,"counts_all_exact_baseline_gains_losses":counts,"excluded":excluded,
        "source":"SheetJS/bessel2e2331329891dc5c16bda3461db0b9fe590ff2e9 bessel.js Apache-2.0; DLMF10.6.1",
        "limits":"Discovery only. Every original scoped outcome is compared, including domain errors; unsupported records remain obligations. Current invoke seeds are conditional leaves. No coefficient, cutoff or store boundary fitted."});
    json_write(destination.join("summary.json"), &result);
    println!("{result}");
}

// PhpSpreadsheet BesselI.php, MIT, copyright (c)2019-2026
// PhpSpreadsheet Authors, revision60cc826911f1a1a125d756f22f0476a88504ab1f.
// See notices/PhpSpreadsheet.txt. This is the complete fixed published
// numerical program; only EXP is supplied by the existing worksheet chain.
const I0_LOW: &[f64] = &[
    1.0,
    3.5156229,
    3.0899424,
    1.2067492,
    0.2659732,
    0.360768e-1,
    0.45813e-2,
];
const I0_HIGH: &[f64] = &[
    0.39894228,
    0.1328592e-1,
    0.225319e-2,
    -0.157565e-2,
    0.916281e-2,
    -0.2057706e-1,
    0.2635537e-1,
    -0.1647633e-1,
    0.392377e-2,
];
const I1_LOW: &[f64] = &[
    0.5,
    0.87890594,
    0.51498869,
    0.15084934,
    0.2658733e-1,
    0.301532e-2,
    0.32411e-3,
];
const I1_HIGH: &[f64] = &[
    0.39894228,
    -0.3988024e-1,
    -0.362018e-2,
    0.163801e-2,
    -0.1031555e-1,
    0.2282967e-1,
    -0.2895312e-1,
    0.1787654e-1,
    -0.420059e-2,
];

fn polynomial(x: f64, coefficients: &[f64]) -> f64 {
    let strict = Program {
        id: 0,
        polynomial: EvalModel::Strict,
        reconstruction: EvalModel::Strict,
        assignment_stores: false,
        helpers: 0,
        association: 0,
        reduction: 0,
        group_y0_log: false,
    };
    Machine::<Hard>::new(&strict)
        .poly(Val::f64(x), coefficients)
        .to_f64::<Hard>()
}

fn published_i0(x: f64) -> f64 {
    let a = x.abs();
    if a < 3.75 {
        let y = x / 3.75;
        polynomial(y * y, I0_LOW)
    } else {
        (excel_exp(a) / a.sqrt()) * polynomial(3.75 / a, I0_HIGH)
    }
}

// Copyright NumFOCUS, Apache-2.0. ITK v5.4.0, commit311b7060ef39e371f3cd209ec135284ff5fde735,
// Modules/Core/Common/include/itkGaussianOperator.hxx:ModifiedBesselI.
// Rust scalar port; the published I0 polynomial uses the existing worksheet EXP.
fn published_i_itk(x: f64, order: usize) -> f64 {
    assert!(order >= 2);
    if x == 0.0 { return 0.0; }
    let twice_over_x = 2.0 / x.abs();
    let mut next = 0.0;
    let mut current = 1.0;
    let mut answer = 0.0;
    let start = 2 * (order + (40.0 * order as f64).sqrt() as usize);
    for j in (1..=start).rev() {
        let previous = next + (j as f64 * twice_over_x) * current;
        next = current;
        current = previous;
        if current.abs() > 1e10 {
            answer *= 1e-10;
            current *= 1e-10;
            next *= 1e-10;
        }
        if j == order { answer = next; }
    }
    answer *= published_i0(x) / current;
    if x < 0.0 && order % 2 == 1 { -answer } else { answer }
}

pub fn replay_i_itk(input: &str, out: &str) {
    use calc_func::value::CalcValue;
    let dest = Path::new(out);
    fs::create_dir(dest).unwrap();
    let rows = lines::<Value>(input);
    let mut counts = BTreeMap::<String, [usize; 5]>::new();
    let mut file = BufWriter::new(fs::File::create_new(dest.join("predictions.jsonl")).unwrap());
    for (index, record) in rows.iter().enumerate() {
        let row = &record["raw"];
        assert_eq!(row["function"], "BESSELI");
        let (raw_x, raw_n) = (decode(&row["args"][0]).unwrap(), decode(&row["args"][1]).unwrap());
        let (baseline, _) = crate::invoke_control("BESSELI", vec![CalcValue::number(raw_x), CalcValue::number(raw_n)]);
        let x = if raw_x.is_subnormal() { 0.0 } else { raw_x };
        let admitted = x.is_finite() && x.abs() <= 709.7 && raw_n.is_finite()
            && raw_n >= 2.0 && raw_n.trunc() <= 10_000_000.0;
        let (prediction, region) = if admitted {
            let (_, seed) = crate::invoke_control("BESSELI", vec![CalcValue::number(x), CalcValue::number(0.0)]);
            assert_eq!(seed.unwrap().to_bits(), published_i0(x).to_bits(), "ITK I0/application seed differs");
            let value = published_i_itk(x, raw_n.trunc() as usize);
            let encoded = if value.is_finite() {
                calc_verdict::encode_result(&CalcValue::number(if value.is_subnormal() { 0.0 } else { value })).unwrap()
            } else { "error:Num".into() };
            (encoded, "published_higher_I")
        } else { (baseline.clone(), "unchanged_outer_or_lower_order") };
        let expected = row["expected"].as_str().unwrap();
        let matched = prediction == expected;
        let was = baseline == expected;
        for k in ["all", region] {
            for (slot, inc) in counts.entry(k.into()).or_default().iter_mut().zip([
                1, usize::from(matched), usize::from(was), usize::from(matched && !was), usize::from(!matched && was)]) { *slot += inc; }
        }
        writeln!(file,"{}",json!({"trace_index":index,"id":row["id"],"args":row["args"],"baseline":baseline,"expected":expected,"prediction":prediction,"region":region,"matched":matched})).unwrap();
    }
    let result = json!({"complete":true,"input_rows":rows.len(),"counts_all_exact_baseline_gains_losses":counts,
        "source":"ITK v5.4.0 311b7060ef39e371f3cd209ec135284ff5fde735 itkGaussianOperator.hxx:ModifiedBesselI and I0, Apache-2.0; Copyright NumFOCUS",
        "limits":"One complete fixed published higher-order numerical program, strict published I0 with worksheet EXP, unchanged outer order/flush/admission and existing order0/1 product routes. All original records checked; no fitted constant, grouping, cutoff or row-specific selector. Discovery, not held-out."});
    json_write(dest.join("summary.json"), &result);
    println!("{result}");
}

fn published_i(x: f64, order: usize) -> f64 {
    let a = x.abs();
    if order == 0 {
        return published_i0(x);
    }
    let result = if order == 1 {
        if a < 3.75 {
            let y = x / 3.75;
            a * polynomial(y * y, I1_LOW)
        } else {
            polynomial(3.75 / a, I1_HIGH) * (excel_exp(a) / a.sqrt())
        }
    } else if x == 0.0 {
        0.0
    } else {
        let twice_over_x = 2.0 / a;
        let mut next = 0.0;
        let mut current = 1.0;
        let mut answer = 0.0;
        let start = 2 * (order + (40.0 * order as f64).sqrt() as usize);
        for j in (1..=start).rev() {
            let previous = next + (j as f64 * twice_over_x) * current;
            next = current;
            current = previous;
            if current.abs() > 1e12 {
                answer *= 1e-12;
                current *= 1e-12;
                next *= 1e-12;
            }
            if j == order {
                answer = next;
            }
        }
        answer * (published_i0(x) / current)
    };
    if x < 0.0 && order % 2 == 1 {
        -result
    } else {
        result
    }
}

#[derive(Clone, Debug, Serialize)]
struct IProgram {
    id: usize,
    model: EvalModel,
    store_state: bool,
    kornia: bool,
    itk: bool,
    seed_model: EvalModel,
    seed_store: bool,
    seed_return_store: bool,
}

fn i0_model<F: Fpu>(x: f64, p: &IProgram) -> Val {
    let program = Program { id: 0, polynomial: p.seed_model, reconstruction: p.seed_model,
        assignment_stores: p.seed_store, helpers: 0, association: 0, reduction: 0, group_y0_log: false };
    let m = Machine::<F>::new(&program);
    let mode = p.seed_model;
    let a = Val::f64(x.abs());
    let result = if x.abs() < 3.75 {
        let z = m.assign(m.div(Val::f64(x), Val::f64(3.75), mode));
        let square = m.assign(m.mul(z, z, mode));
        m.poly(square, I0_LOW)
    } else {
        let z = m.assign(m.div(Val::f64(3.75), a, mode));
        let scale = m.div(Val::f64(excel_exp(x.abs())), m.sqrt(a), mode);
        m.mul(scale, m.poly(z, I0_HIGH), mode)
    };
    if p.seed_return_store { Val::f64(result.to_f64::<F>()) } else { result }
}

fn positive_above(value: Val, bound: f64) -> bool {
    // Recurrence operands are positive. Compare their exact significands and
    // exponents, including values too large to publish in binary64.
    let key = |v: Val| {
        let bytes = v.to_ext().0;
        (
            u16::from_le_bytes([bytes[8], bytes[9]]) & 0x7fff,
            u64::from_le_bytes(bytes[..8].try_into().unwrap()),
        )
    };
    let a = key(value);
    if a.0 == 0x7fff && a.1 != 1 << 63 {
        return false;
    }
    a > key(Val::f64(bound))
}

fn published_i_model<F: Fpu>(x: f64, order: usize, p: &IProgram) -> f64 {
    let seed = if p.itk { i0_model::<F>(x, p) } else { Val::f64(published_i0(x)) };
    if p.itk && order == 0 { return seed.to_f64::<F>(); }
    if order <= 1 || x == 0.0 {
        // Kornia v0.6.4 kernels.py, Apache-2.0. Its I0 and low I1
        // expressions agree with the PHP source; high I1 groups the final
        // polynomial*EXP before division. The outer parity contract remains.
        if p.kornia && order == 1 && x.abs() >= 3.75 {
            let a = x.abs();
            let result = (polynomial(3.75 / a, I1_HIGH) * excel_exp(a)) / a.sqrt();
            return if x < 0.0 { -result } else { result };
        }
        return published_i(x, order);
    }
    let strict = Program {
        id: 0,
        polynomial: EvalModel::Strict,
        reconstruction: p.model,
        assignment_stores: p.store_state,
        helpers: 0,
        association: 0,
        reduction: 0,
        group_y0_log: false,
    };
    let m = Machine::<F>::new(&strict);
    let mode = p.model;
    let twice_over_x = m.assign(m.div(Val::f64(2.0), Val::f64(x.abs()), mode));
    let mut next = Val::f64(0.0);
    let mut current = Val::f64(1.0);
    let mut answer = Val::f64(0.0);
    let start = 2 * (order + (40.0 * order as f64).sqrt() as usize);
    for j in (1..=start).rev() {
        let previous = m.add(
            next,
            m.mul(m.mul(Val::f64(j as f64), twice_over_x, mode), current, mode),
            mode,
        );
        next = current;
        current = m.assign(previous);
        if positive_above(current, if p.kornia || p.itk { 1e10 } else { 1e12 }) {
            let factor = Val::f64(if p.kornia || p.itk { 1e-10 } else { 1e-12 });
            answer = m.assign(m.mul(answer, factor, mode));
            current = m.assign(m.mul(current, factor, mode));
            next = m.assign(m.mul(next, factor, mode));
        }
        if j == order {
            answer = next;
        }
    }
    let value = if p.kornia {
        m.div(
            m.mul(answer, seed, mode),
            current,
            mode,
        )
    } else {
        m.mul(
            answer,
            m.div(seed, current, mode),
            mode,
        )
    };
    let result = value.to_f64::<F>();
    if x < 0.0 && order % 2 == 1 {
        -result
    } else {
        result
    }
}

pub fn replay_i_models<F: Fpu>(input: &str, out: &str) {
    run_i_models::<F>(input, out, 0);
}

pub fn replay_i_kornia<F: Fpu>(input: &str, out: &str) {
    run_i_models::<F>(input, out, 1);
}

pub fn replay_i_itk_models<F: Fpu>(input: &str, out: &str) {
    run_i_models::<F>(input, out, 2);
}

pub fn replay_i_itk_complete<F: Fpu>(input: &str, out: &str) {
    run_i_models::<F>(input, out, 3);
}

/// Independently rounded sites in the two fixed published higher-I bodies.
/// The ratio-first normalization is separately labelled as a mathematical identity.
#[derive(Clone, Debug, Serialize)]
struct ISites {
    id: usize,
    source: usize,
    ratio_first: bool,
    stores: u8,
    models: [EvalModel; 7],
}

fn i_sites_value<F: Fpu>(x: f64, order: usize, seed: f64, p: &ISites) -> f64 {
    use crate::{dsl::Op, eval::binop};
    let op = |kind, site: usize, a, b| binop::<F>(kind, p.models[site], a, b);
    let mul = |site, a, b| op(Op::Mul(0, 1), site, a, b);
    let div = |site, a, b| op(Op::Div(0, 1), site, a, b);
    let store = |v: Val, mask| if p.stores & mask != 0 { Val::f64(v.to_f64::<F>()) } else { v };
    let tox = div(0, Val::f64(2.0), Val::f64(x.abs()));
    let mut next = Val::f64(0.0);
    let mut current = Val::f64(1.0);
    let mut answer = Val::f64(0.0);
    let start = 2 * (order + (40.0 * order as f64).sqrt() as usize);
    let (limit, factor) = if p.source == 0 { (1e12, 1e-12) } else { (1e10, 1e-10) };
    for j in (1..=start).rev() {
        let first = mul(1, Val::f64(j as f64), tox);
        let second = mul(2, first, current);
        let previous = op(Op::Add(0, 1), 3, next, second);
        next = store(current, 2);
        current = store(previous, 1);
        if positive_above(current, limit) {
            current = store(mul(4, current, Val::f64(factor)), 1);
            next = store(mul(4, next, Val::f64(factor)), 2);
            answer = store(mul(4, answer, Val::f64(factor)), 4);
        }
        if j == order { answer = store(next, 4); }
    }
    let seed = Val::f64(seed);
    let value = if p.ratio_first {
        mul(6, div(5, answer, current), seed)
    } else if p.source == 0 {
        mul(6, answer, div(5, seed, current))
    } else {
        div(6, mul(5, answer, seed), current)
    };
    let result = value.to_f64::<F>();
    if x < 0.0 && order % 2 == 1 { -result } else { result }
}

pub fn replay_i_sites<F: Fpu>(input: &str, out: &str) {
    use crate::dsl::WIDE_MODELS;
    let destination = Path::new(out);
    fs::create_dir(destination).unwrap();
    let original = lines::<Value>(input);
    let mut rows = Vec::new();
    for (index, record) in original.iter().enumerate() {
        let row = &record["raw"];
        assert_eq!(row["function"], "BESSELI");
        let (x, n) = (decode(&row["args"][0]).unwrap(), decode(&row["args"][1]).unwrap());
        if x.is_finite() && !x.is_subnormal() && x != 0.0 && x.abs() <= 709.7
            && n.is_finite() && n >= 2.0 && n.trunc() <= 10000.0
            && row["expected"].as_str().unwrap().starts_with("0x") {
            let seed = published_i0(x);
            assert!(seed.is_finite());
            rows.push((index, x, n.trunc() as usize, seed, decode(&row["expected"]).unwrap().to_bits()));
        }
    }
    assert!(!rows.is_empty());
    // Deterministic discovery order: the first previously lost observations,
    // then every other original observation. No answer-dependent program selection.
    rows.sort_by_key(|(index, _, _, _, _)| {
        let r = &original[*index]["raw"];
        (r["actual"] != r["expected"], *index)
    });
    let scope: Vec<_> = rows.iter().map(|r| r.0).collect();
    json_write(destination.join("scope.json"), &scope);
    let per_models = 5usize.pow(7);
    let programs = 2 * 2 * 8 * per_models;
    let mut kills = vec![u32::MAX; programs];
    let mut survivors = Vec::new();
    let mut killed_by = vec![0usize; rows.len()];
    let mut evaluations = 0u64;
    for source in 0..2 {
        for ratio_first in [false, true] {
            for stores in 0..8 {
                for code in 0..per_models {
                    let id = (((source * 2 + usize::from(ratio_first)) * 8 + stores) * per_models) + code;
                    let mut digits = code;
                    let mut models = [EvalModel::Strict; 7];
                    for m in &mut models { *m = WIDE_MODELS[digits % 5]; digits /= 5; }
                    let p = ISites { id, source, ratio_first, stores: stores as u8, models };
                    for (order, &(_, x, n, seed, expected)) in rows.iter().enumerate() {
                        let got = i_sites_value::<F>(x, n, seed, &p);
                        let got = if got.is_subnormal() || got == 0.0 { 0.0 } else { got };
                        evaluations += 1;
                        if got.to_bits() != expected {
                            kills[id] = order as u32;
                            killed_by[order] += 1;
                            break;
                        }
                    }
                    if kills[id] == u32::MAX { survivors.push(p); }
                }
            }
        }
    }
    let mut file = BufWriter::new(fs::File::create_new(destination.join("first-kill-order.u32le")).unwrap());
    for k in &kills { file.write_all(&k.to_le_bytes()).unwrap(); }
    file.flush().unwrap();
    json_write(destination.join("survivors.json"), &survivors);
    let summary = json!({"complete":true,"original_rows":original.len(),"scoped_rows":rows.len(),"programs":programs,
        "survivors":survivors.len(),"evaluations":evaluations,"killed_by_order":killed_by,
        "grammar":{"models":WIDE_MODELS,"sites":["2/x","j*tox","coefficient*current","next+product","rescale","normalization_first","normalization_second"],"state_store_bits":["current","next","selected"],"source_programs":["PhpSpreadsheet1e12 and source normalization","Kornia0.6.4 1e10 and source normalization"],"ratio_first":[false,true],"program_id":"(((source*2+ratio_first)*8+stores)*5^7)+little_endian_base5_models"},
        "limits":"Declared higher-order bodies, fixed published starts, rescale factors and strict published I0 seed with worksheet EXP. Only admitted nonzero normal x and numeric results with truncated orders2..10000. Independent precision/store sites and a separately labelled ratio identity; no fitted coefficient, cutoff or row-specific policy. A full survivor would still require application-route and full-bank checks; first contradictions exclude only these programs."});
    json_write(destination.join("summary.json"), &summary);
    println!("{summary}");
}

#[cfg(test)]
mod i_site_controls {
    use super::*;

    #[test]
    fn strict_sites_reproduce_both_published_bodies() {
        for source in 0..2 {
            let p = ISites { id: 0, source, ratio_first: false, stores: 0,
                models: [EvalModel::Strict; 7] };
            let old = IProgram { id: 0, model: EvalModel::Strict,
                store_state: false, kornia: source == 1, itk: false,
                seed_model: EvalModel::Strict, seed_store: false, seed_return_store: false };
            for x in [0.001, 0.3, 2.0, 3.75, 5.2, 8.0, 12.3, 25.6, 100.0, -9.1] {
                for n in [2, 3, 10, 50, 100] {
                    let a = i_sites_value::<Hard>(x, n, published_i0(x), &p);
                    let b = published_i_model::<Hard>(x, n, &old);
                    assert_eq!(a.to_bits(), b.to_bits(), "source{source}, x{x}, n{n}");
                }
            }
        }
    }
}

fn run_i_models<F: Fpu>(input: &str, out: &str, source: usize) {
    use calc_func::value::CalcValue;
    let destination = Path::new(out);
    fs::create_dir(destination).unwrap();
    let rows = lines::<Value>(input);
    let mut programs = Vec::new();
    for seed_model in if source == 3 { WIDE_MODELS.to_vec() } else { vec![EvalModel::Strict] } {
      for seed_store in if source == 3 { vec![false, true] } else { vec![false] } {
       for seed_return_store in if source == 3 { vec![false, true] } else { vec![false] } {
    for model in WIDE_MODELS {
        for store_state in [false, true] {
            programs.push(IProgram {
                id: programs.len(),
                model,
                store_state,
                kornia: source == 1,
                itk: source >= 2,
                seed_model,
                seed_store,
                seed_return_store,
            });
        }
    }
       }
      }
    }
    json_write(destination.join("programs.json"), &programs);
    let mut counts = vec![[0_usize; 5]; programs.len()];
    let mut first = vec![None; programs.len()];
    let mut file =
        BufWriter::new(fs::File::create_new(destination.join("predictions.jsonl")).unwrap());
    for (index, record) in rows.iter().enumerate() {
        let row = &record["raw"];
        assert_eq!(row["function"], "BESSELI");
        let a = row["args"].as_array().unwrap();
        let (x, n) = (decode(&a[0]).unwrap(), decode(&a[1]).unwrap());
        let (baseline, _) =
            crate::invoke_control("BESSELI", vec![CalcValue::number(x), CalcValue::number(n)]);
        let expected = row["expected"].as_str().unwrap();
        let x = if x.is_subnormal() { 0.0 } else { x };
        let admitted = x.is_finite()
            && n.is_finite()
            && n >= 0.0
            && n.trunc() <= 10_000_000.0
            && x.abs() <= 709.7;
        for p in &programs {
            let prediction = if !admitted {
                "error:Num".into()
            } else {
                let v = published_i_model::<F>(x, n.trunc() as usize, p);
                if p.id == 0 && source != 1 {
                    assert_eq!(
                        v.to_bits(),
                        if source >= 2 && n.trunc() >= 2.0 { published_i_itk(x, n.trunc() as usize) } else { published_i(x, n.trunc() as usize) }.to_bits(),
                        "strict source control"
                    );
                }
                if !v.is_finite() {
                    "error:Num".into()
                } else {
                    calc_verdict::encode_result(&CalcValue::number(if v.is_subnormal() {
                        0.0
                    } else {
                        v
                    }))
                    .unwrap()
                }
            };
            let matched = prediction == expected;
            let was = baseline == expected;
            for (slot, increment) in counts[p.id].iter_mut().zip([
                1,
                usize::from(matched),
                usize::from(was),
                usize::from(matched && !was),
                usize::from(!matched && was),
            ]) {
                *slot += increment;
            }
            if !matched && first[p.id].is_none() {
                first[p.id] = Some(
                    json!({"trace_index":index,"id":row["id"],"expected":expected,"prediction":prediction}),
                );
            }
            writeln!(file, "{}", json!({"trace_index":index,"id":row["id"],"program":p.id,"prediction":prediction,"expected":expected,"baseline":baseline,"matched":matched})).unwrap();
        }
    }
    let result = json!({"complete":true,"input_rows":rows.len(),"programs":programs,"counts_all_exact_baseline_gains_losses":counts,"first_contradictions":first,
        "source":if source >= 2 {"ITK v5.4.0 311b7060ef39e371f3cd209ec135284ff5fde735 itkGaussianOperator.hxx, Apache-2.0; existing worksheet EXP and unchanged outer admissions; IEEE754/Intel x87 models"} else if source == 1 {"Kornia v0.6.4 kernels.py5dc3c71a5385aa1cdebb39b7aa4f47a0663c9661 (Apache-2.0); binary64 caller dtype, existing worksheet EXP, outer Excel parity/admission; IEEE754/Intel x87 models"} else {"Fixed PhpSpreadsheet BesselI.php60cc826911f1a1a125d756f22f0476a88504ab1f (MIT); IEEE754 and Intel x87 arithmetic/store models"},
        "limits":"Coherent arithmetic/state-store realizations of the unchanged published higher-order body; source3 also gives its complete I0 helper the declared coherent arithmetic, local and return stores, including its public order0 observations. Worksheet EXP and outer policy remain fixed. No coefficient, start, rescale or row-selector search. All answers are discovery."});
    json_write(destination.join("summary.json"), &result);
    println!("{result}");
}

pub fn replay_i_source(input: &str, out: &str) {
    use calc_func::value::CalcValue;
    let destination = Path::new(out);
    fs::create_dir(destination).unwrap();
    let rows = lines::<Value>(input);
    let mut groups = BTreeMap::<String, [usize; 5]>::new();
    let mut file =
        BufWriter::new(fs::File::create_new(destination.join("predictions.jsonl")).unwrap());
    for (index, record) in rows.iter().enumerate() {
        let row = &record["raw"];
        assert_eq!(row["function"], "BESSELI");
        let a = row["args"].as_array().unwrap();
        let (x, n) = (decode(&a[0]).unwrap(), decode(&a[1]).unwrap());
        let (baseline, _) =
            crate::invoke_control("BESSELI", vec![CalcValue::number(x), CalcValue::number(n)]);
        let expected = row["expected"].as_str().unwrap();
        let x = if x.is_subnormal() { 0.0 } else { x };
        // These are the existing product's outer admissions, not a new
        // published numerical branch or an observed-input selector.
        let admitted = x.is_finite()
            && n.is_finite()
            && n >= 0.0
            && n.trunc() <= 10_000_000.0
            && x.abs() <= 709.7;
        let region = if !admitted {
            "outer_policy".into()
        } else {
            format!(
                "I{}:{}",
                if n.trunc() <= 1.0 {
                    n.trunc().to_string()
                } else {
                    "n".into()
                },
                if x.abs() < 3.75 { "low" } else { "high" }
            )
        };
        let prediction = if !admitted {
            "error:Num".into()
        } else {
            let v = published_i(x, n.trunc() as usize);
            if !v.is_finite() {
                "error:Num".into()
            } else {
                calc_verdict::encode_result(&CalcValue::number(if v.is_subnormal() {
                    0.0
                } else {
                    v
                }))
                .unwrap()
            }
        };
        let matched = prediction == expected;
        let was = baseline == expected;
        for key in [
            "all".to_owned(),
            region.clone(),
            format!("{}|{}", record["cohort"].as_str().unwrap(), region),
        ] {
            for (slot, increment) in groups.entry(key).or_default().iter_mut().zip([
                1,
                usize::from(matched),
                usize::from(was),
                usize::from(matched && !was),
                usize::from(!matched && was),
            ]) {
                *slot += increment;
            }
        }
        writeln!(file, "{}", json!({"trace_index":index,"id":row["id"],"prediction":prediction,"expected":expected,"baseline":baseline,"matched":matched,"region":region})).unwrap();
    }
    let result = json!({"complete":true,"input_rows":rows.len(),"counts_all_exact_baseline_gains_losses":groups,
        "source":"PhpSpreadsheet BesselI.php60cc826911f1a1a125d756f22f0476a88504ab1f (MIT), existing worksheet EXP and outer policy",
        "limits":"Whole published numerical program, fixed start/rescale/constants; no fitted program parameter or region selection. Discovery, no held-out claim. Losses reject this whole program."});
    json_write(destination.join("summary.json"), &result);
    println!("{result}");
}
