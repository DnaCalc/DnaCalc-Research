//! Bounded normal-CDF sketches from Cody's published rational tables.
//! CALERF: netlib SPECFUN/erf; direct tables: SciPy1.14.1 cdflib.c:cumnor.
//! See notices/SciPy-cdflib.txt. The later inverse studies use the current
//! cleared MIT Acklam body and the separately pinned CPython AS241 source.
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

// Complete Wichura AS241 coefficient sets as published in CPython3.13.0,
// Modules/_statisticsmodule.c,60403a5409ff2c3f3b07dd2ca91a7a3e096839c7.
// All literals load binary64; no fitted coefficient or region boundary.
const AS241:[[f64;8];6]=[
 [2.5090809287301226727e3,3.3430575583588128105e4,6.7265770927008700853e4,4.5921953931549871457e4,1.3731693765509461125e4,1.9715909503065514427e3,1.3314166789178437745e2,3.3871328727963666080],
 [5.2264952788528545610e3,2.8729085735721942674e4,3.9307895800092710610e4,2.1213794301586595867e4,5.3941960214247511077e3,6.8718700749205790830e2,4.2313330701600911252e1,1.0],
 [7.74545014278341407640e-4,2.27238449892691845833e-2,2.41780725177450611770e-1,1.27045825245236838258,3.64784832476320460504,5.76949722146069140550,4.63033784615654529590,1.42343711074968357734],
 [1.05075007164441684324e-9,5.47593808499534494600e-4,1.51986665636164571966e-2,1.48103976427480074590e-1,6.89767334985100004550e-1,1.67638483018380384940,2.05319162663775882187,1.0],
 [2.01033439929228813265e-7,2.71155556874348757815e-5,1.24266094738807843860e-3,2.65321895265761230930e-2,2.96560571828504891230e-1,1.78482653991729133580,5.46378491116411436990,6.65790464350110377720],
 [2.04426310338993978564e-15,1.42151175831644588870e-7,1.84631831751005468180e-5,7.86869131145613259100e-4,1.48753612908506148525e-2,1.36929880922735805310e-1,5.99832206555887937690e-1,1.0]];

#[derive(Clone,Copy,Debug,Serialize)]
struct InverseProgram {model:EvalModel,local_store:bool,log:usize}
fn inverse_store<F:Fpu>(v:Val,p:InverseProgram)->Val {
    if p.local_store {Val::f64(v.to_f64::<F>())}else{v}
}
fn inverse_abs_le(v:Val,b:f64)->bool {
    let key=|e:Extended|(u16::from_le_bytes(e.0[8..10].try_into().unwrap())&0x7fff,
        u64::from_le_bytes(e.0[..8].try_into().unwrap()));
    key(v.to_ext())<=key(research_load(b))
}
fn inverse_as241<F:Fpu>(probability:f64,p:InverseProgram)->f64 {
    let op=|kind,a,b|binop::<F>(kind,p.model,a,b);
    let add=|a,b|op(Op::Add(0,1),a,b);let mul=|a,b|op(Op::Mul(0,1),a,b);
    let q=inverse_store::<F>(add(Val::f64(probability),Val::f64(-0.5)),p);
    let (r,tables,central)=if inverse_abs_le(q,0.425) {
        (inverse_store::<F>(add(Val::f64(0.180625),neg(mul(q,q))),p),0,true)
    }else{
        let tail=if q.to_f64::<F>()<=0.0{Val::f64(probability)}else{add(Val::f64(1.0),Val::f64(-probability))};
        let tail=inverse_store::<F>(tail,p);
        let log=match p.log {
            0=>Val::f64(crt_log(tail.to_f64::<F>())),
            1=>Val::f64(excel_log(tail.to_f64::<F>())),
            2=>Val::Ext(research_log::<F>(&tail.to_ext())),
            _=>unreachable!(),
        };
        let cw=if matches!(p.model,EvalModel::StoredPc53|EvalModel::ResidentPc53){(CW&!0x300)|0x200}else{CW};
        let root=if p.model==EvalModel::Strict {Val::f64((-log.to_f64::<F>()).sqrt())}
            else {let v=Val::Ext(F::sqrt(&neg(log).to_ext(),cw));if matches!(p.model,EvalModel::StoredPc64|EvalModel::StoredPc53){Val::f64(v.to_f64::<F>())}else{v}};
        let root=inverse_store::<F>(root,p);
        let (shift,tables)=if inverse_abs_le(root,5.0){(1.6,2)}else{(5.0,4)};
        (inverse_store::<F>(add(root,Val::f64(-shift)),p),tables,false)
    };
    let polynomial=|table:usize| {
        AS241[table][1..].iter().fold(Val::f64(AS241[table][0]),|value,&c|add(mul(value,r),Val::f64(c)))
    };
    let mut numerator=polynomial(tables);
    if central {numerator=mul(numerator,q);}
    numerator=inverse_store::<F>(numerator,p);
    let denominator=inverse_store::<F>(polynomial(tables+1),p);
    let value=inverse_store::<F>(op(Op::Div(0,1),numerator,denominator),p);
    let value=if !central&&q.to_f64::<F>()<0.0{neg(value)}else{value};
    value.to_f64::<F>()
}

fn inverse_acklam<F:Fpu>(probability:f64,p:InverseProgram)->f64 {
    // Published Acklam tables from the current MIT cleared implementation.
    const A:&[f64]=&[-3.969683028665376e1,2.209460984245205e2,-2.759285104469687e2,1.38357751867269e2,-3.066479806614716e1,2.506628277459239];
    const B:&[f64]=&[-5.447609879822406e1,1.615858368580409e2,-1.556989798598866e2,6.680131188771972e1,-1.328068155288572e1,1.0];
    const C:&[f64]=&[-7.784894002430293e-3,-3.223964580411365e-1,-2.400758277161838,-2.549732539343734,4.374664141464968,2.938163982698783];
    const D:&[f64]=&[7.784695709041462e-3,3.224671290700398e-1,2.445134137142996,3.754408661907416,1.0];
    let op=|kind,a,b|binop::<F>(kind,p.model,a,b);
    let add=|a,b|op(Op::Add(0,1),a,b);let mul=|a,b|op(Op::Mul(0,1),a,b);let div=|a,b|op(Op::Div(0,1),a,b);
    let cw=if matches!(p.model,EvalModel::StoredPc53|EvalModel::ResidentPc53){(CW&!0x300)|0x200}else{CW};
    let root=|v:Val|if p.model==EvalModel::Strict{Val::f64(v.to_f64::<F>().sqrt())}else{
        let v=Val::Ext(F::sqrt(&v.to_ext(),cw));if matches!(p.model,EvalModel::StoredPc64|EvalModel::StoredPc53){Val::f64(v.to_f64::<F>())}else{v}};
    let poly=|table:&[f64],r|table[1..].iter().fold(Val::f64(table[0]),|v,&c|add(mul(v,r),Val::f64(c)));
    let log=|v:Val|match p.log{0=>Val::f64(crt_log(v.to_f64::<F>())),1=>Val::f64(excel_log(v.to_f64::<F>())),2=>Val::Ext(research_log::<F>(&v.to_ext())),_=>unreachable!()};
    let x=if probability<0.02425 || probability>1.0-0.02425 {
        let tail=if probability<0.02425{Val::f64(probability)}else{add(Val::f64(1.0),Val::f64(-probability))};
        let q=inverse_store::<F>(root(mul(Val::f64(-2.0),log(tail))),p);
        let value=div(poly(C,q),poly(D,q));if probability<0.02425{value}else{neg(value)}
    }else{
        let q=inverse_store::<F>(add(Val::f64(probability),Val::f64(-0.5)),p);
        let r=inverse_store::<F>(mul(q,q),p);div(mul(poly(A,r),q),poly(B,r))
    };
    let x=inverse_store::<F>(x,p);let xf=x.to_f64::<F>();
    // The erfc and PHI operands cross the same binary64 function boundaries
    // as the original cleared refinement. Their producers are conditional.
    let argument=div(Val::f64(xf.abs()),Val::f64(std::f64::consts::SQRT_2)).to_f64::<F>();
    let tail=inverse_store::<F>(mul(Val::f64(0.5),Val::f64(libm::erfc(argument))),p);
    let error=if xf<0.0{add(tail,Val::f64(-probability))}else{add(add(Val::f64(1.0),Val::f64(-probability)),neg(tail))};
    let error=inverse_store::<F>(error,p);
    let sq=excel_x87_mul(xf,xf);let density=excel_x87_mul(excel_exp(-sq/2.0),f64::from_bits(0x3fd9884533d43651));
    if density.abs()<f64::MIN_POSITIVE {xf}else{add(x,neg(div(error,Val::f64(density)))).to_f64::<F>()}
}

pub fn inverse_programs<F:Fpu>(input:&str,output:&str,acklam:bool) {
    let destination=Path::new(output);fs::create_dir(destination).unwrap();
    let rows=crate::lines::<Value>(input);let mut programs=Vec::new();
    for model in WIDE_MODELS {for local_store in [true,false] {for log in 0..3 {programs.push(InverseProgram{model,local_store,log});}}}
    let mut counts=vec![[0usize;5];programs.len()];let mut first=vec![None;programs.len()];let mut excluded=BTreeMap::<String,usize>::new();
    let mut out=BufWriter::new(fs::File::create_new(destination.join("predictions.jsonl")).unwrap());
    for (index,row) in rows.iter().enumerate() {
        let name=row["function"].as_str().unwrap();
        if !matches!(name,"NORM.S.INV"|"NORMSINV"|"NORM.INV"|"NORMINV"|"LOGNORM.INV"|"LOGINV"|"CONFIDENCE"|"CONFIDENCE.NORM") || !row["result_index"].is_null() {
            *excluded.entry("surface_or_selection".into()).or_default()+=1;continue;
        }
        let a=row["args"].as_array().unwrap();
        let numeric:Option<Vec<_>>=a.iter().map(|v|v.as_str().and_then(|s|s.strip_prefix("0x")).and_then(|s|u64::from_str_radix(s,16).ok()).map(f64::from_bits)).collect();
        let Some(a)=numeric else{*excluded.entry("not_scalar_numeric".into()).or_default()+=1;continue;};
        if !a.iter().all(|v|v.is_finite()) || !(0.0<a[0]&&a[0]<1.0) || (a.len()==3&&a[2]<=0.0) {
            *excluded.entry("outside_finite_domain".into()).or_default()+=1;continue;
        }
        let arguments=row["args"].as_array().unwrap().iter().enumerate().map(|(i,v)|calc_verdict::decode_argument(v,&format!("args[{i}]")).unwrap()).collect();
        let (baseline,_)=crate::invoke_control(name,arguments);assert_eq!(row["actual"],baseline,"changed frozen inverse baseline");
        let expected=row["expected"].as_str().unwrap();let was=baseline==expected;
        for (id,p) in programs.iter().enumerate() {
            let confidence=matches!(name,"CONFIDENCE"|"CONFIDENCE.NORM");
            let probability=if confidence{1.0-a[0]/2.0}else{a[0]};
            let z=if acklam{inverse_acklam::<F>(probability,*p)}else{inverse_as241::<F>(probability,*p)};
            // The existing public wrappers keep their binary64 boundaries.
            let value=if confidence {z*a[1]/a[2].trunc().sqrt()}
                else if a.len()==1 {z}
                else if matches!(name,"LOGNORM.INV"|"LOGINV") {excel_exp(a[1]+a[2]*z)}
                else {a[1]+a[2]*z};
            let got=if value.is_finite(){format!("0x{:016x}",value.to_bits())}else{"error:Num".into()};
            if acklam&&id==0 {assert_eq!(got,baseline,"strict original Acklam control {}",row["id"]);}
            let now=got==expected;
            for(slot,v)in counts[id].iter_mut().zip([true,now,was,now&&!was,was&&!now]){*slot+=usize::from(v);}
            if !now&&first[id].is_none(){first[id]=Some(json!({"trace_index":index,"id":row["id"],"got":got,"expected":expected}));}
            writeln!(out,"{}",json!({"trace_index":index,"id":row["id"],"function":name,"program":id,"prediction":got,"expected":expected,"baseline":baseline})).unwrap();
        }
    }
    assert!(counts[0][0]>0);assert_eq!(counts[0][0]+excluded.values().sum::<usize>(),rows.len());
    let safe:Vec<_>=counts.iter().enumerate().filter(|(_,c)|c[3]>0&&c[4]==0).map(|(i,_)|i).collect();
    let result=json!({"complete":true,"input_rows":rows.len(),"programs":programs,"counts_all_exact_baseline_gains_losses":counts,
        "first_contradictions":first,"excluded":excluded,"zero_loss_improvements":safe,
        "algorithm":if acklam{"Acklam with original conditional refinement"}else{"AS241"},
        "source":if acklam{"Published Acklam tables and cleared MIT product program at1ebdbe85 normal_log_family.rs; original libm erfc/current PHI leaves."}else{"CPython3.13.0 commit60403a5409ff2c3f3b07dd2ca91a7a3e096839c7 _statisticsmodule.c, PSF2; complete published Wichura AS241, fixed binary64 literals and region boundaries."},
        "limits":"Five coherent arithmetic models, local stores on/off, three declared conditional log adapters; existing public affine/EXP/confidence wrappers retained. Discovery corpus only; any survivor needs application tests/full bank and independent arithmetic verification. No fitted coefficient or threshold."});
    json_write(destination.join("summary.json"),&result);println!("{result}");
}

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
