//! Trace the registered/workbook routes. Explicit conditional producer studies
//! use a generated research copy and require an unchanged registered baseline.
use crate::ported::complex_family::{ParsedComplex, format_complex_text, parse_complex_text};
use crate::{Given, json_write, ported::financial_time_value_family as tv};
use calc_func::locale::{CalcEnvironment, DateSystem};
use calc_func::{EvalServices, FunctionCallScratch, Notation, catalogue, invoke};
use calc_host::value::{CalcValue, CellKey, CoreRef, NumericMode};
use calc_host::{CellEdit, CellInput, CellText, FiniteNumber, ProjectedValue, SkinIntent};
use dna_calc_worker::workbook::Workbook;
use serde_json::{Value, json};

/// Coherent stored-PC64 arithmetic in the existing A&S K program. Two
/// templates retain or stage SQRT; EXP/LOG and coefficient tables stay fixed.
#[cfg(bessel_k_audit)]
pub fn bessel_k_precision(input:&str,output:&str) {
    let destination=std::path::Path::new(output);std::fs::create_dir(destination).unwrap();
    let rows=crate::lines::<Value>(input);let mut counts=[[0_usize;5];2];
    let mut changes=Vec::new();let mut excluded=0;
    for row in &rows {
        assert_eq!(row["function"],"BESSELK");
        let args:Vec<_>=row["args"].as_array().unwrap().iter()
            .map(|v|calc_verdict::decode_argument(v,"K precision").unwrap()).collect();
        let Some(numbers)=args.iter().map(|a|a.as_number()).collect::<Option<Vec<_>>>() else {excluded+=1;continue;};
        let baseline=crate::invoke_control("BESSELK",args).0;
        assert_eq!(baseline,row["actual"],"original registered K return");
        let format=|answer:Result<f64,calc_func::value::WorksheetErrorCode>| {
            answer.map(|x|format!("0x{:016x}",x.to_bits())).unwrap_or_else(|_|"error:Num".into())
        };
        assert_eq!(format(crate::ported::bessel::besselk_kernel(numbers[0],numbers[1])),baseline,"unchanged K body");
        let predictions=[crate::ported::bessel_k_pc64::besselk_kernel(numbers[0],numbers[1]),
                         crate::ported::bessel_k_pc64_sqrt::besselk_kernel(numbers[0],numbers[1])];
        for (id,prediction)in predictions.into_iter().map(format).enumerate() {
            let now=prediction==row["expected"].as_str().unwrap();let old=baseline==row["expected"].as_str().unwrap();
            for (j,v) in [true,now,old,now&&!old,old&&!now].into_iter().enumerate(){counts[id][j]+=usize::from(v);}
            if prediction!=baseline {changes.push(json!({"program":id,"original":row,"prediction":prediction}));}
        }
    }
    crate::json_write(destination.join("changed-predictions.json").to_str().unwrap(),&changes);
    let summary=json!({"complete":true,"program_labels":2,"original_rows":rows.len(),"excluded_still_fail":excluded,
        "count_order":["all","exact","baseline_exact","gains","losses"],"counts":counts,
        "sources":"Existing original A&S9.8 K program and original Excel readings; shared IEEE754/x87 arithmetic. No coefficient or row selector changed.",
        "limits":"All D add/subtract/multiply/divide operations store from PC64; square root is either the current correctly-rounded binary64 leaf or shared staged SQRT. Existing fixed EXP/LOG and seed tables remain conditional numerical leaves. Application admission/flush/order/recurrence rules retained. Complete baseline-change records are retained, including losses. This is a finite candidate comparison, not whole-family exclusion."});
    crate::json_write(destination.join("summary.json").to_str().unwrap(),&summary);println!("{summary}");
}
use std::{fs, path::Path};

/// Fixed cleared CRT leaves for the uncancelled shape-one density expression.
/// These are discovery hypotheses; only `baseline` uses the registered application route.
pub fn gamma_unit_crt(input:&str,out:&str) {
    use crate::ported::excel_numeric::{crt_log,crt_exp,excel_x87_sub,excel_x87_div,excel_x87_mul,excel_pow_chain,excel_pow_x87_direct};
    let doc:Value=serde_json::from_slice(&fs::read(input).unwrap()).unwrap();
    assert_eq!(doc["complete"],true,"complete original study required");
    let rows=doc["records"].as_array().unwrap();assert!(rows.len()<=2048,"bounded leaf cohort");
    let mut records=Vec::new();let mut counts=std::collections::BTreeMap::<String,[usize;2]>::new();
    for r in rows {
        let parse=|s:&str|f64::from_bits(u64::from_str_radix(s.trim_start_matches("0x"),16).unwrap());
        let x=parse(r["input"].as_str().unwrap());assert!(x>0.0&&x.is_finite());
        let public=parse(r["public_LN"]["bits"].as_str().unwrap());
        let logarithm=crt_log(x);let target=r["density"]["bits"].as_str().unwrap();
        let mut predictions=std::collections::BTreeMap::new();
        for (log_name,log) in [("crt_LN",logarithm),("public_LN",public)] {
            for (arg_name,arg) in [("strict_argument",log-x),("pc64_argument",excel_x87_sub(log,x))] {
                let prefix=crt_exp(arg);
                for (div_name,value) in [("strict_division",prefix/x),("pc64_division",excel_x87_div(prefix,x))] {
                    let name=format!("{log_name}/{arg_name}/{div_name}");let bits=format!("{:016x}",value.to_bits());
                    let count=counts.entry(name.clone()).or_default();count[0]+=1;count[1]+=usize::from(bits==target);
                    predictions.insert(name,bits);
                }
            }
        }
        let observed_exp=parse(r["public_EXP"]["bits"].as_str().unwrap());
        for (power_name,power) in [("raw_log_power",excel_pow_chain(x,1.0)),("raw_log2_power",excel_pow_x87_direct(x,1.0))] {
            for (model_name,stored) in [("strict",false),("stored_pc64",true)] {
                let mul=|a,b|if stored { excel_x87_mul(a,b) } else { a*b };
                let div=|a,b|if stored { excel_x87_div(a,b) } else { a/b };
                for (association,value) in [("product_first",div(mul(power,observed_exp),x)),
                    ("power_ratio_first",mul(div(power,x),observed_exp)),
                    ("exp_ratio_first",mul(power,div(observed_exp,x)))] {
                    let name=format!("{power_name}/{model_name}/{association}");let bits=format!("{:016x}",value.to_bits());
                    let count=counts.entry(name.clone()).or_default();count[0]+=1;count[1]+=usize::from(bits==target);
                    predictions.insert(name,bits);
                }
            }
        }
        let baseline=crate::invoke_control("GAMMA.DIST",[x,1.0,1.0,0.0].into_iter().map(CalcValue::number).collect()).0;
        records.push(json!({"input":r["input"],"density":r["density"],"public_LN":r["public_LN"],
            "crt_LN":format!("{:016x}",logarithm.to_bits()),"crt_exp_minus_x":format!("{:016x}",crt_exp(-x).to_bits()),
            "baseline":baseline,"hypotheses":predictions}));
    }
    let result=json!({"complete":true,"records":records,"counts_rows_matched":counts,"not_excel_evidence":true,
        "limits":"Eight fixed cleared-CRT compositions of EXP(LN(x)-x)/x and twelve raw-power/product associations using independently measured public EXP. Source-shape order only. No coefficient, normalizer or selector fitted. A survivor identifies only this sampled functional expression. The separate baseline is the registered application calculation; hypothesis matches are not product gains."});
    crate::json_write(out,&result);println!("{}",result["counts_rows_matched"]);
}

/// Replay the unchanged cleared prefix only to locate recorded consumer-step losses.
/// Every reported before/after return is checked against its canonical application evidence.
pub fn gammaln_low_stages(input: &str, out: &str) {
    use crate::ported::excel_numeric::{crt_log, excel_log, excel_x87_sub, gammaln_excel};
    let doc: Value = serde_json::from_slice(&fs::read(input).unwrap()).unwrap();
    let mut records = Vec::new();
    let mut generated = Vec::new();
    if let Some(count) = doc["rng_count"].as_u64() {
        assert!(count <= 4_194_304, "bounded primitive mining only");
        let mut state = doc["rng_seed"].as_u64().unwrap();
        for _ in 0..count {
            state ^= state << 13; state ^= state >> 7; state ^= state << 17;
            let exponent = 1023 - (1 + (state >> 52) % 52);
            let bits = (exponent << 52) | (state & ((1u64 << 52) - 1));
            let x = f64::from_bits(bits);
            if x < 0.7 { generated.push(json!({"x": format!("{bits:016x}")})); }
        }
    }
    let inputs = doc["records"].as_array().unwrap_or(&generated);
    for row in inputs {
        let x = f64::from_bits(u64::from_str_radix(row["x"].as_str().unwrap(), 16).unwrap());
        assert!(x > 0.0 && x < 0.7 && x.is_normal());
        let runtime = crt_log(x);
        let worksheet = excel_log(x);
        if doc["only_log_disagreements"].as_bool() == Some(true) && runtime == worksheet { continue; }
        let upper = gammaln_excel(1.0 + x);
        let (actual, _) = crate::invoke_control("GAMMALN", vec![CalcValue::number(x)]);
        records.push(json!({"original": row, "upper_argument": format!("{:016x}", (1.0+x).to_bits()),
            "upper": format!("{:016x}", upper.to_bits()),
            "runtime_log": format!("{:016x}", runtime.to_bits()),
            "worksheet_log": format!("{:016x}", worksheet.to_bits()),
            "runtime_subtract": format!("{:016x}", (upper-runtime).to_bits()),
            "worksheet_subtract": format!("{:016x}", (upper-worksheet).to_bits()),
            "worksheet_pc64_subtract": format!("{:016x}", excel_x87_sub(upper,worksheet).to_bits()),
            "registered": actual}));
    }
    crate::json_write(out, &json!({"format":"dnacalc-gammaln-low-stage-replay/1", "complete":true,
        "records":records, "eligible_primitive_inputs":inputs.len(), "selection":doc,
        "unit": serde_json::from_str::<Value>(include_str!(concat!(env!("OUT_DIR"),"/research-unit.json"))).unwrap(),
        "limits":"Current cleared stage replay and actual registered returns. Intermediate values are hypotheses, never additional Excel observations."}));
}

pub fn gamma_density_normalizers(input: &str, out: &str) {
    use crate::ported::excel_numeric::{excel_exp_compensated, excel_log, excel_x87_div, excel_x87_mul, gammaln_excel};
    let doc: Value = serde_json::from_slice(&fs::read(input).unwrap()).unwrap();
    let mut records = Vec::new();
    for row in doc["records"].as_array().unwrap() {
        let number = |key: &str| f64::from_bits(u64::from_str_radix(row[key].as_str().unwrap(),16).unwrap());
        let (x,a,b) = (number("x"),number("alpha"),number("beta"));
        let y = excel_x87_div(x,b);
        assert!(x>0.0 && a>0.0 && a!=1.0 && b>0.0 && (a<=1.0 || y<1.1));
        let current = gammaln_excel(a+1.0);
        let public = f64::from_bits(u64::from_str_radix(row["public_normalizer"]["bits"].as_str().unwrap(),16).unwrap());
        let term = excel_x87_mul(a,excel_log(y));
        let factor = excel_x87_div(a,x);
        let replay = |normalizer| {
            let mass=excel_exp_compensated(term,&[y,normalizer]);
            let result=excel_x87_mul(mass,factor);
            if result.is_subnormal() {0.0} else {result}
        };
        let actual = crate::invoke_control("GAMMA.DIST",[x,a,b,0.0].into_iter().map(CalcValue::number).collect()).0;
        let predicted = format!("0x{:016x}",replay(current).to_bits());
        assert_eq!(actual,predicted,"current gamma prefix differs from registered route");
        records.push(json!({"original":row,"registered":actual,
            "current_normalizer":format!("{:016x}",current.to_bits()),
            "public_normalizer_prediction":format!("{:016x}",replay(public).to_bits()),
            "hypothetical_only":true}));
    }
    crate::json_write(out,&json!({"format":"dnacalc-gamma-density-normalizer-replay/1","complete":true,
        "records":records,"unit":serde_json::from_str::<Value>(include_str!(concat!(env!("OUT_DIR"),"/research-unit.json"))).unwrap(),
        "limits":"Actual registered baseline and one diagnostic replacement by an independent public normalizer. The replacement is not an application calculation or a stored Excel-answer implementation."}));
}

pub fn gamma_density_prefixes(input:&str,out:&str) {
    use crate::ported::excel_numeric::{excel_log,excel_exp_compensated,gammaln_excel,excel_x87_div,excel_x87_mul};
    let doc:Value=serde_json::from_slice(&fs::read(input).unwrap()).unwrap();
    let rows=doc["records"].as_array().unwrap();let mut records=Vec::new();
    let encode=|x:f64|if !x.is_finite(){"error:Num".into()}else if x!=0.0&&x.is_subnormal(){"fault:NotImplemented".into()}else{format!("0x{:016x}",if x==0.0 {0.0}else{x}.to_bits())};
    for row in rows {
        let aa=row["args"].as_array().unwrap();let mut ns=Vec::new();let mut unresolved=false;
        for v in aa.iter().take(3) {
            if let Some(s)=v.as_str().and_then(|s|s.strip_prefix("0x")) {
                ns.push(f64::from_bits(u64::from_str_radix(s,16).unwrap()));
            } else if let Some(b)=v.as_bool() {ns.push(if b {1.0}else{0.0});}
            else {unresolved=true;break;}
        }
        if unresolved||ns.len()!=3 {records.push(json!({"original":row,"classification":"nonnumeric_argument_needs_wrapper_audit"}));continue;}
        let (x,a,b)=(ns[0],ns[1],ns[2]);assert!(x>0.0&&a!=1.0&&b>0.0,"changed generic-density route");
        let y=excel_x87_div(x,b);let log=excel_log(y);let normalizer=gammaln_excel(a+1.0);
        let r=excel_exp_compensated(excel_x87_mul(a,log),&[y,normalizer]);let factor=excel_x87_div(a,x);
        let old=encode(excel_x87_div(excel_x87_mul(r,a),x));let now=encode(excel_x87_mul(r,factor));let strict=encode(r*factor);
        let args=aa.iter().map(|v|calc_verdict::decode_argument(v,"original gamma argument").unwrap()).collect();
        let actual=crate::invoke_control(row["function"].as_str().unwrap(),args).0;
        assert_eq!(actual,row["actual"].as_str().unwrap(),"registered loss return changed");
        assert_eq!(now,actual,"prefix replay differs from registered application");
        assert_eq!(old,row["expected"].as_str().unwrap(),"old prefix replay differs from original matched return");
        records.push(json!({"original":row,"prefix":format!("{:016x}",r.to_bits()),"stored_factor":format!("{:016x}",factor.to_bits()),
            "scaled_input":format!("{:016x}",y.to_bits()),"normalizer":format!("{:016x}",normalizer.to_bits()),
            "old_replayed":old,"new_registered":actual,"strict_multiply":strict,"classification":"prefix_replay_verified"}));
    }
    crate::json_write(out,&json!({"complete":true,"records":records,"not_excel_evidence":true,
        "limits":"Cleared-prefix diagnostic replay, checked against each old matched answer and current registered return. Hypothetical intermediate values are not additional engine passes or private Excel values. Final exact inverse classification is recorded separately."}));
}

/// Fixed product-root and natural-log helper boundaries. No coefficients,
/// thresholds or per-observation choices occur in these mathematical paths.
pub fn geomean_stored<F:crate::ported::excel_numeric::Fpu>(input:&str,out:&str) {
    use crate::{dsl::{EvalModel,Op},eval::{Val,binop}};
    use crate::ported::excel_numeric::*;
    let rows=crate::lines::<Value>(input);
    let mut counts=[[0usize;5];20];let mut excluded=std::collections::BTreeMap::<String,usize>::new();
    let mut answers=Vec::new();
    fn flatten(v:&Value,values:&mut Vec<f64>)->bool {
        if let Some(array)=v.as_array(){return array.iter().all(|v|flatten(v,values));}
        let Some(s)=v.as_str()else{return false;};
        let Some(s)=s.strip_prefix("0x")else{return false;};
        let Ok(bits)=u64::from_str_radix(s,16)else{return false;};
        let x=f64::from_bits(bits);if !(x>0.0&&x.is_finite()){return false;}
        values.push(x);true
    }
    for (index,row)in rows.iter().enumerate() {
        if row["function"]!="GEOMEAN"{continue;}
        let mut xs=Vec::new();
        if !row["result_index"].is_null()||!flatten(&row["args"],&mut xs)||xs.is_empty(){
            *excluded.entry("nonpositive_or_unresolved_argument".into()).or_default()+=1;continue;
        }
        let baseline=crate::invoke_control("GEOMEAN",xs.iter().map(|&x|CalcValue::number(x)).collect()).0;
        assert_eq!(baseline,row["actual"],"registered GEOMEAN {}",row["id"]);
        let mut outputs=Vec::new();
        for product_model in [EvalModel::Strict,EvalModel::StoredPc64] {
            let product=xs.iter().fold(Val::f64(1.0),|p,&x|binop::<F>(Op::Mul(0,1),product_model,p,Val::f64(x))).to_f64::<F>();
            if !product.is_normal(){break;}
            let p=Val::f64(product).to_ext();let n=Val::f64(xs.len() as f64).to_ext();
            let store=|x|Val::f64(F::store(&x,0x133f)).to_ext();
            // A separately stored log(product) or quotient is a helper-call
            // boundary, not an adjustment of the mathematical logarithm.
            for (log_store,quotient_store)in [(false,false),(true,false),(false,true),(true,true)] {
                let log=research_log::<F>(&p);let log=if log_store{store(log)}else{log};
                let q=F::div(&log,&n,0x133f);let q=if quotient_store{store(q)}else{q};
                outputs.push(F::store(&research_exp::<F>(&q),0x133f));
            }
            let base2=F::fyl2x(&Val::f64(1.0).to_ext(),&p,0x133f);
            for quotient_store in [false,true] {
                let q=F::div(&base2,&n,0x133f);let q=if quotient_store{store(q)}else{q};
                outputs.push(F::store(&research_pow2::<F>(&q),0x133f));
            }
            let reciprocal=F::store(&F::div(&Val::f64(1.0).to_ext(),&n,0x133f),0x133f);
            outputs.push(F::store(&research_pow::<F>(product,reciprocal),0x133f));
            outputs.push(crt_pow(product,1.0/xs.len() as f64));
            outputs.push(crt_exp(crt_log(product)/xs.len() as f64));
            let logs=xs.iter().map(|&x|F::store(&research_log::<F>(&Val::f64(x).to_ext()),0x133f));
            let sum=logs.fold(Val::f64(0.0),|s,x|binop::<F>(Op::Add(0,1),product_model,s,Val::f64(x)));
            let q=F::store(&F::div(&sum.to_ext(),&n,0x133f),0x133f);
            outputs.push(F::store(&research_exp::<F>(&Val::f64(q).to_ext()),0x133f));
        }
        if outputs.len()!=20{*excluded.entry("product_range".into()).or_default()+=1;continue;}
        let outputs:Vec<_>=outputs.into_iter().map(|x|format!("0x{:016x}",x.to_bits())).collect();
        for (i,got)in outputs.iter().enumerate(){
            let now=*got==row["expected"];let old=baseline==row["expected"];
            for(j,v)in [true,now,old,now&&!old,old&&!now].into_iter().enumerate(){counts[i][j]+=usize::from(v);}
        }
        answers.push(json!({"trace_index":index,"id":row["id"],"expected":row["expected"],"baseline":baseline,"predictions":outputs}));
    }
    fs::create_dir(out).unwrap();crate::json_write(Path::new(out).join("predictions.json"),&answers);
    let summary=json!({"complete":true,"counts_all_exact_baseline_gains_losses":counts,"excluded":excluded,
        "programs":"Two product policies (strict,storedPC64), each with natural product log and quotient store flags00,10,01,11; base2 quotient store0,1; stored reciprocal POW; CRT POW; CRT EXP(LOG(P)/n); stored sum of LN inputs and quotient.",
        "sources":"Defining product-root and exponential/logarithm identities, Intel SDM; existing cleared numerical helpers and original captures. No fitted constants or answer-dependent selector.",
        "limits":"Positive finite numeric operands and normal binary64 products. Helper microcode and CRT bodies remain conditional leaves. Counts are diagnostics; no changing program is adopted solely from these observations."});
    crate::json_write(Path::new(out).join("summary.json"),&summary);println!("{summary}");
}

#[cfg(regression_audit)]
std::thread_local! {
    static REGRESSION_LEAVES: std::cell::RefCell<Vec<Value>> = const { std::cell::RefCell::new(Vec::new()) };
    static REGRESSION_FACTORS: std::cell::RefCell<Option<Value>> = const { std::cell::RefCell::new(None) };
    static REGRESSION_TRIANGULAR: std::cell::RefCell<Option<Value>> = const { std::cell::RefCell::new(None) };
}

#[cfg(regression_householder_audit)]
std::thread_local! {
    static REGRESSION_HOUSEHOLDER_PROGRAM:std::cell::Cell<Option<(usize,crate::dsl::EvalModel,bool)>>=const{std::cell::Cell::new(None)};
}

/// Existing unnormalized reflector, or LAPACK DLARFG's normalized vector/tau
/// representation. Original pivot/mean/rank rules and binary64 squared norms
/// remain conditional. Every matrix/vector write stores binary64.
#[cfg(regression_householder_audit)]
pub fn regression_householder_step(matrix:&mut[(Option<usize>,Vec<f64>)],target:&mut[f64],j:usize,x0:f64,rest:f64,squares:f64)->f64 {
    use crate::dsl::EvalModel;
    use crate::ported::excel_numeric::{Hard,Soft};
    let (shape,model,software)=REGRESSION_HOUSEHOLDER_PROGRAM.with(|p|p.get()).unwrap_or((0,EvalModel::Strict,false));
    if software {householder_step::<Soft>(matrix,target,j,x0,rest,squares,shape,model)}
        else{householder_step::<Hard>(matrix,target,j,x0,rest,squares,shape,model)}
}

#[cfg(regression_householder_audit)]
fn householder_step<F:crate::ported::excel_numeric::Fpu>(matrix:&mut[(Option<usize>,Vec<f64>)],target:&mut[f64],j:usize,x0:f64,rest:f64,squares:f64,shape:usize,model:crate::dsl::EvalModel)->f64 {
    use crate::{dsl::Op,eval::{Val,binop}};
    let add=|a,b|binop::<F>(Op::Add(0,1),model,a,b);
    let sub=|a,b|add(a,crate::boundary::neg(b));
    let mul=|a,b|binop::<F>(Op::Mul(0,1),model,a,b);
    let div=|a,b|binop::<F>(Op::Div(0,1),model,a,b);
    let alpha=-squares.sqrt().copysign(x0);
    if shape==1&&rest==0.0{return x0;} // DLARFG zero-tail identity reflector.
    let mut vector=matrix[j].1[j..].to_vec();
    let multiplier=if shape==0 {
        vector[0]=sub(Val::f64(x0),Val::f64(alpha)).to_f64::<F>();
        add(Val::f64(rest),mul(Val::f64(vector[0]),Val::f64(vector[0])))
    }else{
        assert_eq!(shape,1);
        let denominator=sub(Val::f64(x0),Val::f64(alpha));
        // DSCAL receives a binary64 reciprocal; tau crosses DLARFG's scalar
        // output boundary before applying the reflector to another column.
        let reciprocal=Val::f64(div(Val::f64(1.0),denominator).to_f64::<F>());
        for value in &mut vector[1..]{*value=mul(Val::f64(*value),reciprocal).to_f64::<F>();}
        vector[0]=1.0;
        Val::f64(div(sub(Val::f64(alpha),Val::f64(x0)),Val::f64(alpha)).to_f64::<F>())
    };
    let reflect=|values:&mut[f64]| {
        let dot=vector.iter().zip(&values[j..]).fold(Val::f64(0.0),|sum,(&v,&a)|add(sum,mul(Val::f64(v),Val::f64(a))));
        let factor=if shape==0{div(mul(Val::f64(2.0),dot),multiplier)}else{mul(multiplier,dot)};
        for(a,&v)in values[j..].iter_mut().zip(&vector){*a=sub(Val::f64(*a),mul(factor,Val::f64(v))).to_f64::<F>();}
    };
    for(_,column)in matrix.iter_mut().skip(j+1){reflect(column);}
    reflect(target);alpha
}

#[cfg(regression_householder_audit)]
pub fn regression_householder(input:&str,out:&str,software:bool) {
    use std::io::Write;
    let destination=Path::new(out);fs::create_dir(destination).unwrap();let rows=crate::lines::<Value>(input);
    let mut programs=Vec::new();for shape in 0..2{for model in crate::dsl::WIDE_MODELS{programs.push((shape,model));}}
    let mut counts=vec![[0usize;5];programs.len()];let mut first_losses=vec![None;programs.len()];let mut excluded=std::collections::BTreeMap::<String,usize>::new();
    let evaluate=regression_evaluate;
    let mut file=std::io::BufWriter::new(fs::File::create_new(destination.join("predictions.jsonl")).unwrap());
    for(index,row)in rows.iter().enumerate(){
        let name=row["function"].as_str().unwrap();if !matches!(name,"LINEST"|"LOGEST"|"TREND"|"GROWTH"){*excluded.entry("other_surface".into()).or_default()+=1;continue;}
        let selected=if row["result_index"].is_null(){None}else{let v=row["result_index"].as_array().unwrap();Some((v[0].as_u64().unwrap()as usize,v[1].as_u64().unwrap()as usize))};
        let args:Vec<_>=row["args"].as_array().unwrap().iter().enumerate().map(|(i,v)|calc_verdict::decode_argument(v,&format!("args[{i}]")).unwrap()).collect();
        REGRESSION_HOUSEHOLDER_PROGRAM.with(|p|p.set(None));
        let Some(observed)=evaluate(name,&args,selected)else{*excluded.entry("outside_selectable_result_shape".into()).or_default()+=1;continue;};
        let(baseline,_)=crate::invoke_selected_control(name,args.clone(),selected);
        assert_eq!(baseline,row["actual"].as_str().unwrap(),"frozen registered QR {}",row["id"]);
        assert_eq!(observed,baseline,"unmodified reflector/registered QR {}",row["id"]);
        let expected=row["expected"].as_str().unwrap();let was=baseline==expected;
        for(id,&(shape,model))in programs.iter().enumerate(){
            REGRESSION_HOUSEHOLDER_PROGRAM.with(|p|p.set(Some((shape,model,software))));
            let got=evaluate(name,&args,selected).unwrap_or_else(||"error:Ref".into());
            if id==0{assert_eq!(got,baseline,"original strict QR {}",row["id"]);}
            let now=got==expected;for(slot,b)in counts[id].iter_mut().zip([true,now,was,now&&!was,was&&!now]){*slot+=usize::from(b);}
            if was&&!now&&first_losses[id].is_none(){first_losses[id]=Some(json!({"trace_index":index,"id":row["id"],"got":got,"expected":expected}));}
            writeln!(file,"{}",json!({"trace_index":index,"id":row["id"],"program":id,"prediction":got,"expected":expected,"baseline":baseline})).unwrap();
        }
    }
    REGRESSION_HOUSEHOLDER_PROGRAM.with(|p|p.set(None));assert!(counts[0][0]>0);assert_eq!(counts[0][0]+excluded.values().sum::<usize>(),rows.len());
    let safe:Vec<_>=counts.iter().enumerate().filter(|(_,v)|v[3]>0&&v[4]==0).map(|(i,_)|i).collect();
    let result=json!({"complete":true,"input_rows":rows.len(),"programs":programs,"counts_all_exact_baseline_gains_losses":counts,"first_losses":first_losses,"excluded":excluded,"zero_loss_improvements":safe,
        "source":"Householder I-2vv'/v'v and LAPACK3.12.1 DLARFG ordinary normalized-vector/tau formulas,6ec7f2bc4ecf4c4a93496aa2fa519575bc0e39ca,modifiedBSD.",
        "limits":"10 fixed reflector labels. Original rank,mean,pivot and binary64 squared norms remain conditional; arrays store binary64, normalized reciprocal/tau cross binary64 helper boundaries. This is not a full LAPACK QR implementation. No fitted threshold/coefficient/row selector or new Excel answer."});
    json_write(destination.join("summary.json"),&result);println!("{result}");
}

#[cfg(regression_audit)]
fn regression_evaluate(name:&str,args:&[CalcValue],selected:Option<(usize,usize)>)->Option<String> {
        REGRESSION_LEAVES.with(|v|v.borrow_mut().clear());REGRESSION_FACTORS.with(|v|*v.borrow_mut()=None);REGRESSION_TRIANGULAR.with(|v|*v.borrow_mut()=None);
        let result=match name{
            "LINEST"=>crate::ported::regression_forecast_family::eval_linest_surface(args),
            "LOGEST"=>crate::ported::regression_forecast_family::eval_logest_surface(args),
            "TREND"=>crate::ported::regression_forecast_family::eval_trend_surface(args),
            "GROWTH"=>crate::ported::regression_forecast_family::eval_growth_surface(args),
            _=>unreachable!(),
        };
        let value=result.unwrap_or_else(CalcValue::error);
        // The original statistics adapter refuses any subnormal array cell,
        // converts a nonfinite selected number to NUM, and normalizes zero.
        // Reuse the existing research publication rule after the array check.
        let subnormal=|cell:&CalcValue|cell.as_number().is_some_and(f64::is_subnormal);
        if value.as_array().map_or_else(||subnormal(&value),|a|a.iter_row_major().any(subnormal)){
            return Some("fault:NotImplemented".into());
        }
        calc_verdict::select_result(&value,selected).and_then(|chosen|chosen.as_number().map_or_else(
            ||calc_verdict::encode_result(&chosen),
            |n|Some(crate::eval::outcome::<crate::ported::excel_numeric::Hard>(crate::eval::Val::f64(n)))))

}

#[cfg(regression_outer_audit)]
std::thread_local! {
    static REGRESSION_OUTER_PROGRAM:std::cell::Cell<Option<(usize,bool)>>=const{std::cell::Cell::new(None)};
    static REGRESSION_FOLD_PROGRAM:std::cell::Cell<Option<(usize,bool)>>=const{std::cell::Cell::new(None)};
    static REGRESSION_NORM_PROGRAM:std::cell::Cell<Option<(usize,bool)>>=const{std::cell::Cell::new(None)};
}

/// A fixed PC64 product/sum loop, with binary64 publication at the loop result.
#[cfg(regression_outer_audit)]
pub fn regression_outer_fold(site:usize,initial:f64,terms:impl Iterator<Item=(f64,f64)>)->f64 {
    use crate::{dsl::{EvalModel,Op},eval::{Val,binop},ported::excel_numeric::{Fpu,Hard,Soft}};
    let bit=match site{0=>0,2=>1,4=>2,5=>3,_=>unreachable!()};
    let (mask,software)=REGRESSION_FOLD_PROGRAM.with(|p|p.get()).unwrap_or((0,false));
    if mask&(1<<bit)==0 {
        return terms.fold(initial,|sum,(a,b)|regression_outer_op(site,if site==5{1}else{0},sum,regression_outer_op(site,2,a,b)));
    }
    fn calculate<F:Fpu>(site:usize,initial:f64,terms:impl Iterator<Item=(f64,f64)>)->f64 {
        let model=EvalModel::ResidentPc64;
        terms.fold(Val::f64(initial),|sum,(a,b)|{
            let product=binop::<F>(Op::Mul(0,1),model,Val::f64(a),Val::f64(b));
            let product=if site==5{match product{Val::F64(b)=>Val::F64(b^(1<<63)),Val::Ext(mut v)=>{v.0[9]^=0x80;Val::Ext(v)}}}else{product};
            binop::<F>(Op::Add(0,1),model,sum,product)
        }).to_f64::<F>()
    }
    if software{calculate::<Soft>(site,initial,terms)}else{calculate::<Hard>(site,initial,terms)}
}

/// Operations outside the reflector: all sites in a group use one precision.
#[cfg(regression_outer_audit)]
pub fn regression_outer_op(site:usize,operation:usize,a:f64,b:f64)->f64 {
    use crate::ported::excel_numeric::{Fpu,Hard,Soft,research_load,excel_sqrt_staged};
    let (mask,software)=REGRESSION_OUTER_PROGRAM.with(|p|p.get()).unwrap_or((0,false));
    let norm=REGRESSION_NORM_PROGRAM.with(|p|p.get());
    let (stored,software)=if site==2 {
        norm.map(|(program,soft)|((operation==2&&program&1!=0)||(operation==0&&program&2!=0),soft))
            .unwrap_or((mask&(1<<site)!=0,software))
    }else{(mask&(1<<site)!=0,software)};
    if !stored {return match operation{0=>a+b,1=>a-b,2=>a*b,3=>a/b,4=>a.sqrt(),_=>unreachable!()};}
    fn calculate<F:Fpu>(operation:usize,a:f64,b:f64)->f64 {
        use crate::eval::CW;
        let (a,b)=(research_load(a),research_load(b));
        let value=match operation{0=>F::add(&a,&b,CW),1=>F::sub(&a,&b,CW),2=>F::mul(&a,&b,CW),3=>F::div(&a,&b,CW),4=>F::sqrt(&a,CW),_=>unreachable!()};
        F::store(&value,CW)
    }
    if software{calculate::<Soft>(operation,a,b)}else{calculate::<Hard>(operation,a,b)}
}

#[cfg(regression_outer_audit)]
pub fn regression_outer(input:&str,out:&str,software:bool) {
    regression_outer_study(input,out,software,false,false);
}
#[cfg(regression_outer_audit)]
pub fn regression_folds(input:&str,out:&str,software:bool) {
    regression_outer_study(input,out,software,true,false);
}
#[cfg(regression_outer_audit)]
pub fn regression_norm_stores(input:&str,out:&str,software:bool) {
    regression_outer_study(input,out,software,false,true);
}
#[cfg(regression_outer_audit)]
fn regression_expectation(row:&Value)->Option<&str> {
    if row["unobserved_question"]==true {
        assert!(row["expected"].is_null(),"an unanswered question must not discard a recorded expectation");
        None
    }else{Some(row["expected"].as_str().expect("recorded QR expectation"))}
}

#[cfg(regression_outer_audit)]
fn regression_outer_study(input:&str,out:&str,software:bool,folds:bool,norms:bool) {
    use std::io::Write;
    let destination=Path::new(out);fs::create_dir(destination).unwrap();let rows=crate::lines::<Value>(input);
    let programs=if norms{4}else if folds{16}else{64};
    let mut counts=vec![[0usize;5];programs];let mut first_losses=vec![None;programs];let mut first_gains=vec![None;programs];let mut excluded=std::collections::BTreeMap::<String,usize>::new();
    let mut unobserved=0usize;
    let mut changes=std::io::BufWriter::new(fs::File::create_new(destination.join("changes.jsonl")).unwrap());
    for(index,row)in rows.iter().enumerate(){
        let name=row["function"].as_str().unwrap();if !matches!(name,"LINEST"|"LOGEST"|"TREND"|"GROWTH"){*excluded.entry("other_surface".into()).or_default()+=1;continue;}
        let selected=if row["result_index"].is_null(){None}else{let v=row["result_index"].as_array().unwrap();Some((v[0].as_u64().unwrap()as usize,v[1].as_u64().unwrap()as usize))};
        let args:Vec<_>=row["args"].as_array().unwrap().iter().enumerate().map(|(i,v)|calc_verdict::decode_argument(v,&format!("args[{i}]")).unwrap()).collect();
        REGRESSION_OUTER_PROGRAM.with(|p|p.set(None));
        REGRESSION_FOLD_PROGRAM.with(|p|p.set(None));
        REGRESSION_NORM_PROGRAM.with(|p|p.set(None));
        let Some(observed)=regression_evaluate(name,&args,selected)else{*excluded.entry("outside_selectable_result_shape".into()).or_default()+=1;continue;};
        let(baseline,_)=crate::invoke_selected_control(name,args.clone(),selected);
        assert_eq!(baseline,row["actual"].as_str().unwrap(),"frozen registered outer QR {}",row["id"]);
        assert_eq!(observed,baseline,"unmodified outer/registered QR {}",row["id"]);
        let expected=regression_expectation(row);
        unobserved+=usize::from(expected.is_none());
        let was=expected.is_some_and(|e|baseline==e);
        for mask in 0..programs {
            if norms{REGRESSION_NORM_PROGRAM.with(|p|p.set(Some((mask,software))));}
            else if folds{REGRESSION_FOLD_PROGRAM.with(|p|p.set(Some((mask,software))));}
            else{REGRESSION_OUTER_PROGRAM.with(|p|p.set(Some((mask,software))));}
            let got=regression_evaluate(name,&args,selected).unwrap_or_else(||"error:Ref".into());
            if mask==0{assert_eq!(got,baseline,"original outer QR {}",row["id"]);}
            let now=expected.is_some_and(|e|got==e);for(slot,b)in counts[mask].iter_mut().zip([expected.is_some(),now,was,now&&!was,was&&!now]){*slot+=usize::from(b);}
            if got!=baseline {
                let change=json!({"trace_index":index,"id":row["id"],"program":mask,"prediction":got,"expected":expected,"baseline":baseline});
                if was&&!now&&first_losses[mask].is_none(){first_losses[mask]=Some(change.clone());}
                if now&&!was&&first_gains[mask].is_none(){first_gains[mask]=Some(change.clone());}
                writeln!(changes,"{}",change).unwrap();
            }
        }
    }
    REGRESSION_OUTER_PROGRAM.with(|p|p.set(None));REGRESSION_FOLD_PROGRAM.with(|p|p.set(None));REGRESSION_NORM_PROGRAM.with(|p|p.set(None));assert!(counts[0][0]+unobserved>0);assert_eq!(counts[0][0]+unobserved+excluded.values().sum::<usize>(),rows.len());
    let safe:Vec<_>=counts.iter().enumerate().filter(|(_,v)|v[3]>0&&v[4]==0).map(|(i,_)|i).collect();
    let result=json!({"complete":true,"input_rows":rows.len(),"unobserved_questions_not_scored":unobserved,"programs":programs,"folds":folds,"norms":norms,"groups":if norms{vec!["strict","stored square products","stored norm additions","stored products and additions"]}else if folds{vec!["mean loop","squared norm loop","back substitution dot","intercept loop"]}else{vec!["mean","centering","squared norms","square root","back substitution","intercept"]},"counts_all_exact_baseline_gains_losses":counts,"first_losses":first_losses,"first_gains":first_gains,"excluded":excluded,"zero_loss_improvements":safe,
        "source":"Existing cleared Householder QR, IEEE754 and Intel SDM stored PC64 arithmetic. All reflector operations, rank, pivot/tie, logarithms and statistics consumers unchanged.",
        "limits":if norms{"Four fixed norm programs: strict, stored products only, stored additions only, stored products and additions. Reflector, square root, rank, pivot and other arithmetic unchanged. Counts do not uniquely identify stores; adoption requires independent square-product evidence and classification of every loss. Unanswered questions are not scored."}else if folds{"16 fixed accumulator-retention policies: binary64 stores or residentPC64 inside each of four product/sum loops, stored at the loop result. Other outer operations remain strict. Reflector, rank, branches and published constants unchanged; no fitted parameter or selector. Unanswered questions are not scored."}else{"64 coherent group policies, no fitted parameter or data-dependent selector. Counts are diagnostic; adoption requires structural evidence and classification of every loss. Unanswered questions are not scored. Discovery only, no new Excel answers."}});
    json_write(destination.join("summary.json"),&result);println!("{}",json!({"complete":true,"rows":counts[0][0],"zero_loss_improvements":safe}));
}

/// Observe the original R input and inverse output without changing either.
#[cfg(regression_audit)]
pub fn regression_triangular(r: &[Vec<f64>], size: usize, inverse: &[Vec<f64>],
    order: &[Option<usize>], predictors: usize) {
    let matrix = |rows: &[Vec<f64>]| rows.iter().take(size)
        .map(|row| row.iter().take(size).map(|v| format!("0x{:016x}", v.to_bits())).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    REGRESSION_TRIANGULAR.with(|slot| *slot.borrow_mut() = Some(json!({
        "r": matrix(r), "inverse": matrix(inverse), "order": &order[..size], "predictors": predictors,
    })));
}

#[cfg(regression_audit)]
pub fn regression_factors(weights:&[f64],inverse:&[Vec<f64>],observations:usize,value:f64)->f64 {
    let bits=|v:f64|format!("0x{:016x}",v.to_bits());
    REGRESSION_FACTORS.with(|factors|*factors.borrow_mut()=Some(json!({"observations":observations,
        "weights":weights.iter().copied().map(bits).collect::<Vec<_>>(),
        "inverse":inverse.iter().map(|row|row.iter().copied().map(bits).collect::<Vec<_>>()).collect::<Vec<_>>()})));
    value
}

/// Observe original binary64 operands and return the unchanged original value.
#[cfg(regression_audit)]
pub fn regression_leaf(ss:f64,covariance:f64,df:f64,row:usize,column:usize,role:usize,value:f64)->f64 {
    let bits=|v:f64|format!("0x{:016x}",v.to_bits());
    let factors=if role==1 {REGRESSION_FACTORS.with(|f|f.borrow().clone())}else{None};
    let triangular = REGRESSION_TRIANGULAR.with(|f| f.borrow().clone());
    REGRESSION_LEAVES.with(|leaves|leaves.borrow_mut().push(json!({"ss":bits(ss),"covariance":bits(covariance),
        "df":bits(df),"row":row,"column":column,"role":role,"value":bits(value),"factors":factors,"triangular":triangular})));
    value
}

/// Freeze standard-error operands only after an exact registered-route control.
#[cfg(regression_audit)]
pub fn regression_stat_leaves(input:&str,out:&str) {
    use std::io::Write;
    let destination=Path::new(out);fs::create_dir(destination).unwrap();
    let rows=crate::lines::<Value>(input);
    let mut scoped=std::io::BufWriter::new(fs::File::create_new(destination.join("scoped.jsonl")).unwrap());
    let mut omitted=std::collections::BTreeMap::<String,usize>::new();
    let mut controls=0;let mut admitted=0;
    for (index,raw) in rows.iter().enumerate() {
        let name=raw["function"].as_str().unwrap();
        let selected=raw["result_index"].as_array();
        if !matches!(name,"LINEST"|"LOGEST") || selected.is_none_or(|s|s.len()!=2 || s[0]!=2 || s[1].as_u64().is_none_or(|c|c==0)) {
            *omitted.entry("outside_standard_error_selection".into()).or_default()+=1;continue;
        }
        let selected=selected.unwrap();let column=selected[1].as_u64().unwrap() as usize;
        let args:Vec<_>=raw["args"].as_array().unwrap().iter().enumerate()
            .map(|(i,v)|calc_verdict::decode_argument(v,&format!("args[{i}]")).unwrap()).collect();
        let (baseline,_)=crate::invoke_selected_control(name,args.clone(),Some((2,column)));
        assert_eq!(baseline,raw["actual"].as_str().unwrap(),"frozen registered regression {}",raw["id"]);
        REGRESSION_LEAVES.with(|leaves|leaves.borrow_mut().clear());
        REGRESSION_FACTORS.with(|factors|*factors.borrow_mut()=None);
        REGRESSION_TRIANGULAR.with(|triangular|*triangular.borrow_mut()=None);
        let observed=if name=="LINEST" {crate::ported::regression_forecast_family::eval_linest_surface(&args)}
            else {crate::ported::regression_forecast_family::eval_logest_surface(&args)};
        let value=observed.unwrap_or_else(CalcValue::error);
        let chosen=calc_verdict::select_result(&value,Some((2,column))).unwrap();
        assert_eq!(calc_verdict::encode_result(&chosen).unwrap(),baseline,"instrumented regression route {}",raw["id"]);
        controls+=1;
        let leaves=REGRESSION_LEAVES.with(|leaves|leaves.borrow().clone());
        let leaf=leaves.iter().find(|leaf|leaf["row"]==2 && leaf["column"]==column);
        let Some(leaf)=leaf else {*omitted.entry("structural_zero_error_or_exact_fit".into()).or_default()+=1;continue;};
        let decoded:Vec<_>=["ss","covariance","df","value"].into_iter().map(|key| {
            let s=leaf[key].as_str().unwrap();f64::from_bits(u64::from_str_radix(&s[2..],16).unwrap())
        }).collect();
        if !(decoded.iter().all(|v|v.is_finite()) && decoded[0]>=0.0 && decoded[1]>=0.0 && decoded[2]>0.0
            && raw["expected"].as_str().is_some_and(|s|s.starts_with("0x"))) {
            *omitted.entry("outside_finite_numeric_consumer".into()).or_default()+=1;continue;
        }
        assert_eq!(leaf["value"].as_str().unwrap(),baseline,"selected leaf mismatch");
        writeln!(scoped,"{}",json!({"trace_index":index,"raw":raw,"leaf":leaf})).unwrap();admitted+=1;
    }
    assert!(admitted>0);assert_eq!(admitted+omitted.values().sum::<usize>(),rows.len());
    let summary=json!({"complete":true,"input_rows":rows.len(),"registered_and_instrumented_controls":controls,"scoped":admitted,"excluded":omitted,
        "limits":"Final standard-error operands only, original QR/covariance and original arithmetic retained. Every selected instrumented value equals the registered invocation and frozen actual result. No new Excel answer or producer identification."});
    json_write(destination.join("summary.json"),&summary);println!("{summary}");
}

/// Research-only translation of ExcelFinancialFunctions Common.newton
/// (Apache-2.0, Luca Bolognese, copyright 2010-2013), changed to return a typed
/// failure on nonfinite evaluation. Its fixed 0.05 seed, 1e-6 difference/step and
/// count>20 termination are retained. This is the Newton phase, not its fallback.
#[cfg(bond_solver_audit)]
fn published_yield_newton(price: &crate::ported::bond_core_family::ResearchPrice, target: f64, staged: bool) -> Option<f64> {
    let mut x = 0.05_f64;
    let precision = 0.000001_f64;
    for count in 0..=21 {
        let f = price.value(x, staged).ok()? - target;
        let above = price.value(x+precision, staged).ok()? - target;
        let below = price.value(x-precision, staged).ok()? - target;
        let derivative = (above-below)/(2.0*precision);
        let next = x-f/derivative;
        if !next.is_finite() { return None; }
        if (next-x).abs() < precision { return Some(next); }
        if count > 20 { return None; }
        x = next;
    }
    None
}

#[cfg(bond_solver_audit)]
pub fn yield_newton(input: &str, out: &str) {
    use crate::ported::bond_core_family::ResearchPrice;
    let destination = Path::new(out);
    fs::create_dir(destination).unwrap();
    let rows = crate::lines::<Value>(input);
    let mut counts = [[0_usize; 5]; 2]; // pass, fail, gains, losses, failed Newton phase
    let mut output = String::new();
    for (index, row) in rows.iter().enumerate() {
        assert_eq!(row["function"], "YIELD");
        let args: Vec<_> = row["args"].as_array().unwrap().iter()
            .map(|v| calc_verdict::decode_argument(v, "YIELD research").unwrap()).collect();
        let numeric: Vec<_> = args.iter().map(|v| match v.core_ref() {
            CoreRef::Number(n) => n,
            _ => panic!("numeric YIELD arguments required"), // ubs:ignore — malformed diagnostic inputs must stop, never synthesize a measurement.
        }).collect();
        let environment = CalcEnvironment::dates(DateSystem::System1900);
        let baseline = crate::invoke_selected_control_in("YIELD", args, None, Some(&environment)).0;
        assert_eq!(baseline, row["actual"]);
        let price = ResearchPrice::new(&numeric).unwrap();
        for staged in [false, true] {
            let program = usize::from(staged);
            let prediction = if price.count() <= 1 {
                Some(baseline.clone())
            } else {
                published_yield_newton(&price, numeric[3], staged)
                    .map(|v| calc_verdict::encode_result(&CalcValue::number(v)).unwrap())
            };
            let got = prediction.clone().unwrap_or_else(|| "error:Num".to_owned());
            let matched = got == row["expected"];
            let was = baseline == row["expected"];
            counts[program][usize::from(!matched)] += 1;
            counts[program][2] += usize::from(matched && !was);
            counts[program][3] += usize::from(was && !matched);
            counts[program][4] += usize::from(prediction.is_none());
            output.push_str(&serde_json::to_string(&json!({"trace_index":index,"program":program,
                "got":got,"matched":matched,"baseline_matched":was,"phase_failed":prediction.is_none()})).unwrap());
            output.push('\n');
        }
    }
    fs::write(destination.join("predictions.jsonl"), output).unwrap();
    let result = json!({"complete":true,"observations":rows.len(),"programs":2,
        "counts_pass_fail_gains_losses_failed_phase":counts,
        "source":"ExcelFinancialFunctions97756a375769e45059c7d58fc51fc7e391bc4984 common.fs/bonds.fs; Apache-2.0",
        "limits":"Published Newton phase only, with unchanged existing product date/price bodies and two existing discount leaves. Failed phases are NUM hypotheses, not a faithful port of the reference fallback. Single-coupon controls retain the current product result. Fixed constants from the reference, no fitted solver parameter or tolerance claim."});
    json_write(destination.join("summary.json").to_str().unwrap(), &result);
    println!("{result}");
}

/// Ordinary algebraic reciprocal programs. Every program applies to all four
/// reciprocal functions, so a function or residual cannot choose its graph.
fn reciprocal(z: ParsedComplex, graph: usize) -> ParsedComplex {
    let (a, b) = (z.re, z.im);
    let (re, im) = match graph {
        0 => {
            let d = a * a + b * b;
            (a / d, -b / d)
        }
        1 => {
            let d = 1.0 / (a * a + b * b);
            (a * d, -b * d)
        }
        2 => {
            if a.abs() >= b.abs() {
                let r = b / a;
                let t = 1.0 / (a + b * r);
                (t, -r * t)
            } else {
                let r = a / b;
                let t = 1.0 / (b + a * r);
                (r * t, -t)
            }
        }
        3 => {
            if a.abs() >= b.abs() {
                let r = b / a;
                let t = (1.0 / a) / (1.0 + r * r);
                (t, -r * t)
            } else {
                let r = a / b;
                let t = (1.0 / b) / (1.0 + r * r);
                (r * t, -t)
            }
        }
        _ => unreachable!(),
    };
    ParsedComplex::new(re, im, None)
}

pub fn complex_reciprocals(input: &str, out: &str) {
    let mut counts = vec![[0usize; 4]; 4]; // matches, failures, gains, losses
    let mut first = vec![None; 4];
    let mut rows = crate::lines::<Value>(input);
    rows.push(
        json!({"id":"catalogue-prefix-metadata-A128", "function":"IMSEC", "args":["4+3i"],
        "expected":"-0.0652940278579471-0.0752249603027732i"}),
    );
    let mut examples = vec![];
    for row in &rows {
        let function = row["function"].as_str().unwrap();
        let argument = calc_verdict::decode_argument(&row["args"][0], "complex").unwrap();
        let z = match argument.core_ref() {
            CoreRef::Number(n) => ParsedComplex::new(n, 0.0, None),
            CoreRef::Text(t) => match parse_complex_text(&t.to_string_lossy()) {
                Some(z) => z,
                None => continue,
            },
            _ => continue,
        };
        let value = match function {
            "IMSEC" => z.cos(),
            "IMCSC" => z.sin(),
            "IMSECH" => z.cosh(),
            "IMCSCH" => z.sinh(),
            _ => unreachable!(),
        };
        let (baseline, _) = crate::invoke_control(function, vec![argument]);
        let expected = if row["id"] == "catalogue-prefix-metadata-A128" {
            encoded(&CalcValue::text(
                crate::value::ExcelText::from_interop_assignment(row["expected"].as_str().unwrap()),
            ))
        } else {
            row["expected"].as_str().unwrap().to_owned()
        };
        for graph in 0..4 {
            let candidate = reciprocal(value, graph);
            let got = if value.re == 0.0 && value.im == 0.0 {
                "error:Num".to_owned()
            } else {
                match format_complex_text(candidate, z.suffix.unwrap_or('i')) {
                    Ok(s) => encoded(&CalcValue::text(
                        crate::value::ExcelText::from_interop_assignment(&s),
                    )),
                    Err(e) => encoded(&CalcValue::error(e)),
                }
            };
            let matched = got == expected;
            let old = baseline == expected;
            counts[graph][usize::from(!matched)] += 1;
            counts[graph][2] += usize::from(matched && !old);
            counts[graph][3] += usize::from(!matched && old);
            if !matched && first[graph].is_none() {
                first[graph] =
                    Some(json!({"witness":row,"baseline":baseline,"got":got,"expected":expected}));
            }
            if row["id"] == "catalogue-prefix-metadata-A128" {
                examples.push(json!({"graph":graph,"re_bits":format!("{:016x}",candidate.re.to_bits()),"im_bits":format!("{:016x}",candidate.im.to_bits()),"outcome":got,"matched":matched}));
            }
        }
    }
    json_write(
        out,
        &json!({"rows":rows.len(),"counts_matches_failures_gains_losses":counts,
        "first_failures":first,"corpus_example":examples,
        "graphs":["conjugate / squared norm","conjugate * reciprocal squared norm","scaled reciprocal then multiplication","reciprocal component divided by 1+ratio squared"],
        "scope":"All parseable original scalar inputs of IMSEC/IMCSC/IMSECH/IMCSCH; existing observations only; no fitted quantities."}),
    );
    println!("reciprocal graph counts {counts:?}");
}

fn encoded(value: &CalcValue) -> String {
    calc_verdict::encode_result(value).expect("scalar audit value")
}

/// Fixed DLMF Hankel programs with universal decimal-constant quantizers.
/// This is a numerical hypothesis study, not a replacement application route.
pub fn bessel_asymptotic(input: &str, tables: &str, out: &str) {
    bessel_asymptotic_search(input, tables, out, false);
}

/// The same Hankel family with nine algebraically equal amplitude associations.
pub fn bessel_amplitudes(input: &str, tables: &str, out: &str) {
    bessel_asymptotic_search(input, tables, out, true);
}

fn bessel_asymptotic_search(input: &str, tables: &str, out: &str, amplitudes: bool) {
    use crate::dsl::{EvalModel, Op, WIDE_MODELS};
    use crate::eval::{binop, Val, CW};
    use crate::ported::excel_numeric::{Fpu, Hard, excel_cos, excel_sin};
    use std::io::Write;
    fn arithmetic(kind: usize, model: EvalModel, a: Val, b: Val) -> Val {
        binop::<Hard>(match kind { 0 => Op::Add(0,0), 1 => Op::Mul(0,0), _ => Op::Div(0,0) },model,a,b)
    }
    fn negative(a: Val) -> Val {
        arithmetic(1,EvalModel::ResidentPc64,Val::f64(-1.0),a)
    }
    fn extended_unary(a: Val, model: EvalModel, kind: usize) -> Val {
        if matches!(model,EvalModel::Strict) {
            let x=a.to_f64::<Hard>();
            if kind==0 {return Val::f64(x.sqrt());}
            let x=Val::f64(x).to_ext();
            let v=if kind==1 {Hard::sin(&x,CW)}else{Hard::cos(&x,CW)};
            return Val::f64(Hard::store(&v,CW));
        }
        let cw=if matches!(model,EvalModel::StoredPc53|EvalModel::ResidentPc53) {0x123f}else{CW};
        let x=a.to_ext();
        let v=match kind {0=>Hard::sqrt(&x,cw),1=>Hard::sin(&x,cw),_=>Hard::cos(&x,cw)};
        if matches!(model,EvalModel::StoredPc53|EvalModel::StoredPc64) {Val::f64(Hard::store(&v,CW))}else{Val::Ext(v)}
    }
    fn seed(function:&str,x:f64,n:usize,phase:f64,scale:f64,model:EvalModel,trig:usize,terms:usize,amplitude:usize)->Val {
        let a=Val::f64(x);let one=Val::f64(1.0);
        let z=arithmetic(2,model,one,a);
        let z2=arithmetic(1,model,z,z);
        let phase=arithmetic(0,model,a,Val::f64(-phase));
        let (s,c)=match trig {
            0=>(Val::f64(excel_sin(phase.to_f64::<Hard>())),Val::f64(excel_cos(phase.to_f64::<Hard>()))),
            1=>(extended_unary(phase,model,1),extended_unary(phase,model,2)),
            _=>(Val::f64(libm::sin(phase.to_f64::<Hard>())),Val::f64(libm::cos(phase.to_f64::<Hard>()))),
        };
        // Integer-order DLMF a1,a2,a3 are exact dyadic rationals here.
        let (a1,a2,a3)=if n==0 {(-0.125,9.0/128.0,-75.0/1024.0)}else{(0.375,-15.0/128.0,105.0/1024.0)};
        let p=if terms>=2 {arithmetic(0,model,one,negative(arithmetic(1,model,Val::f64(a2),z2)))}else{one};
        let q=if terms==0 {Val::f64(0.0)}else {
            let coefficient=if terms>=2 {arithmetic(0,model,Val::f64(a1),negative(arithmetic(1,model,Val::f64(a3),z2)))}else{Val::f64(a1)};
            arithmetic(1,model,coefficient,z)
        };
        let body=if function=="BESSELJ" {
            arithmetic(0,model,arithmetic(1,model,c,p),negative(arithmetic(1,model,s,q)))
        }else{arithmetic(0,model,arithmetic(1,model,s,p),arithmetic(1,model,c,q))};
        let c=Val::f64(scale);
        match amplitude {
            0=>arithmetic(1,model,extended_unary(arithmetic(2,model,c,a),model,0),body),
            1=>arithmetic(1,model,arithmetic(2,model,extended_unary(c,model,0),extended_unary(a,model,0)),body),
            2=>arithmetic(1,model,extended_unary(c,model,0),arithmetic(2,model,body,extended_unary(a,model,0))),
            3=>arithmetic(2,model,arithmetic(1,model,extended_unary(c,model,0),body),extended_unary(a,model,0)),
            4=>arithmetic(1,model,arithmetic(1,model,extended_unary(c,model,0),extended_unary(z,model,0)),body),
            5=>arithmetic(1,model,extended_unary(c,model,0),arithmetic(1,model,extended_unary(z,model,0),body)),
            6=>arithmetic(1,model,arithmetic(1,model,extended_unary(c,model,0),body),extended_unary(z,model,0)),
            7=>arithmetic(2,model,body,extended_unary(arithmetic(2,model,a,c),model,0)),
            8=>arithmetic(1,model,arithmetic(2,model,one,extended_unary(arithmetic(2,model,a,c),model,0)),body),
            _=>unreachable!(),
        }
    }
    let destination=Path::new(out);fs::create_dir(destination).unwrap();
    let original=crate::lines::<Value>(input);
    let tables:Value=serde_json::from_slice(&fs::read(tables).unwrap()).unwrap();
    let phases=tables["phases"].as_array().unwrap();let scales=tables["scales"].as_array().unwrap();
    assert_eq!((phases.len(),scales.len()),(10,10));
    let decode=|v:&Value| f64::from_bits(u64::from_str_radix(v.as_str().unwrap().trim_start_matches("0x"),16).unwrap());
    let mut counts=Vec::new();let mut programs=Vec::new();
    for phase in 0..10 {for scale in 0..10 {for body in 0..5 {for recurrence in 0..5 {for trig in 0..3 {for terms in 0..3 {for association in 0..3 {for amplitude in 0..if amplitudes {9}else{1} {
        programs.push([phase,scale,body,recurrence,trig,terms,association,amplitude]);
    }}}}}}}}
    let mut axes=json!([["phase",10],["scale",10],["body",5],["recurrence",5],["trig",3],["terms",3],["association",3]]);
    if amplitudes {axes.as_array_mut().unwrap().push(json!(["amplitude",9]));}
    let mut grammar=json!({
        "programs":programs.len(),"axes_in_loop_order":axes,
        "id":"Mixed-radix rank, association fastest. Model IDs follow WIDE_MODELS. Trig0=current worksheet kernel,1=native instruction,2=libm.",
        "terms":"0 leading sin/cos;1 includes a1/x;2 includes a2/x^2 and a3/x^3."});
    if amplitudes {
        grammar["id"]=json!("Mixed-radix rank, amplitude fastest. Other axes retain the original Hankel meanings; model IDs follow WIDE_MODELS.");
        grammar["amplitudes"]=json!(["sqrt(c/x)*b","(sqrt(c)/sqrt(x))*b","sqrt(c)*(b/sqrt(x))","(sqrt(c)*b)/sqrt(x)","(sqrt(c)*sqrt(1/x))*b","sqrt(c)*(sqrt(1/x)*b)","(sqrt(c)*b)*sqrt(1/x)","b/sqrt(x/c)","(1/sqrt(x/c))*b"]);
    }
    json_write(destination.join("grammar.json").to_str().unwrap(),&grammar);
    for function in ["BESSELJ","BESSELY"] {
        let rows:Vec<_>=original.iter().enumerate().filter_map(|(i,r)| {
            if r["function"]!=function||!r["result_index"].is_null()||!r["expected"].as_str()?.starts_with("0x") {return None;}
            let args=r["args"].as_array()?;if args.len()!=2||!args.iter().all(|a|a.as_str().is_some_and(|s|s.starts_with("0x"))) {return None;}
            let x=decode(&args[0]);let n=decode(&args[1]);
            (x.abs()>=2_f64.powi(24)&&x.is_finite()&&n==n.trunc()&&(0.0..=4.0).contains(&n)&&(function=="BESSELJ"||x>0.0))
                .then_some((i,x,n as usize,r))
        }).collect();
        assert!(!rows.is_empty());
        let mut first=std::io::BufWriter::new(fs::OpenOptions::new().write(true).create_new(true).open(destination.join(format!("{function}-first-kills.tsv"))).unwrap());
        writeln!(first,"program\ttrace_index\tgot").unwrap();
        let mut alive=Vec::new();let mut executions=0;
        for (id,p) in programs.iter().enumerate() {
            let index=|key:&str|match key {"phase"=>p[0],"scale"=>p[1],"body"=>p[2],"recurrence"=>p[3],"trig"=>p[4],"terms"=>p[5],"association"=>p[6],"amplitude"=>p[7],"id"=>id,_=>unreachable!()};
            let phase=&phases[index("phase")]["values"];
            let phase0=decode(&phase[0]["binary64"]);let phase1=decode(&phase[1]["binary64"]);
            let scale=decode(&scales[index("scale")]["value"]["binary64"]);
            let body=WIDE_MODELS[index("body")];let model=WIDE_MODELS[index("recurrence")];
            let mut killed=false;
            for &(at,x,n,row) in &rows {
                let ax=x.abs();
                let mut previous=seed(function,ax,0,phase0,scale,body,index("trig"),index("terms"),index("amplitude"));
                let mut current=seed(function,ax,1,phase1,scale,body,index("trig"),index("terms"),index("amplitude"));
                let tox=arithmetic(2,model,Val::f64(2.0),Val::f64(ax));
                for k in 1..n {
                    let step=match index("association") {
                        0=>arithmetic(1,model,arithmetic(1,model,Val::f64(k as f64),tox),current),
                        1=>arithmetic(1,model,arithmetic(2,model,Val::f64(2.0*k as f64),Val::f64(ax)),current),
                        _=>arithmetic(2,model,arithmetic(1,model,Val::f64(2.0*k as f64),current),Val::f64(ax)),
                    };
                    let next=arithmetic(0,model,step,negative(previous));previous=current;current=next;
                }
                let result=if n==0 {previous}else{current};
                let mut got=result.to_f64::<Hard>();if x<0.0&&n%2==1 {got=-got;}
                let encoded=format!("0x{:016x}",got.to_bits());executions+=1;
                if row["expected"]!=encoded {
                    writeln!(first,"{}\t{at}\t{encoded}",index("id")).unwrap();killed=true;break;
                }
            }
            if !killed {alive.push(index("id"));}
        }
        json_write(destination.join(format!("{function}-survivors.json")).to_str().unwrap(),&json!(alive));
        counts.push(json!({"function":function,"programs":programs.len(),"observations":rows.len(),"executions":executions,"survivors":alive.len()}));
    }
    let result=json!({"complete":true,"counts":counts,"sources":"DLMF10.17.1--4 and generated universal-constant quantizers; no Excel output used as a component.",
        "limits":"Numeric observations only, |x|>=2^24, integer orders0..4. Exact input/answer traces retained; domain outcomes excluded as separate obligations. Five coherent arithmetic/store models; current worksheet trig, native instruction and libm trig leaves. No candidate is an Excel kernel unless every obligation in its region is satisfied; this study alone cannot promote a body."});
    json_write(destination.join("summary.json").to_str().unwrap(),&result);println!("{result}");
}

/// Freeze ordinary Workbook predictions for one question sheet. No expected
/// answer is read or synthesized; these predictions never enter an Excel bank.
pub fn sheet_predictions(input: &str, out: &str) {
    use std::io::Write;
    let doc: Value = serde_json::from_slice(&fs::read(input).unwrap()).unwrap();
    let mut book = Workbook::new_blank();
    let sheet = book.sheets()[0].id;
    let questions = doc["questions"].as_array().unwrap();
    for question in questions {
        assert!(question.get("expected").is_none(), "a question contains an answer");
        let at = CellKey::new(sheet, question["row"].as_u64().unwrap() as u32,
            question["column"].as_u64().unwrap() as u32).unwrap();
        assert!(book.enter(at, &format!("={}", question["formula"].as_str().unwrap())).accepted());
    }
    let mut output = std::io::BufWriter::new(fs::OpenOptions::new().write(true).create_new(true).open(out).unwrap());
    for question in questions {
        let at = CellKey::new(sheet, question["row"].as_u64().unwrap() as u32,
            question["column"].as_u64().unwrap() as u32).unwrap();
        let cell = book.cell(at).unwrap();
        let (prediction, provenance) = match cell.value() {
            ProjectedValue::Calculated { value, .. } => (encoded(value), "calculated"),
            ProjectedValue::NotModeled { refusal, .. } if doc["allow_refusals"] == true => {
                (format!("fault:{refusal:?}"), "not_modeled")
            }
            _ => panic!("question prediction was not calculated"), // ubs:ignore — missing provenance aborts; only explicitly requested typed refusals can be frozen
        };
        writeln!(output, "{}", json!({"id":question["id"],"formula":question["formula"],
            "row":question["row"],"column":question["column"],"calc_func_prediction":prediction,
            "prediction_provenance":provenance,
            "not_excel_evidence":true,"route":"ordinary Workbook::enter on the whole question sheet"})).unwrap();
    }
}

#[cfg(test)]
mod question_tests {
    use super::*;

    #[cfg(regression_outer_audit)]
    #[test]
    fn unanswered_qr_questions_have_no_scored_expectation() {
        assert_eq!(regression_expectation(&json!({"expected":null,"unobserved_question":true})),None);
        assert_eq!(regression_expectation(&json!({"expected":"0x3ff0000000000000"})),Some("0x3ff0000000000000"));
    }

    #[cfg(regression_outer_audit)]
    #[test]
    #[should_panic(expected="must not discard a recorded expectation")]
    fn a_question_flag_cannot_hide_a_recorded_qr_answer() {
        regression_expectation(&json!({"expected":"0x3ff0000000000000","unobserved_question":true}));
    }

    /// Native-only: the research command writes its prediction files.
    #[test]
    fn question_sheet_keeps_refusal_provenance_and_strict_default() {
        let folder = std::env::temp_dir().join(format!("dnacalc-question-test-{}-{}",
            std::process::id(), std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()));
        fs::create_dir(&folder).unwrap();
        let input = folder.join("questions.json");
        let output = folder.join("predictions.jsonl");
        // Entry validation now rejects missing IRR arguments. Preserve that
        // check and exercise an admitted hyperlink with an unsupported friendly name.
        let mut admission = Workbook::new_blank();
        let admission_sheet = admission.sheets()[0].id;
        let invalid = CellKey::new(admission_sheet, 1, 1).unwrap();
        assert!(!admission.enter(invalid, "=IRR()").accepted());
        let mut document = json!({"questions":[
            {"id":"synthetic-calculated","row":1,"column":4,"formula":"1+1"},
            {"id":"synthetic-refusal","row":2,"column":4,"formula":"HYPERLINK(\"#A1\",LAMBDA(x,x))"}
        ],"allow_refusals":true});
        fs::write(&input, document.to_string()).unwrap();
        sheet_predictions(input.to_str().unwrap(), output.to_str().unwrap());
        let rows = crate::lines::<Value>(output.to_str().unwrap());
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0]["calc_func_prediction"], "0x4000000000000000");
        assert_eq!(rows[0]["prediction_provenance"], "calculated");
        assert_eq!(rows[1]["prediction_provenance"], "not_modeled");
        assert!(rows[1]["calc_func_prediction"].as_str().unwrap().starts_with("fault:"));
        document.as_object_mut().unwrap().remove("allow_refusals");
        fs::write(&input, document.to_string()).unwrap();
        assert!(std::panic::catch_unwind(|| sheet_predictions(input.to_str().unwrap(),
            folder.join("strict.jsonl").to_str().unwrap())).is_err());
        for name in ["questions.json", "predictions.jsonl", "strict.jsonl"] {
            fs::remove_file(folder.join(name)).unwrap();
        }
        fs::remove_dir(folder).unwrap();
    }
}

/// Check literal operands and whole calls through the ordinary workbook route.
/// This extends the existing audit; it never supplies replacement Excel answers.
pub fn function_routes(input: &str, out: &str) {
    use std::io::Write;
    let mut output = std::io::BufWriter::new(
        fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(out)
            .unwrap(),
    );
    let rows = crate::lines::<Value>(input);
    for raw in &rows {
        let mut book = Workbook::new_blank();
        let sheet = book.sheets()[0].id;
        let at = |row, column| CellKey::new(sheet, row, column).unwrap();
        let mut echoes = vec![];
        let numeric_cells = raw["numeric_cell_inputs"].as_bool() == Some(true);
        let formula;
        if numeric_cells {
            assert!(
                raw["result_index"].is_null(),
                "numeric cell audit requires a scalar result"
            );
            let arguments = raw["args"].as_array().unwrap();
            let mut references = vec![];
            for (index, encoded_argument) in arguments.iter().enumerate() {
                let argument =
                    calc_verdict::decode_argument(encoded_argument, "numeric cell input").unwrap();
                assert!(
                    matches!(argument.core_ref(), CoreRef::Number(_)),
                    "explicit numeric cells required"
                );
                let row = 2 + 40 * index as u32;
                assert!(
                    book.dispatch(SkinIntent::Edit {
                        edits: vec![CellEdit::Set {
                            at: at(row, 1),
                            input: literal(&argument),
                        }]
                    })
                    .accepted()
                );
                let cell = book.cell(at(row, 1)).unwrap();
                let ProjectedValue::Literal(value) = cell.value() else {
                    panic!("authored numeric input changed kind") // ubs:ignore — changed inputs cannot establish a calculation-route comparison
                };
                assert_eq!(
                    encoded(value),
                    encoded(&argument),
                    "numeric cell input changed"
                );
                echoes.push(encoded(value));
                references.push(format!("A{row}"));
            }
            formula = format!(
                "={}({})",
                raw["function"].as_str().unwrap(),
                references.join(",")
            );
        } else {
            formula = raw["formula"].as_str().unwrap().to_owned();
            for (index, source) in raw["argument_sources"]
                .as_array()
                .unwrap()
                .iter()
                .enumerate()
            {
                let argument =
                    calc_verdict::decode_argument(&raw["args"][index], "original argument")
                        .unwrap();
                let top = 2 + 40 * index as u32;
                assert!(
                    book.enter(at(top, 1), &format!("={}", source.as_str().unwrap()))
                        .accepted()
                );
                if let Some(array) = argument.as_array() {
                    for r in 0..array.shape().rows {
                        for c in 0..array.shape().cols {
                            let cell = book.cell(at(top + r as u32, 1 + c as u32)).unwrap();
                            let ProjectedValue::Calculated { value, .. } = cell.value() else {
                                panic!("argument echo was not calculated") // ubs:ignore — fail closed on a missing parser observation
                            };
                            let wanted = encoded(array.get(r, c).unwrap());
                            assert_eq!(encoded(value), wanted, "operand changed: {}", raw["id"]);
                            echoes.push(wanted);
                        }
                    }
                } else {
                    let cell = book.cell(at(top, 1)).unwrap();
                    let ProjectedValue::Calculated { value, .. } = cell.value() else {
                        panic!("scalar echo was not calculated") // ubs:ignore — cached values cannot establish parser parity
                    };
                    assert_eq!(
                        encoded(value),
                        encoded(&argument),
                        "operand changed: {}",
                        raw["id"]
                    );
                    echoes.push(encoded(value));
                }
            }
        }
        let args = raw["args"]
            .as_array()
            .unwrap()
            .iter()
            .map(|a| calc_verdict::decode_argument(a, "original argument").unwrap())
            .collect();
        let selection = raw["result_index"].as_array().map(|v| {
            (
                v[0].as_u64().unwrap() as usize,
                v[1].as_u64().unwrap() as usize,
            )
        });
        let environment = (raw["date_system_1900"].as_bool() == Some(true))
            .then(|| CalcEnvironment::dates(DateSystem::System1900));
        let (registered, _) = crate::invoke_selected_control_in(
            raw["function"].as_str().unwrap(), args, selection, environment.as_ref());
        let result = at(1, 30);
        assert!(book.enter(result, &formula).accepted());
        let cell = book.cell(result).unwrap();
        let ProjectedValue::Calculated { value, .. } = cell.value() else {
            panic!("function call was not calculated") // ubs:ignore — failed route cannot count as a semantic observation
        };
        let workbook = encoded(value);
        assert_eq!(
            workbook, registered,
            "workbook/registered discrepancy: {}",
            raw["id"]
        );
        assert_eq!(
            workbook, raw["actual"],
            "current/frozen product discrepancy: {}",
            raw["id"]
        );
        writeln!(output, "{}", json!({"raw":raw,"workbook":workbook,"registered":registered,"formula":formula,
            "argument_transport":if numeric_cells {"authored numeric cells"}else {"formula literals"},
            "echoes":echoes,"new_excel_readings":0,
            "limit":"Checks the current application's literal and typed routes. Original Excel effective-input transport remains an assumption."})).unwrap();
    }
    println!("{} complete function-route and operand audits", rows.len());
}

/// Verify the numeric operand of each conversion observation through the
/// workbook, before using it as a digit-extraction constraint.
pub fn number_inputs(input: &str, out: &str) {
    let rows: Vec<Value> = serde_json::from_slice(&fs::read(input).unwrap()).unwrap();
    let mut book = Workbook::new_blank();
    let sheet = book.sheets()[0].id;
    let at = |row, column| CellKey::new(sheet, row, column).unwrap();
    let output = at(1, 6);
    let mut records = vec![];
    for row in &rows {
        let bits = u64::from_str_radix(row["input_bits"].as_str().unwrap(), 16).unwrap();
        let wanted = f64::from_bits(bits);
        let source = row["formula"]
            .as_str()
            .unwrap()
            .strip_suffix("&\"\"")
            .unwrap();
        let reference = source.strip_prefix('A').and_then(|s| s.parse::<u32>().ok());
        if let Some(r) = reference {
            assert!(
                book.dispatch(SkinIntent::Edit {
                    edits: vec![CellEdit::Set {
                        at: at(r, 1),
                        input: CellInput::Number(FiniteNumber::new(wanted).unwrap())
                    }]
                })
                .accepted()
            );
        }
        assert!(book.enter(output, &format!("={source}")).accepted());
        let cell = book.cell(output).unwrap();
        let ProjectedValue::Calculated { value, .. } = cell.value() else {
            panic!("number operand was not calculated") // ubs:ignore — refuse an audit that did not run the engine
        };
        let actual = value.as_number().expect("a numeric operand");
        assert!(book.enter(output, &format!("={source}&\"\"")).accepted());
        let text_cell = book.cell(output).unwrap();
        let ProjectedValue::Calculated { value, .. } = text_cell.value() else {
            panic!("text was not calculated") // ubs:ignore — refuse an audit that did not calculate the conversion
        };
        let text = value.as_text().unwrap().to_string_lossy();
        records.push(json!({"source":row,"admitted_bits":format!("{:016x}",actual.to_bits()),
            "input_kept":actual.to_bits()==bits,"calculated_text":text,
            "matches_original_text":text==row["expected_text"].as_str().unwrap(),
            "binding":if reference.is_some(){"exact authored numeric input"}else{"actual formula arithmetic"}}));
    }
    let changed = records.iter().filter(|r| r["input_kept"] != true).count();
    let summary = json!({"rows":rows.len(),"unchanged_numeric_operands":rows.len()-changed,"changed_operands":changed,
        "matches_original_text":records.iter().filter(|r|r["matches_original_text"]==true).count(),
        "scope":"Existing conversion observations, inspected through actual workbook calculation; no new Excel reading."});
    fs::create_dir(out).unwrap();
    json_write(Path::new(out).join("records.json"), &records);
    json_write(Path::new(out).join("summary.json"), &summary);
    println!("{summary}");
    assert_eq!(
        changed, 0,
        "conversion-program constraints used an unverified input"
    );
}

fn literal(value: &CalcValue) -> CellInput {
    match value.core_ref() {
        CoreRef::Number(n) => CellInput::Number(FiniteNumber::new(n).unwrap()),
        CoreRef::Logical(b) => CellInput::Logical(b),
        CoreRef::Text(t) => CellInput::Text(CellText::new(t.clone()).unwrap()),
        _ => panic!("not a scalar literal"), // ubs:ignore — the frozen audit admits scalar literals only
    }
}

/// Compare the original literal surface and echo cells without treating an
/// echo as proof of Excel's unobservable effective callee arguments.
pub fn literal_inputs(request: &str, observation: &str, out: &str) {
    let submitted: Value = serde_json::from_slice(&fs::read(request).unwrap()).unwrap();
    let observed: Value = serde_json::from_slice(&fs::read(observation).unwrap()).unwrap();
    let mut book = Workbook::new_blank();
    let sheet = book.sheets()[0].id;
    let at = |column| CellKey::new(sheet, 1, column).unwrap();
    let mut records = vec![];
    for row in submitted["rows"].as_array().unwrap() {
        let Some(literals) = row.get("argument_literals").and_then(Value::as_array) else {
            continue;
        };
        let reading = observed["rows"]
            .as_array()
            .unwrap()
            .iter()
            .find(|r| r["id"] == row["id"])
            .unwrap();
        assert_eq!(
            row["argument_literals"],
            reading["input_audit"]["argument_literals"]
        );
        assert_eq!(row["formula"], reading["input_audit"]["formula"]);
        let mut echoes = vec![];
        for (i, literal) in literals.iter().enumerate() {
            assert!(
                book.enter(at(i as u32 + 1), &format!("={}", literal.as_str().unwrap()))
                    .accepted()
            );
            let cell = book.cell(at(i as u32 + 1)).unwrap();
            let ProjectedValue::Calculated { value, .. } = cell.value() else {
                panic!("literal echo did not calculate") // ubs:ignore — failed calculation cannot supply a parser observation
            };
            echoes.push(encoded(value));
        }
        let mut outputs = vec![];
        for formula in [row["formula"].as_str().unwrap(), "=PMT(A1,B1,C1,D1,E1)"] {
            assert!(book.enter(at(6), formula).accepted());
            let cell = book.cell(at(6)).unwrap();
            let ProjectedValue::Calculated { value, .. } = cell.value() else {
                panic!("literal route did not calculate") // ubs:ignore — cached outputs cannot supply an arithmetic observation
            };
            outputs.push(encoded(value));
        }
        let excel_echoes: Vec<_> = reading["input_audit"]["read_back"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| format!("0x{}", v["number"].as_str().unwrap()))
            .collect();
        let excel_result = format!("0x{}", reading["value"]["number"].as_str().unwrap());
        records.push(json!({"id":row["id"],"request":row,"original_observation":reading,
            "application_echoes":echoes,"excel_echoes":excel_echoes,"echoes_agree":echoes==excel_echoes,
            "literal_output":outputs[0],"reference_output":outputs[1],"routes_agree":outputs[0]==outputs[1],
            "matches_recorded_literal_result":outputs[0]==excel_result,
            "limits":"Literal-surface audit only. Echoes do not establish Excel's effective direct-callee operands; changed original inputs are not semantic matches."}));
    }
    let summary = json!({"rows":records.len(),"matching_echoes":records.iter().filter(|r|r["echoes_agree"]==true).count(),
        "application_literal_reference_agreement":records.iter().filter(|r|r["routes_agree"]==true).count(),
        "raw_literal_output_coincidences":records.iter().filter(|r|r["matches_recorded_literal_result"]==true).count(),
        "new_excel_answers":0,"limits":"Original literal audits only; no changed-input output counts as a same-input semantic match."});
    fs::create_dir(out).unwrap();
    json_write(Path::new(out).join("records.json"), &records);
    json_write(Path::new(out).join("summary.json"), &summary);
    println!("{summary}");
}

pub fn run(input: &str, out: &str) {
    let doc: Value = serde_json::from_slice(&fs::read(input).unwrap()).unwrap();
    assert_eq!(doc["function"], "PMT");
    let mut book = Workbook::new_blank();
    let sheet = book.sheets()[0].id;
    let at = |column| CellKey::new(sheet, 1, column).unwrap();
    let environment = CalcEnvironment {
        profile: book.profile().cloned(),
        ..CalcEnvironment::dates(DateSystem::System1900)
    };
    let mode = NumericMode::Fast;
    let target = catalogue().resolve_surface(Notation::Call, "PMT").unwrap();
    let mut records = Vec::new();
    for row in doc["witnesses"].as_array().unwrap() {
        // Literal transport audits are deliberately not semantic observations.
        if row.get("input_transport").and_then(Value::as_str) == Some("formula_literals") {
            continue;
        }
        let args: Vec<_> = row["args"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| calc_verdict::decode_argument(v, "arg").unwrap())
            .collect();
        assert_eq!(args.len(), 5);
        let edit = book.dispatch(SkinIntent::Edit {
            edits: args
                .iter()
                .enumerate()
                .map(|(i, v)| CellEdit::Set {
                    at: at(i as u32 + 1),
                    input: literal(v),
                })
                .collect(),
        });
        assert!(edit.accepted(), "{:?}", edit.error);
        assert!(book.enter(at(6), "=PMT(A1,B1,C1,D1,E1)").accepted());
        let echoes: Vec<_> = (1..=5)
            .map(|column| {
                let cell = book.cell(at(column)).unwrap();
                let ProjectedValue::Literal(value) = cell.value() else {
                    panic!("literal changed kind") // ubs:ignore — changed input kinds cannot support a route comparison
                };
                encoded(value)
            })
            .collect();
        assert_eq!(
            echoes,
            args.iter().map(encoded).collect::<Vec<_>>(),
            "authored inputs changed"
        );
        let result = book.cell(at(6)).unwrap();
        let ProjectedValue::Calculated { value, .. } = result.value() else {
            panic!("not calculated: {:?}", result.value()) // ubs:ignore — refuse cached, stale or missing audit outputs
        };
        let worksheet = encoded(value);
        let registered = invoke(
            target,
            &mut Given(args.clone()),
            &mut EvalServices::new(&mode).with_environment(&environment),
            &mut FunctionCallScratch::new(),
        )
        .map(|v| encoded(&v))
        .unwrap();
        assert_eq!(
            worksheet, registered,
            "worksheet route differs from registered call"
        );
        let admitted: Vec<_> = args
            .iter()
            .map(|v| match v.core_ref() {
                CoreRef::Number(n) => n,
                CoreRef::Logical(b) => {
                    if b {
                        1.0
                    } else {
                        0.0
                    }
                }
                CoreRef::Text(t) => calc_func::read_number_text(t.utf16_code_units(), &environment)
                    .unwrap()
                    .unwrap(),
                _ => panic!("unsupported PMT audit argument"), // ubs:ignore — frozen audit argument types are checked before comparison
            })
            .collect();
        let timing = if admitted[4].trunc() == 0.0 {
            tv::PaymentTiming::EndOfPeriod
        } else {
            tv::PaymentTiming::BeginningOfPeriod
        };
        let direct = tv::pmt(admitted[0], admitted[1], admitted[2], admitted[3], timing).unwrap();
        assert_eq!(
            registered,
            encoded(&CalcValue::number(direct)),
            "adapter differs from direct cleared body"
        );
        let (saved, _) = book.save().unwrap();
        let reopened = Workbook::open("route-audit.xlsx", saved).unwrap();
        let re_sheet = reopened.sheets()[0].id;
        let re_at = |column| CellKey::new(re_sheet, 1, column).unwrap();
        for (i, wanted) in echoes.iter().enumerate() {
            let cell = reopened.cell(re_at(i as u32 + 1)).unwrap();
            let ProjectedValue::Literal(v) = cell.value() else {
                panic!("saved input changed kind") // ubs:ignore — a changed saved input invalidates the audit
            };
            assert_eq!(&encoded(v), wanted, "save/reopen changed authored bits");
        }
        let re_result = reopened.cell(re_at(6)).unwrap();
        let ProjectedValue::Calculated { value, .. } = re_result.value() else {
            panic!("reopen did not calculate") // ubs:ignore — a file cache cannot stand in for a calculated result
        };
        assert_eq!(
            encoded(value),
            worksheet,
            "save/reopen changed calculated result"
        );
        let variant = if row["args"]
            .as_array()
            .unwrap()
            .iter()
            .any(Value::is_boolean)
        {
            "logical_timing"
        } else if row["args"]
            .as_array()
            .unwrap()
            .iter()
            .any(|x| x.as_str().is_some_and(|s| !s.starts_with("0x")))
        {
            "text_rate"
        } else {
            "numeric"
        };
        records.push(json!({"id":row["id"], "args":row["args"], "expected":row["expected_bits"],
            "authored_echoes":echoes, "admitted_numbers":admitted.iter().map(|&v|encoded(&CalcValue::number(v))).collect::<Vec<_>>(),
            "canonical_timing": if timing==tv::PaymentTiming::EndOfPeriod {0} else {1},
            "workbook":worksheet,"registered":registered,"direct_body":encoded(&CalcValue::number(direct)),
            "saved_reopened":encoded(value),"variant":variant,
            "matches_excel":worksheet==row["expected_bits"].as_str().unwrap(),
            "provenance":"engine Calculated; original expected value retained separately"}));
    }
    fs::create_dir_all(out).unwrap();
    json_write(Path::new(out).join("route-records.json"), &records);
    let mut cohorts = std::collections::BTreeMap::<String, [usize; 2]>::new();
    for r in &records {
        cohorts
            .entry(r["variant"].as_str().unwrap().into())
            .or_default()[usize::from(r["matches_excel"] != true)] += 1;
    }
    let summary = json!({"calls":records.len(), "unchanged_literal_echoes":records.len(),
        "worksheet_registered_direct_agreement":records.len(),
        "save_reopen_preserves_inputs_and_calculated_results":records.len(),
        "cohorts_matches_failures":cohorts,
        "matches_excel":records.iter().filter(|r|r["matches_excel"]==true).count(),
        "historical_route":"historical predictions are recorded Excel readings; no historical application executable was supplied",
        "finding":"For these original inputs the typed adapter and worksheet references do not explain the mismatch. The current cleared arithmetic body already returns the differing bits. Text rates use the workbook locale here.",
        "limitations":"This is a route audit on existing observations, not identification of Excel's internal operations."});
    json_write(Path::new(out).join("summary.json"), &summary);
    println!("{summary}");
}
