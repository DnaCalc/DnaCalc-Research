//! Fixed Loader binomial arithmetic programs, using the current cleared
//! Stirling table and near-mean series. No fitted constants or row selectors.
use crate::{dsl::{EvalModel, Op}, eval::{Val, binop, CW}, ported::{discrete_dist as current, excel_numeric::*}};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::{fs, io::{BufRead, BufReader, BufWriter, Write}, path::Path};

#[derive(Clone, Debug, Serialize, Deserialize)]
struct Profile {
    id: usize,
    direct_model: EvalModel,
    log: usize,
    association: usize,
    store_ratio: bool,
    correction_model: EvalModel,
    store_correction: bool,
    normalizer_association: usize,
    normalizer_log: usize,
}

fn profiles() -> Vec<Profile> {
    let mut ps = Vec::new();
    for direct_model in crate::dsl::WIDE_MODELS {
        for log in 0..3 {
            for association in 0..5 {
                for store_ratio in [false, true] {
                    for correction_model in crate::dsl::WIDE_MODELS {
                        for store_correction in [false, true] {
                            for normalizer_association in 0..3 {
                                for normalizer_log in 0..2 {
                                    ps.push(Profile {id:ps.len(), direct_model, log, association,
                                        store_ratio, correction_model, store_correction,
                                        normalizer_association, normalizer_log});
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    ps
}

fn negative(v: Val) -> Val {
    match v { Val::F64(b) => Val::F64(b ^ (1<<63)), Val::Ext(mut e) => {e.0[9] ^= 0x80; Val::Ext(e)} }
}
fn add<F:Fpu>(m:EvalModel,a:Val,b:Val)->Val {binop::<F>(Op::Add(0,1),m,a,b)}
fn sub<F:Fpu>(m:EvalModel,a:Val,b:Val)->Val {add::<F>(m,a,negative(b))}
fn mul<F:Fpu>(m:EvalModel,a:Val,b:Val)->Val {binop::<F>(Op::Mul(0,1),m,a,b)}
fn div<F:Fpu>(m:EvalModel,a:Val,b:Val)->Val {binop::<F>(Op::Div(0,1),m,a,b)}

fn deviance<F:Fpu>(x:f64,mean:f64,p:&Profile)->Val {
    let difference=(x-mean).abs();
    let boundary=0.1*(x+mean);
    if difference<=boundary {
        return Val::f64(current::deviance(x,mean));
    }
    let m=p.direct_model;
    let a=Val::f64(x);let b=Val::f64(mean);
    let mut ratio=div::<F>(m,a,b);
    if p.store_ratio {ratio=Val::f64(ratio.to_f64::<F>());}
    let log=match p.log {
        0=>Val::f64(excel_log(ratio.to_f64::<F>())),
        1=>Val::f64(crt_log(ratio.to_f64::<F>())),
        2=>Val::Ext(research_log::<F>(&ratio.to_ext())),
        _=>unreachable!(),
    };
    match p.association {
        0=>sub::<F>(m,add::<F>(m,mul::<F>(m,a,log),b),a),
        1=>add::<F>(m,mul::<F>(m,a,log),sub::<F>(m,b,a)),
        2=>add::<F>(m,sub::<F>(m,mul::<F>(m,a,log),a),b),
        3=>add::<F>(m,mul::<F>(m,a,sub::<F>(m,log,Val::f64(1.0))),b),
        4=>sub::<F>(m,mul::<F>(m,a,log),sub::<F>(m,a,b)),
        _=>unreachable!(),
    }
}

fn mass<F:Fpu>(x:f64,n:f64,prob:f64,p:&Profile)->String {
    let m=p.correction_model;
    let rest=n-x;let q=1.0-prob;
    let st=|v| Val::f64(current::stirling_error(v));
    let lc=sub::<F>(m,sub::<F>(m,sub::<F>(m,sub::<F>(m,st(n),st(x)),st(rest)),
        deviance::<F>(x,n*prob,p)),deviance::<F>(rest,n*q,p));
    let lc=if p.store_correction {Val::f64(lc.to_f64::<F>())}else{lc};
    let tau=2.0*std::f64::consts::PI;
    let z=match p.normalizer_association {
        0=>tau*x*rest/n,
        1=>tau*(x*rest/n),
        2=>tau*x*(rest/n),
        _=>unreachable!(),
    };
    let lf=if p.normalizer_log==0 {crt_log(z)}else{excel_log(z)};
    let half=F::div(&research_load(lf),&research_load(2.0),CW);
    let exponent=F::sub(&lc.to_ext(),&half,CW);
    let output=F::store(&research_exp::<F>(&exponent),CW);
    if !output.is_finite() {"error:Num".into()}
    else if output.is_subnormal() || output==0.0 {"0x0000000000000000".into()}
    else {format!("0x{:016x}",output.to_bits())}
}

fn numeric(v:&Value)->Option<f64> {
    let s=v.as_str()?;Some(f64::from_bits(u64::from_str_radix(s.strip_prefix("0x")?,16).ok()?))
}

fn coordinates(row:&Value)->Option<(f64,f64,f64)> {
    if !matches!(row["function"].as_str(),Some("BINOM.DIST"|"BINOMDIST")) || !row["result_index"].is_null() {return None;}
    let a=row["args"].as_array()?;
    if a.len()!=4 || !(a[3]==false || numeric(&a[3])==Some(0.0)) {return None;}
    let (x,n,p)=(numeric(&a[0])?,numeric(&a[1])?,numeric(&a[2])?);
    if !(x.is_finite() && n.is_finite() && 0.0<x && x<n && x==x.trunc() && n==n.trunc() && n<=1_000_000.0 && p>0.0 && p<1.0) {return None;}
    let direct=|v:f64,mean:f64| (v-mean).abs()>0.1*(v+mean);
    if !(direct(x,n*p)||direct(n-x,n*(1.0-p))) {return None;}
    if !row["expected"].as_str()?.starts_with("0x") {return None;}
    Some((x,n,p))
}

pub fn run<F:Fpu>(input:&str,output:&str) {
    let destination=Path::new(output);fs::create_dir(destination).unwrap();
    let original:Vec<Value>=BufReader::new(fs::File::open(input).unwrap()).lines()
        .map(|l|serde_json::from_str(&l.unwrap()).unwrap()).collect();
    let scoped:Vec<_>=original.iter().enumerate().filter_map(|(i,r)|coordinates(r).map(|a|(i,r,a))).collect();
    assert!(!scoped.is_empty());
    crate::json_write(destination.join("scope.json"),&json!({"input_rows":original.len(),"selected_indices":scoped.iter().map(|(i,_,_)|i).collect::<Vec<_>>(),
        "rule":"Interior integer PMF counts <=1e6 trials, finite numeric answers, with at least one non-series deviance. Resource limit1e6 is research-only."}));
    let ps=profiles();assert_eq!(ps.len(),9000);
    crate::json_write(destination.join("profiles.json"),&ps);
    let mut kills=BufWriter::new(fs::File::create_new(destination.join("first-kills.jsonl")).unwrap());
    let mut survivors=Vec::new();let mut executions=0;let mut controls=0;
    // A strict baseline must reproduce every current product mass in scope.
    let baseline=ps.iter().find(|p|matches!(p.direct_model,EvalModel::Strict)&&p.log==0&&p.association==0&&p.store_ratio&&matches!(p.correction_model,EvalModel::Strict)&&p.store_correction&&p.normalizer_association==0&&p.normalizer_log==0).unwrap();
    for (_,row,(x,n,p)) in &scoped {
        let current=current::binom_mass(*x,*n,*p,1.0-p);
        let encoded=if current==0.0||current.is_subnormal(){"0x0000000000000000".into()}else{format!("0x{:016x}",current.to_bits())};
        assert_eq!(mass::<F>(*x,*n,*p,baseline),encoded,"baseline {}",row["id"]);controls+=1;
    }
    for profile in &ps {
        let mut alive=true;
        for (order,(index,row,(x,n,p))) in scoped.iter().enumerate() {
            let got=mass::<F>(*x,*n,*p,profile);executions+=1;
            if Some(got.as_str())!=row["expected"].as_str() {
                writeln!(kills,"{}",json!({"candidate_id":profile.id,"order":order,"trace_index":index,"id":row["id"],"got":got,"expected":row["expected"]})).unwrap();alive=false;break;
            }
        }
        if alive {survivors.push(profile.id);}
    }
    crate::json_write(destination.join("survivors.json"),&survivors);
    let result=json!({"complete":true,"profiles":ps.len(),"scoped":scoped.len(),"input_rows":original.len(),"executions":executions,"survivors":survivors.len(),"baseline_controls":controls,
        "limits":"Fixed arithmetic programs, no fitted constants. Current Stirling constants and near-mean deviance series are conditional leaves, not independently identified Excel producers. Endpoints/CDFs and other programs remain open."});
    crate::json_write(destination.join("summary.json"),&result);println!("{result}");
}

/// Test the defining incomplete-beta identities as whole fixed programs.
/// Their equivalence over the reals does not imply Excel uses any of them.
pub fn cdf(input:&str,output:&str) {
    let destination=Path::new(output);fs::create_dir(destination).unwrap();
    let original:Vec<Value>=BufReader::new(fs::File::open(input).unwrap()).lines()
        .map(|l|serde_json::from_str(&l.unwrap()).unwrap()).collect();
    let mut counts=std::collections::BTreeMap::<String,Vec<[usize;5]>>::new();
    let mut first=std::collections::BTreeMap::<String,Vec<Option<Value>>>::new();
    let mut selected=Vec::new();
    for (index,row) in original.iter().enumerate() {
        let Some(name)=row["function"].as_str()else{continue};
        if !matches!(name,"BINOM.DIST"|"BINOMDIST"|"NEGBINOM.DIST"|"NEGBINOMDIST") || !row["result_index"].is_null() {continue;}
        let Some(a)=row["args"].as_array()else{continue};
        if a.len()!=4 || !(a[3]==true || numeric(&a[3]).is_some_and(|n|n!=0.0)) {continue;}
        let Some(x)=numeric(&a[0])else{continue};
        let Some(n)=numeric(&a[1])else{continue};
        let Some(p)=numeric(&a[2])else{continue};
        let is_binom=name.starts_with("BINOM");
        if !(x>=0.0 && n>=1.0 && x==x.trunc() && n==n.trunc() && p>0.0 && p<1.0 && x<=10000.0 && n<=10000.0) || is_binom && x>=n {continue;}
        if !row["expected"].as_str().is_some_and(|s|s.starts_with("0x")) {continue;}
        selected.push(index);
        let q=1.0-p;
        let (left,right)=if is_binom {(n-x,x+1.0)}else{(n,x+1.0)};
        let (z,complement)=if is_binom {(q,p)}else{(p,q)};
        let direct=crate::ported::special_math_common::bratio(left,right,z,complement);
        let reverse=crate::ported::special_math_common::bratio(right,left,complement,z);
        let values=[direct.0,reverse.1,1.0-reverse.0,0.5+(0.5-reverse.0),
            crate::ported::special_math_common::regularized_beta(z,left,right)];
        let stats=counts.entry(name.into()).or_insert_with(||vec![[0;5];values.len()]);
        let killed=first.entry(name.into()).or_insert_with(||vec![None;values.len()]);
        let was=row["actual"]==row["expected"];
        for (program,value) in values.into_iter().enumerate() {
            let got=if !value.is_finite(){"error:Num".into()}else if value.is_subnormal()||value==0.0{"0x0000000000000000".into()}else{format!("0x{:016x}",value.to_bits())};
            let matched=Some(got.as_str())==row["expected"].as_str();
            stats[program][usize::from(!matched)]+=1;
            stats[program][2]+=usize::from(matched&&!was);
            stats[program][3]+=usize::from(!matched&&was);
            stats[program][4]+=usize::from(row["actual"]!=got);
            if !matched && killed[program].is_none(){killed[program]=Some(json!({"trace_index":index,"id":row["id"],"expected":row["expected"],"got":got}));}
        }
    }
    assert!(!selected.is_empty());
    crate::json_write(destination.join("scope.json"),&selected);
    crate::json_write(destination.join("first-kills.json"),&first);
    let result=json!({"complete":true,"scoped":selected.len(),"programs_per_function":5,"counts_pass_fail_gains_losses_changed":counts,
        "programs":["bratio(a,b,z,w).p","bratio(b,a,w,z).q","1-bratio(b,a,w,z).p","0.5+(0.5-bratio(b,a,w,z).p)","regularized_beta(z,a,b)"],
        "limits":"Existing incomplete-beta body; fixed defining identities, no fitted constants or answers. Integer nonendpoint counts/shapes <=10000 are a research resource scope. All rows are discovery. Exact rejection means any recorded mismatch; scores never authorize a row-specific program selector."});
    crate::json_write(destination.join("summary.json"),&result);println!("{result}");
}
