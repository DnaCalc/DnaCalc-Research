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

/// Fixed hypergeometric quotient programs from the conditional-binomial identity.
/// Existing binomial masses are explicit conditional leaves; no table is changed.
pub fn hypergeom_quotients<F:Fpu>(input:&str,output:&str) {
    use calc_func::value::CalcValue;
    let destination=Path::new(output);fs::create_dir(destination).unwrap();
    let rows=crate::lines::<Value>(input);
    let mut profiles=Vec::new();
    for model in crate::dsl::WIDE_MODELS {
        for association in 0..3 {
            for complement in 0..2 {
                profiles.push((model,association,complement));
            }
        }
    }
    crate::json_write(destination.join("profiles.json"),&profiles);
    let mut counts=vec![[0_usize;5];profiles.len()];
    let mut first=vec![None;profiles.len()];
    let mut excluded=std::collections::BTreeMap::<String,usize>::new();
    let mut file=BufWriter::new(fs::File::create_new(destination.join("predictions.jsonl")).unwrap());
    for (index,row) in rows.iter().enumerate() {
        let name=row["function"].as_str().unwrap();
        if !matches!(name,"HYPGEOM.DIST"|"HYPGEOMDIST") || !row["result_index"].is_null() {
            *excluded.entry("surface_or_selection".into()).or_default()+=1;continue;
        }
        let Some(args)=row["args"].as_array() else {*excluded.entry("shape".into()).or_default()+=1;continue;};
        let values:Option<Vec<_>>=args.iter().enumerate().map(|(i,v)| {
            if i==4 {v.as_bool().map(|b|if b {1.0}else{0.0}).or_else(||numeric(v))} else {numeric(v)}
        }).collect();
        let Some(a)=values else {*excluded.entry("non_numeric_argument".into()).or_default()+=1;continue};
        if a.len()!=if name=="HYPGEOMDIST" {4}else{5} {*excluded.entry("shape".into()).or_default()+=1;continue;}
        let (s,n,k,total)=(a[0],a[1],a[2],a[3]);
        if !a[..4].iter().all(|v|v.is_finite() && *v>=0.0 && *v==v.trunc())
            || n>total || k>total || total>10000.0 || s>10000.0 {
            *excluded.entry("outside_integer_domain_or_resource_bound".into()).or_default()+=1;continue;
        }
        let cumulative=a.get(4).is_some_and(|v|*v!=0.0);
        let original_args=args.iter().enumerate().map(|(i,v)|calc_verdict::decode_argument(v,&format!("args[{i}]")).unwrap()).collect();
        let (baseline,_)=crate::invoke_control(name,original_args);
        assert_eq!(baseline,row["actual"].as_str().unwrap(),"frozen baseline {}",row["id"]);
        let expected=row["expected"].as_str().unwrap();
        for (id,&(model,association,complement)) in profiles.iter().enumerate() {
            let mass=|x:f64| {
                if x>n || x>k || n-x>total-k {return Val::f64(0.0);}
                if n==0.0 || n==total {return Val::f64(if x==if n==0.0 {0.0}else{k} {1.0}else{0.0});}
                let p=n/total;
                let q=if complement==0 {(total-n)/total}else{1.0-p};
                let left=Val::f64(current::binom_mass(x,k,p,q));
                let right=Val::f64(current::binom_mass(n-x,total-k,p,q));
                let denominator=Val::f64(current::binom_mass(n,total,p,q));
                match association {
                    0=>div::<F>(model,mul::<F>(model,left,right),denominator),
                    1=>mul::<F>(model,left,div::<F>(model,right,denominator)),
                    2=>mul::<F>(model,div::<F>(model,left,denominator),right),
                    _=>unreachable!(),
                }
            };
            let value=if cumulative && s>=n.min(k) {1.0}
                else if cumulative {
                    let mut sum=Val::f64(0.0);
                    for x in 0..=s as usize {sum=add::<F>(model,sum,mass(x as f64));}
                    sum.to_f64::<F>().min(1.0)
                } else {mass(s).to_f64::<F>()};
            let prediction=if !value.is_finite() {"error:Num".into()}
                else {calc_verdict::encode_result(&CalcValue::number(if value.is_subnormal(){0.0}else{value})).unwrap()};
            if id==0 {assert_eq!(prediction,baseline,"strict quotient control {}",row["id"]);}
            let was=baseline==expected;let matched=prediction==expected;
            for (slot,increment) in counts[id].iter_mut().zip([1,usize::from(matched),usize::from(was),usize::from(matched&&!was),usize::from(was&&!matched)]) {*slot+=increment;}
            if !matched && first[id].is_none() {first[id]=Some(json!({"trace_index":index,"id":row["id"],"prediction":prediction,"expected":expected}));}
            writeln!(file,"{}",json!({"trace_index":index,"id":row["id"],"program":id,"prediction":prediction,"expected":expected,"baseline":baseline})).unwrap();
        }
    }
    assert!(counts[0][0]>0,"empty hypergeometric quotient scope");
    assert_eq!(counts[0][0]+excluded.values().sum::<usize>(),rows.len(),"unaccounted input rows");
    let summary=json!({"complete":true,"input_rows":rows.len(),"profiles":profiles,"counts_all_exact_baseline_gains_losses":counts,"first_contradictions":first,"excluded":excluded,
        "source":"Hypergeometric conditional-binomial identity; IEEE754/Intel x87 coherent arithmetic models; current cleared binomial mass leaves.",
        "limits":"30 fixed quotient/complement programs with coherent CDF sums; integer domain bounded at10000 for research resources. Conditional binomial leaves, discovery answers, no fitted parameter or general kernel identification."});
    crate::json_write(destination.join("summary.json"),&summary);
    println!("{summary}");
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

/// Q(2,mean)=exp(-mean)*(1+mean), with two fixed algebraic forms.
/// EXP is the existing cleared worksheet helper, not an inferred private leaf.
pub fn poisson_one_cdf<F:Fpu>(input:&str,output:&str) {
    let destination=Path::new(output);fs::create_dir(destination).unwrap();
    let rows=crate::lines::<Value>(input);
    let profiles:Vec<_>=crate::dsl::WIDE_MODELS.into_iter().flat_map(|model|[(model,0),(model,1)]).collect();
    let mut counts=vec![[0_usize;5];profiles.len()];
    let mut first=vec![None;profiles.len()];
    let mut excluded=std::collections::BTreeMap::<String,usize>::new();
    let mut file=BufWriter::new(fs::File::create_new(destination.join("predictions.jsonl")).unwrap());
    for (index,row) in rows.iter().enumerate() {
        let name=row["function"].as_str().unwrap();
        if !matches!(name,"POISSON.DIST"|"POISSON") || !row["result_index"].is_null() {*excluded.entry("surface_or_selection".into()).or_default()+=1;continue;}
        let Some(a)=row["args"].as_array()else{*excluded.entry("shape".into()).or_default()+=1;continue;};
        if a.len()!=3 {*excluded.entry("shape".into()).or_default()+=1;continue;}
        let (Some(x),Some(mean))=(numeric(&a[0]),numeric(&a[1]))else{*excluded.entry("non_numeric".into()).or_default()+=1;continue;};
        let cumulative=a[2]==true || numeric(&a[2]).is_some_and(|v|v!=0.0);
        if !(x>=1.0 && x<2.0 && mean>=0.0 && mean.is_finite() && cumulative) {*excluded.entry("outside_one_count_positive_cdf".into()).or_default()+=1;continue;}
        let original_args=a.iter().enumerate().map(|(i,v)|calc_verdict::decode_argument(v,&format!("args[{i}]")).unwrap()).collect();
        let (baseline,_)=crate::invoke_control(name,original_args);
        assert_eq!(baseline,row["actual"].as_str().unwrap(),"frozen Poisson baseline {}",row["id"]);
        let expected=row["expected"].as_str().unwrap();
        let exponential=Val::f64(excel_exp(-mean));
        for (id,&(model,form)) in profiles.iter().enumerate() {
            let v=if form==0 {mul::<F>(model,exponential,add::<F>(model,Val::f64(1.0),Val::f64(mean)))}
                else {add::<F>(model,exponential,mul::<F>(model,Val::f64(mean),exponential))};
            let value=v.to_f64::<F>();
            let prediction=if !value.is_finite(){"error:Num".into()}
                else{calc_verdict::encode_result(&calc_func::value::CalcValue::number(if value.is_subnormal(){0.0}else{value})).unwrap()};
            let was=baseline==expected;let matched=prediction==expected;
            for (slot,increment) in counts[id].iter_mut().zip([1,usize::from(matched),usize::from(was),usize::from(matched&&!was),usize::from(was&&!matched)]) {*slot+=increment;}
            if !matched && first[id].is_none(){first[id]=Some(json!({"trace_index":index,"id":row["id"],"prediction":prediction,"expected":expected}));}
            writeln!(file,"{}",json!({"trace_index":index,"id":row["id"],"program":id,"prediction":prediction,"baseline":baseline,"expected":expected})).unwrap();
        }
    }
    assert!(counts[0][0]>0,"empty Poisson CDF scope");
    assert_eq!(counts[0][0]+excluded.values().sum::<usize>(),rows.len());
    let summary=json!({"complete":true,"input_rows":rows.len(),"profiles":profiles,"counts_all_exact_baseline_gains_losses":counts,"first_contradictions":first,"excluded":excluded,
        "source":"Defining Poisson finite sum at count1; DLMF8.4.8; IEEE754 and Intel x87 coherent arithmetic models; existing cleared EXP helper.",
        "limits":"Ten fixed expressions, positive finite means, original count in[1,2), cumulative only. EXP leaf is conditional. Discovery, no parameter fitting, row selection or complete gamma/Poisson algorithm identification."});
    crate::json_write(destination.join("summary.json"),&summary);println!("{summary}");
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
