//! Finite Poisson mass programs from the defining PMF and its recurrence.
//! Tier 1 NIST e-Handbook 1.3.6.6.19 / IEEE 754, own algebra. No nmath,
//! discrete_dist source/history or other restricted implementation is read.
use crate::{dsl::*, eval::*, boundary::{neg,store}, json_write, ported::excel_numeric::*, value};
use serde::{Serialize,Deserialize};
use serde_json::{json,Value};
use std::{fs,io::{BufWriter,Write},path::Path,sync::atomic::{AtomicUsize,Ordering},time::Instant};

#[derive(Clone,Serialize,Deserialize)]
pub struct MassRow {pub trace_index:usize,pub raw:Value,pub count:usize,pub mean:u64}
impl MassRow {fn expected(&self)->&str {self.raw["expected"].as_str().unwrap()}}

pub fn prepare(trace:&str,output:&str) {
    fs::create_dir_all(output).unwrap();let out=Path::new(output);
    let mut scoped=BufWriter::new(fs::File::create(out.join("scoped.jsonl")).unwrap());
    let mut controls=BufWriter::new(fs::File::create(out.join("controls.jsonl")).unwrap());
    let mut counts=std::collections::BTreeMap::<String,usize>::new();
    for (trace_index,raw) in crate::lines::<Value>(trace).into_iter().enumerate() {
        assert_eq!(raw["function"],"POISSON.DIST");
        let decoded:Vec<_>=raw["args"].as_array().unwrap().iter().enumerate().map(|(i,a)|calc_verdict::decode_argument(a,&format!("args[{i}]"))).collect();
        let mut leaf=None;
        let reason=if decoded.len()!=3 || !raw["result_index"].is_null() {"unsupported_shape"}
        else if decoded.iter().any(|v|v.is_err()) {"failed_argument_decode"}
        else {
            let vals:Vec<_>=decoded.into_iter().map(Result::unwrap).collect();
            let cumulative=match vals[2].core() {value::CoreValue::Logical(b)=>Some(*b),value::CoreValue::Number(n)=>Some(*n!=0.0),_=>None};
            if cumulative==Some(true) {"cdf_control"}
            else if cumulative.is_none(){"argument_type_control"}
            else if let (value::CoreValue::Number(k),value::CoreValue::Number(mean))=(vals[0].core(),vals[1].core()) {
                if !k.is_finite()||!mean.is_finite()||k.is_subnormal()||mean.is_subnormal(){"input_admission_control"}
                else if *k<1.0||k.trunc()>170.0||*mean<=0.0 {"domain_or_other_count_control"}
                else {leaf=Some((k.trunc()as usize,mean.to_bits()));"scoped"}
            }else{"argument_type_control"}
        };
        *counts.entry(reason.into()).or_default()+=1;
        if let Some((count,mean))=leaf {
            calc_verdict::decode_result(raw["expected"].as_str().unwrap()).expect("unknown scoped outcome is failed load");
            serde_json::to_writer(&mut scoped,&MassRow{trace_index,raw,count,mean}).unwrap();writeln!(scoped).unwrap();
        }else{serde_json::to_writer(&mut controls,&json!({"trace_index":trace_index,"reason":reason,"raw":raw})).unwrap();writeln!(controls).unwrap();}
    }
    json_write(out.join("coverage.json"),&json!({"decoded_trace_counts":counts,"unread_source_rows":4,"unread_policy":"all four remain failures in application totals; no invented trace row","fresh_validation":0}));
    println!("Poisson prepare {}",json!(counts));
}

#[derive(Clone,Debug,Serialize,Deserialize)]
pub struct Program {
    pub id:usize,
    // 0 recurrence, 1 factorial with linear power, 2 factorial with binary power,
    // 3 exponential of k*ln(mean) - mean - sum(ln(i)).
    pub family:usize,pub association:usize,pub exp:usize,pub log:usize,
    pub models:Vec<EvalModel>,pub final_store:bool,
}
pub fn programs()->Vec<Program> {
    let mut out=vec![];
    for family in 0..4 {let sites=if family==0 {2}else{4};
        for association in 0..if family<3 {3}else{2} {
            for exp in 0..5 {for log in 0..if family==3 {4}else{1} {for mut code in 0..3usize.pow(sites) {
                let mut models=vec![];for _ in 0..sites{models.push(MODELS[code%3]);code/=3;}
                for final_store in [true,false]{out.push(Program{id:out.len(),family,association,exp,log,models:models.clone(),final_store});}
            }}}
        }
    }out
}
pub fn exp<F:Fpu>(x:Val,kind:usize)->Val {
    match kind {0|1=>store::<F>(Val::Ext(research_exp::<F>(&x.to_ext())),kind==0),
        2=>Val::f64(exp_portable(x.to_f64::<F>())),3=>Val::f64(libm::exp(x.to_f64::<F>())),4=>Val::f64(crt_exp(x.to_f64::<F>())),_=>unreachable!()}
}
pub fn log<F:Fpu>(x:f64,kind:usize)->Val {
    log_value::<F>(Val::f64(x),kind)
}
pub fn log_value<F:Fpu>(x:Val,kind:usize)->Val {
    match kind {0|1=>store::<F>(Val::Ext(research_log::<F>(&x.to_ext())),kind==0),
        2=>Val::f64(crt_log(x.to_f64::<F>())),3=>Val::f64(libm::log(x.to_f64::<F>())),_=>unreachable!()}
}
fn op<F:Fpu>(op:Op,model:EvalModel,a:Val,b:Val)->Val {binop::<F>(op,model,a,b)}
fn mul<F:Fpu>(model:EvalModel,a:Val,b:Val)->Val {op::<F>(Op::Mul(0,1),model,a,b)}
fn div<F:Fpu>(model:EvalModel,a:Val,b:Val)->Val {op::<F>(Op::Div(0,1),model,a,b)}
fn add<F:Fpu>(model:EvalModel,a:Val,b:Val)->Val {op::<F>(Op::Add(0,1),model,a,b)}
fn power<F:Fpu>(x:Val,mut k:usize,model:EvalModel,binary:bool)->Val {
    let mut out=Val::f64(1.0);
    if !binary {for _ in 0..k {out=mul::<F>(model,out,x)}return out}
    let mut base=x;
    while k!=0 {if k&1!=0{out=mul::<F>(model,out,base)}k>>=1;if k!=0{base=mul::<F>(model,base,base)}}out
}
pub struct Constants {factorials:Vec<[Val;3]>,log_factorials:Vec<[[Val;3];4]>}
fn constants<F:Fpu>()->Constants {
    let one=[Val::f64(1.0);3];let zero=[[Val::f64(0.0);3];4];
    let mut factorials=vec![one];let mut log_factorials=vec![zero];
    for i in 1..=170 {let mut f=one;let mut l=zero;
        for m in 0..3 {f[m]=mul::<F>(MODELS[m],factorials[i-1][m],Val::f64(i as f64));}
        for kind in 0..4 {let li=log::<F>(i as f64,kind);for m in 0..3 {l[kind][m]=add::<F>(MODELS[m],log_factorials[i-1][kind][m],li);}}
        factorials.push(f);log_factorials.push(l);
    }Constants{factorials,log_factorials}
}
struct Cache {exp:[Val;5],powers:[[Val;3];2],logs:[Val;4]}
fn cache<F:Fpu>(r:&MassRow)->Cache {
    let mean=Val::F64(r.mean);let x=f64::from_bits(r.mean);
    Cache{exp:std::array::from_fn(|i|exp::<F>(neg(mean),i)),powers:std::array::from_fn(|i|std::array::from_fn(|m|power::<F>(mean,r.count,MODELS[m],i==1))),logs:std::array::from_fn(|i|log::<F>(x,i))}
}
fn mi(m:EvalModel)->usize {match m{EvalModel::Strict=>0,EvalModel::StoredPc64=>1,EvalModel::ResidentPc64=>2}}
fn evaluate<F:Fpu>(p:&Program,r:&MassRow,c:&Constants,v:&Cache)->Val {
    let mean=Val::F64(r.mean);let k=Val::f64(r.count as f64);let ms=&p.models;
    let result=if p.family==0 {
        let mut term=v.exp[p.exp];for i in 1..=r.count {let i=Val::f64(i as f64);term=match p.association {
            0=>div::<F>(ms[1],mul::<F>(ms[0],term,mean),i),
            1=>mul::<F>(ms[1],term,div::<F>(ms[0],mean,i)),
            2=>mul::<F>(ms[1],div::<F>(ms[0],term,i),mean),_=>unreachable!(),
        }}term
    }else if p.family<3 {
        let power=v.powers[p.family-1][mi(ms[0])];let fact=c.factorials[r.count][mi(ms[1])];let exp=v.exp[p.exp];
        match p.association {
            0=>div::<F>(ms[3],mul::<F>(ms[2],exp,power),fact),
            1=>mul::<F>(ms[3],exp,div::<F>(ms[2],power,fact)),
            2=>mul::<F>(ms[3],div::<F>(ms[2],exp,fact),power),_=>unreachable!(),
        }
    }else{
        let lf=c.log_factorials[r.count][p.log][mi(ms[0])];let kl=mul::<F>(ms[1],k,v.logs[p.log]);
        let exponent=if p.association==0 {add::<F>(ms[3],add::<F>(ms[2],kl,neg(mean)),neg(lf))}
            else {add::<F>(ms[3],kl,neg(add::<F>(ms[2],mean,lf)))};
        exp::<F>(exponent,p.exp)
    };store::<F>(result,p.final_store)
}

pub fn run<F:Fpu+'static>(rows:&[MassRow],output:&str) {
    fs::create_dir_all(output).unwrap();let out=Path::new(output);let ps=programs();let c=constants::<F>();let start=Instant::now();
    json_write(out.join("programs.json"),&ps);let grammar=hash_bytes(&serde_json::to_vec(&ps).unwrap());
    let mut alive:Vec<_>=(0..ps.len()).collect();let mut evaluations=0;let mut summary_kills=vec![];
    let mut kills=BufWriter::new(fs::File::create(out.join("first-kills.jsonl")).unwrap());let mut stats=BufWriter::new(fs::File::create(out.join("witness-kills.jsonl")).unwrap());
    for (order,r) in rows.iter().enumerate(){let before=alive.len();let mut next=vec![];let v=cache::<F>(r);
        for id in alive {let got=outcome::<F>(evaluate::<F>(&ps[id],r,&c,&v));evaluations+=1;if got==r.expected(){next.push(id)}else{
            serde_json::to_writer(&mut kills,&json!({"candidate_id":id,"order":order,"trace_index":r.trace_index,"witness":r.raw,"expected":r.expected(),"got":got})).unwrap();writeln!(kills).unwrap();
        }}
        let stat=json!({"order":order,"trace_index":r.trace_index,"id":r.raw["id"],"killed":before-next.len(),"remaining":next.len()});
        if before!=next.len(){println!("mass {stat}");summary_kills.push(stat.clone())}
        serde_json::to_writer(&mut stats,&stat).unwrap();writeln!(stats).unwrap();alive=next;
    }kills.flush().unwrap();stats.flush().unwrap();
    let survivors:Vec<_>=alive.iter().map(|id|&ps[*id]).collect();json_write(out.join("survivors.json"),&survivors);
    let cursor=AtomicUsize::new(0);let progress=AtomicUsize::new(0);
    let workers=std::env::var("DNA_RESEARCH_WORKERS").ok().map(|v|v.parse::<usize>().unwrap()).unwrap_or(8);
    assert!(workers>=1 && workers<=8);
    let counts:Vec<_>=std::thread::scope(|scope|{let mut tasks=vec![];
        for worker in 0..workers {let cursor=&cursor;let progress=&progress;let c=&c;let ps=&ps;tasks.push(scope.spawn(move||{
            let mut f=BufWriter::new(fs::File::create(out.join(format!("reachability-{worker}.jsonl"))).unwrap());
            let mut u=BufWriter::new(fs::File::create(out.join(format!("unreachable-roots-{worker}.jsonl"))).unwrap());let mut counts=[0usize;4];
            loop {let order=cursor.fetch_add(1,Ordering::Relaxed);if order>=rows.len(){break}let r=&rows[order];let v=cache::<F>(r);let mut roots=vec![];
                for p in ps {let x=evaluate::<F>(p,r,c,&v);if !roots.contains(&x){roots.push(x)}}
                let outputs=publish_roots::<F>(&roots);let reachable=outputs.iter().any(|b|bits_outcome(*b)==r.expected());
                counts[0]+=1;counts[1]+=(!reachable)as usize;counts[2]+=outputs.len();counts[3]+=roots.len();
                serde_json::to_writer(&mut f,&json!({"order":order,"trace_index":r.trace_index,"id":r.raw["id"],"expected":r.expected(),"reachable":reachable,"outputs":outputs.iter().map(|b|bits_outcome(*b)).collect::<Vec<_>>()})).unwrap();writeln!(f).unwrap();
                if !reachable {serde_json::to_writer(&mut u,&json!({"order":order,"trace_index":r.trace_index,"expected":r.expected(),"raw_roots":roots.iter().map(|v|v.exact_encoding()).collect::<Vec<_>>()})).unwrap();writeln!(u).unwrap();}
                let done=progress.fetch_add(1,Ordering::Relaxed)+1;if done%2000==0{println!("mass reachability {done}/{} {:.1}s",rows.len(),start.elapsed().as_secs_f64());}
            }counts
        }));}tasks.into_iter().map(|h|h.join().unwrap()).collect()});
    let mut totals=[0usize;4];for c in counts{for i in 0..4{totals[i]+=c[i]}}assert_eq!(totals[0],rows.len());
    let summary=json!({"programs":ps.len(),"grammar_hash":grammar,"rows":rows.len(),"survivors":survivors.len(),"evaluations":evaluations,"kills":summary_kills,"whole_family_unreachable":totals[1],"whole_family_reachable":rows.len()-totals[1],"published_outcomes":totals[2],"tagged_roots":totals[3],"seconds":start.elapsed().as_secs_f64(),"complete":true,"backend":std::any::type_name::<F>(),"fresh_validation":0,"workers":workers});
    json_write(out.join("summary.json"),&summary);println!("{summary}");
}

pub fn self_test(rows:&[MassRow]) {
    let ps=programs();let h=constants::<Hard>();let s=constants::<Soft>();let mut comparisons=0;let mut raw_differences=0;
    for r in rows.iter().step_by((rows.len()/64).max(1)).take(64) {
        let hv=cache::<Hard>(r);let sv=cache::<Soft>(r);
        for p in &ps {let a=evaluate::<Hard>(p,r,&h,&hv);let b=evaluate::<Soft>(p,r,&s,&sv);comparisons+=1;raw_differences+=(a!=b)as usize;}
    }
    println!("{}",json!({"explicit_backend_comparisons":comparisons,"raw_differences":raw_differences,"note":"Transcendental substrate differences are retained; no assertion of full native/software equivalence"}));
}
