//! Count-one Poisson normalization identities, with independently derived constants.
use crate::{dsl::*,eval::*,boundary::{neg,store},poisson::{self,MassRow},json_write,ported::excel_numeric::*};
use serde::{Serialize,Deserialize};
use serde_json::{json,Value};
use std::{fs,io::{BufWriter,Write},path::Path,collections::HashMap,sync::atomic::{AtomicUsize,Ordering},time::Instant};

#[derive(Clone,Debug,PartialEq,Eq,Serialize,Deserialize)]
#[serde(rename_all="snake_case")]
pub enum Node { Mean, One, Number{bits:u64}, Constant{index:usize,extended:bool}, Neg(Box<Node>),
    ReciprocalLog{provider:usize,reciprocal_model:EvalModel}, DerivedG{tau_extended:bool,model:EvalModel},
    Add{site:usize,a:Box<Node>,b:Box<Node>}, Mul{site:usize,a:Box<Node>,b:Box<Node>},
    Div{site:usize,a:Box<Node>,b:Box<Node>}, Log{provider:usize},
    Exp{provider:usize,entry_store:bool,input:Box<Node>}, Power{provider:usize,exponent:i32},
    ScaledUnderflow{scaled:Box<Node>,denominator:Box<Node>,fallback:Box<Node>,site:usize,reciprocal_site:Option<usize>} }
#[derive(Clone,Debug,Serialize,Deserialize)]
pub struct Template {pub id:usize,pub deviance:usize,pub normalization:usize,pub log_provider:usize,
    pub exp_provider:usize,pub exp_entry_store:bool,pub s_extended:bool,pub other_extended:bool,pub graph:Node}
fn add(site:usize,a:Node,b:Node)->Node {Node::Add{site,a:Box::new(a),b:Box::new(b)}}
fn mul(site:usize,a:Node,b:Node)->Node {Node::Mul{site,a:Box::new(a),b:Box::new(b)}}
fn div(site:usize,a:Node,b:Node)->Node {Node::Div{site,a:Box::new(a),b:Box::new(b)}}
fn minus(a:Node)->Node {Node::Neg(Box::new(a))}
fn sub(site:usize,a:Node,b:Node)->Node {add(site,a,minus(b))}
fn count_sites(node:&Node,counts:&mut [usize;8]) {
    match node {Node::Add{site,a,b}|Node::Mul{site,a,b}|Node::Div{site,a,b}=>{counts[*site]+=1;count_sites(a,counts);count_sites(b,counts)},
        Node::Neg(a)|Node::Exp{input:a,..}=>count_sites(a,counts),
        Node::ScaledUnderflow{scaled,denominator,fallback,site,reciprocal_site}=>{
            counts[*site]+=1;if let Some(r)=reciprocal_site {counts[*r]+=1}
            count_sites(scaled,counts);count_sites(denominator,counts);count_sites(fallback,counts);
        },_=>{}}
}
pub fn templates()->Vec<Template> {
    let mut out=vec![];
    for deviance in 0..3 {for normalization in 0..5 {for log_provider in 0..4 {for exp_provider in 0..5 {
        for exp_entry_store in if exp_provider<2 {vec![true,false]}else{vec![false]} {
            for s_extended in [false,true] {for other_extended in [false,true] {
                let log=Node::Log{provider:log_provider};
                let d=match deviance {0=>sub(1,sub(0,Node::Mean,Node::One),log),
                    1=>sub(1,Node::Mean,add(0,Node::One,log)),2=>sub(1,sub(0,Node::Mean,log),Node::One),_=>unreachable!()};
                let s=Node::Constant{index:1,extended:s_extended};
                let other=Node::Constant{index:if normalization==1||normalization==4 {0}else{2},extended:other_extended};
                let exp=|input|Node::Exp{provider:exp_provider,entry_store:exp_entry_store,input:Box::new(input)};
                let graph=match normalization {
                    0=>div(3,exp(minus(add(2,s,d))),other),
                    1=>exp(sub(3,minus(add(2,s,d)),other)),
                    2=>div(3,mul(2,exp(minus(d)),exp(minus(s))),other),
                    3=>div(3,exp(minus(d)),mul(2,other,exp(s))),
                    4=>mul(3,exp(minus(d)),exp(minus(add(2,s,other)))),_=>unreachable!()};
                let mut sites=[0;8];count_sites(&graph,&mut sites);assert_eq!(sites,[1,1,1,1,0,0,0,0]);
                out.push(Template{id:out.len(),deviance,normalization,log_provider,exp_provider,exp_entry_store,s_extended,other_extended,graph});
            }}
        }
    }}}}
    assert_eq!(out.len(),1680);
    let mut hashes=HashMap::<String,Vec<&Node>>::new();
    for t in &out {let hash=hash_bytes(&serde_json::to_vec(&t.graph).unwrap());let bucket=hashes.entry(hash).or_default();assert!(!bucket.contains(&&t.graph));bucket.push(&t.graph);}
    out
}
pub fn expanded_templates()->Vec<Template> {
    let mut out=templates();
    for variant in 0..10 {for log_provider in 0..4 {for exp_provider in 0..5 {
        for entry in if exp_provider<2 {vec![true,false]}else{vec![false]} {
            let l=Node::Log{provider:log_provider};let two=||Node::Number{bits:2.0f64.to_bits()};
            let exp=|input|Node::Exp{provider:exp_provider,entry_store:entry,input:Box::new(input)};
            let scaled=||exp(sub(1,mul(0,two(),l.clone()),Node::Mean));
            let graph=match variant {
                0=>exp(sub(0,l,Node::Mean)),
                1=>mul(0,exp(l),exp(minus(Node::Mean))),
                2=>exp(sub(2,sub(1,mul(0,two(),l.clone()),Node::Mean),l)),
                3=>exp(sub(2,mul(0,two(),l.clone()),add(1,Node::Mean,l))),
                4=>exp(sub(2,sub(1,mul(0,two(),l.clone()),l),Node::Mean)),
                5=>div(2,scaled(),Node::Mean),
                6=>mul(3,scaled(),div(2,Node::One,Node::Mean)),
                7=>mul(2,div(1,exp(mul(0,two(),l)),Node::Mean),exp(minus(Node::Mean))),
                8=>Node::ScaledUnderflow{scaled:Box::new(scaled()),denominator:Box::new(Node::Mean),fallback:Box::new(exp(sub(3,l,Node::Mean))),site:2,reciprocal_site:None},
                9=>Node::ScaledUnderflow{scaled:Box::new(scaled()),denominator:Box::new(Node::Mean),fallback:Box::new(exp(sub(4,l,Node::Mean))),site:3,reciprocal_site:Some(2)},
                _=>unreachable!(),
            };
            out.push(Template{id:out.len(),deviance:0,normalization:5+variant,log_provider,exp_provider,exp_entry_store:entry,s_extended:false,other_extended:false,graph});
        }
    }}}
    for variant in 0..4 {for power_provider in 0..4 {for exp_provider in 0..5 {
        for entry in if exp_provider<2 {vec![true,false]}else{vec![false]} {
            let exp=Node::Exp{provider:exp_provider,entry_store:entry,input:Box::new(minus(Node::Mean))};
            let pow=Node::Power{provider:power_provider,exponent:match variant{0=>1,1=>-1,_=>2}};
            let graph=match variant {0=>mul(0,pow,exp),1=>div(0,exp,pow),
                2=>mul(1,exp,div(0,pow,Node::Mean)),3=>div(1,mul(0,exp,pow),Node::Mean),_=>unreachable!()};
            out.push(Template{id:out.len(),deviance:power_provider,normalization:15+variant,log_provider:0,exp_provider,exp_entry_store:entry,s_extended:false,other_extended:false,graph});
        }
    }}}
    let mut hashes=HashMap::<String,Vec<&Node>>::new();
    for t in &out {let hash=hash_bytes(&serde_json::to_vec(&t.graph).unwrap());let bucket=hashes.entry(hash).or_default();assert!(!bucket.contains(&&t.graph),"duplicate graph {}",t.id);bucket.push(&t.graph);variants(t);}
    assert_eq!(out.len(),2072);assert_eq!(out.iter().map(variants).sum::<usize>(),152040);out
}
fn models(mut code:usize)->[EvalModel;8] {std::array::from_fn(|_|{let m=MODELS[code%3];code/=3;m})}
fn model_index(m:EvalModel)->usize {match m{EvalModel::Strict=>0,EvalModel::StoredPc64=>1,EvalModel::ResidentPc64=>2}}
fn rewrite(node:&Node,mapper:&impl Fn(&Node)->Option<Node>)->Node {
    if let Some(n)=mapper(node){return n}
    match node {
        Node::Neg(a)=>minus(rewrite(a,mapper)),
        Node::Add{site,a,b}=>add(*site,rewrite(a,mapper),rewrite(b,mapper)),
        Node::Mul{site,a,b}=>mul(*site,rewrite(a,mapper),rewrite(b,mapper)),
        Node::Div{site,a,b}=>div(*site,rewrite(a,mapper),rewrite(b,mapper)),
        Node::Exp{provider,entry_store,input}=>Node::Exp{provider:*provider,entry_store:*entry_store,input:Box::new(rewrite(input,mapper))},
        Node::ScaledUnderflow{scaled,denominator,fallback,site,reciprocal_site}=>Node::ScaledUnderflow{
            scaled:Box::new(rewrite(scaled,mapper)),denominator:Box::new(rewrite(denominator,mapper)),
            fallback:Box::new(rewrite(fallback,mapper)),site:*site,reciprocal_site:*reciprocal_site},
        _=>node.clone()
    }
}
fn has(node:&Node,p:&impl Fn(&Node)->bool)->bool {
    if p(node){return true}
    match node {Node::Neg(a)|Node::Exp{input:a,..}=>has(a,p),
        Node::Add{a,b,..}|Node::Mul{a,b,..}|Node::Div{a,b,..}=>has(a,p)||has(b,p),
        Node::ScaledUnderflow{scaled,denominator,fallback,..}=>has(scaled,p)||has(denominator,p)||has(fallback,p),_=>false}
}
pub fn refined_templates()->Vec<Template> {
    let base=expanded_templates();let mut out=base.clone();
    for t in base.iter().filter(|t|has(&t.graph,&|n|matches!(n,Node::Log{..}))) {
        for reciprocal_model in MODELS {let mut c=t.clone();c.id=out.len();
            c.graph=rewrite(&t.graph,&|n|if let Node::Log{provider}=n {Some(minus(Node::ReciprocalLog{provider:*provider,reciprocal_model}))}else{None});out.push(c);
        }
    }
    assert_eq!(out.len(),7952);assert_eq!(out.iter().map(variants).sum::<usize>(),606144);
    let mut hashes=HashMap::<String,Vec<usize>>::new();
    for (i,t) in out.iter().enumerate(){hashes.entry(hash_bytes(&serde_json::to_vec(&t.graph).unwrap())).or_default().push(i);}
    let original=out.clone();
    for t in original.iter().filter(|t|has(&t.graph,&|n|matches!(n,Node::Constant{index:2,..}))) {
        for tau_extended in [false,true] {for model in MODELS {
            let mut c=t.clone();c.graph=rewrite(&t.graph,&|n|if matches!(n,Node::Constant{index:2,..}){Some(Node::DerivedG{tau_extended,model})}else{None});
            let hash=hash_bytes(&serde_json::to_vec(&c.graph).unwrap());let bucket=hashes.entry(hash).or_default();
            if !bucket.iter().any(|i|out[*i].graph==c.graph){c.id=out.len();bucket.push(c.id);out.push(c);}
        }}
    }
    assert_eq!(out.len(),20048);assert_eq!(out.iter().map(variants).sum::<usize>(),1585920);
    for t in &out {variants(t);}out
}
fn variants(t:&Template)->usize {
    let mut sites=[0;8];count_sites(&t.graph,&mut sites);let arity=sites.iter().filter(|c|**c!=0).count();
    assert!(sites[..arity].iter().all(|c|*c==1));assert!(sites[arity..].iter().all(|c|*c==0));3usize.pow(arity as u32)
}
fn is_zero(v:Val)->bool {let e=v.to_ext();e.0[..8].iter().all(|b|*b==0)&&(e.0[8]==0)&&(e.0[9]&0x7f==0)}
fn power<F:Fpu>(mean:u64,provider:usize,exponent:i32)->Val {
    let x=f64::from_bits(mean);let e=exponent as f64;
    match provider {0|1=>store::<F>(Val::Ext(research_pow::<F>(x,e)),provider==1),
        2=>Val::f64(crt_pow(x,e)),3=>Val::f64(libm::pow(x,e)),_=>unreachable!()}
}
type Constants=[[Val;2];7];
pub fn constants<F:Fpu>(path:&str)->Constants {
    let data:Value=serde_json::from_slice(&fs::read(path).unwrap()).unwrap();
    let mut result=[[Val::f64(0.0);2];7];
    for (idx,name) in ["C","S","G","TAU"].iter().enumerate(){
        let c=&data["constants"][*name];if c.is_null(){assert_eq!(idx,3);continue}
        let b=u64::from_str_radix(c["binary64"].as_str().unwrap(),16).unwrap();
        let raw=c["ext80_le"].as_str().unwrap();let mut ext=[0u8;10];for i in 0..10{ext[i]=u8::from_str_radix(&raw[2*i..2*i+2],16).unwrap()}
        result[idx]=[Val::F64(b),Val::Ext(Extended(ext))];
    }
    if !data["constants"]["TAU"].is_null(){for t in 0..2{for m in 0..3 {
        let tau=result[3][t];let v=if m==0 {Val::f64(tau.to_f64::<F>().sqrt())}
            else {store::<F>(Val::Ext(F::sqrt(&tau.to_ext(),CW)),m==1)};
        result[4+m][t]=v;
    }}}
    result
}
fn evaluate<F:Fpu>(node:&Node,mean:u64,constants:&Constants,models:&[EvalModel;8])->Val {
    match node {Node::Mean=>Val::F64(mean),Node::One=>Val::f64(1.0),Node::Number{bits}=>Val::F64(*bits),Node::Constant{index,extended}=>constants[*index][*extended as usize],
        Node::Power{provider,exponent}=>power::<F>(mean,*provider,*exponent),
        Node::ReciprocalLog{provider,reciprocal_model}=>poisson::log_value::<F>(binop::<F>(Op::Div(0,1),*reciprocal_model,Val::f64(1.0),Val::F64(mean)),*provider),
        Node::DerivedG{tau_extended,model}=>{assert!(!is_zero(constants[3][0]),"derived sqrt(2*pi) requires a TAU constant certificate");constants[4+model_index(*model)][*tau_extended as usize]},
        Node::ScaledUnderflow{scaled,denominator,fallback,site,reciprocal_site}=>{
            let s=evaluate::<F>(scaled,mean,constants,models);
            if is_zero(s) {evaluate::<F>(fallback,mean,constants,models)} else {
                let d=evaluate::<F>(denominator,mean,constants,models);
                if let Some(r)=reciprocal_site {let inv=binop::<F>(Op::Div(0,1),models[*r],Val::f64(1.0),d);binop::<F>(Op::Mul(0,1),models[*site],s,inv)}
                else {binop::<F>(Op::Div(0,1),models[*site],s,d)}
            }
        },
        Node::Neg(a)=>neg(evaluate::<F>(a,mean,constants,models)),Node::Log{provider}=>poisson::log::<F>(f64::from_bits(mean),*provider),
        Node::Exp{provider,entry_store,input}=>poisson::exp::<F>(store::<F>(evaluate::<F>(input,mean,constants,models),*entry_store),*provider),
        Node::Add{site,a,b}|Node::Mul{site,a,b}|Node::Div{site,a,b}=>{
            let op=match node {Node::Add{..}=>Op::Add(0,1),Node::Mul{..}=>Op::Mul(0,1),Node::Div{..}=>Op::Div(0,1),_=>unreachable!()};
            binop::<F>(op,models[*site],evaluate::<F>(a,mean,constants,models),evaluate::<F>(b,mean,constants,models))
        }
    }
}
fn roots<F:Fpu>(node:&Node,mean:u64,constants:&Constants)->Vec<Val> {
    match node {Node::Mean=>vec![Val::F64(mean)],Node::One=>vec![Val::f64(1.0)],Node::Number{bits}=>vec![Val::F64(*bits)],Node::Constant{index,extended}=>vec![constants[*index][*extended as usize]],
        Node::Power{provider,exponent}=>vec![power::<F>(mean,*provider,*exponent)],
        Node::ReciprocalLog{provider,reciprocal_model}=>vec![poisson::log_value::<F>(binop::<F>(Op::Div(0,1),*reciprocal_model,Val::f64(1.0),Val::F64(mean)),*provider)],
        Node::DerivedG{tau_extended,model}=>{assert!(!is_zero(constants[3][0]),"derived sqrt(2*pi) requires a TAU constant certificate");vec![constants[4+model_index(*model)][*tau_extended as usize]]},
        Node::ScaledUnderflow{scaled,denominator,fallback,reciprocal_site,..}=>{
            let mut out=vec![];
            for s in roots::<F>(scaled,mean,constants) {
                if is_zero(s) {for f in roots::<F>(fallback,mean,constants) {if !out.contains(&f){out.push(f)}}}
                else {for d in roots::<F>(denominator,mean,constants) {for m in MODELS {
                    if reciprocal_site.is_some(){for q in MODELS {let inv=binop::<F>(Op::Div(0,1),q,Val::f64(1.0),d);let v=binop::<F>(Op::Mul(0,1),m,s,inv);if !out.contains(&v){out.push(v)}}}
                    else {let v=binop::<F>(Op::Div(0,1),m,s,d);if !out.contains(&v){out.push(v)}}
                }}}
            }out
        },
        Node::Neg(a)=>roots::<F>(a,mean,constants).into_iter().map(neg).collect(),Node::Log{provider}=>vec![poisson::log::<F>(f64::from_bits(mean),*provider)],
        Node::Exp{provider,entry_store,input}=>{let mut out=vec![];for v in roots::<F>(input,mean,constants){let x=poisson::exp::<F>(store::<F>(v,*entry_store),*provider);if !out.contains(&x){out.push(x)}}out},
        Node::Add{a,b,..}|Node::Mul{a,b,..}|Node::Div{a,b,..}=>{
            let op=match node {Node::Add{..}=>Op::Add(0,1),Node::Mul{..}=>Op::Mul(0,1),Node::Div{..}=>Op::Div(0,1),_=>unreachable!()};
            let av=roots::<F>(a,mean,constants);let bv=roots::<F>(b,mean,constants);let mut out=vec![];
            for a in av {for b in &bv {for model in MODELS {let x=binop::<F>(op,model,a,*b);if !out.contains(&x){out.push(x)}}}}out
        }
    }
}
pub fn run<F:Fpu+'static>(all:&[MassRow],const_path:&str,output:&str) {run_inner::<F>(all,const_path,output,0)}
pub fn run_expanded<F:Fpu+'static>(all:&[MassRow],const_path:&str,output:&str) {run_inner::<F>(all,const_path,output,1)}
pub fn run_refined<F:Fpu+'static>(all:&[MassRow],const_path:&str,output:&str) {run_inner::<F>(all,const_path,output,2)}
fn run_inner<F:Fpu+'static>(all:&[MassRow],const_path:&str,output:&str,level:usize) {
    fs::create_dir_all(output).unwrap();let out=Path::new(output);let ts=match level{0=>templates(),1=>expanded_templates(),_=>refined_templates()};let constants=constants::<F>(const_path);let start=Instant::now();
    let rows:Vec<_>=all.iter().filter(|r|r.count==1).collect();
    let mut programs=vec![];let mut layout=vec![];
    for t in &ts {let count=variants(t);layout.push(json!({"template":t.id,"first_candidate_id":programs.len(),"variants":count}));for code in 0..count {programs.push((t.id,code));}}
    let total=programs.len();let grammar=hash_bytes(&serde_json::to_vec(&ts).unwrap());json_write(out.join("program-layout.json"),&layout);
    json_write(out.join("templates.json"),&ts);let mut alive:Vec<_>=(0..total).collect();let mut evaluations=0;let mut kill_summary=vec![];
    let mut kills=BufWriter::new(fs::File::create(out.join("first-kills.jsonl")).unwrap());let mut stats=BufWriter::new(fs::File::create(out.join("witness-kills.jsonl")).unwrap());
    for (order,r) in rows.iter().enumerate(){let before=alive.len();let mut next=vec![];
        for id in alive {let (template_id,model_code)=programs[id];let got=outcome::<F>(evaluate::<F>(&ts[template_id].graph,r.mean,&constants,&models(model_code)));evaluations+=1;
            if got==r.raw["expected"].as_str().unwrap(){next.push(id)}else{
                serde_json::to_writer(&mut kills,&json!({"candidate_id":id,"template_id":template_id,"model_code":model_code,"order":order,"trace_index":r.trace_index,"expected":r.raw["expected"],"got":got})).unwrap();writeln!(kills).unwrap();
            }}
        let record=json!({"order":order,"trace_index":r.trace_index,"id":r.raw["id"],"killed":before-next.len(),"remaining":next.len()});
        if before!=next.len(){println!("count-one {record}");kill_summary.push(record.clone());}
        serde_json::to_writer(&mut stats,&record).unwrap();writeln!(stats).unwrap();alive=next;
    }kills.flush().unwrap();stats.flush().unwrap();json_write(out.join("survivors.json"),&alive);
    let cursor=AtomicUsize::new(0);let progress=AtomicUsize::new(0);let workers=std::env::var("DNA_RESEARCH_WORKERS").ok().map(|v|v.parse::<usize>().unwrap()).unwrap_or(4);assert!(workers>=1&&workers<=8);
    let counts:Vec<_>=std::thread::scope(|scope|{let mut tasks=vec![];
        for worker in 0..workers{let cursor=&cursor;let progress=&progress;let rows=&rows;let ts=&ts;let constants=&constants;tasks.push(scope.spawn(move||{
            let mut f=BufWriter::new(fs::File::create(out.join(format!("reachability-{worker}.jsonl"))).unwrap());let mut u=BufWriter::new(fs::File::create(out.join(format!("unreachable-roots-{worker}.jsonl"))).unwrap());let mut counts=[0usize;4];
            loop {let order=cursor.fetch_add(1,Ordering::Relaxed);if order>=rows.len(){break}let r=rows[order];let mut values=vec![];
                for t in ts {for v in roots::<F>(&t.graph,r.mean,constants) {if !values.contains(&v){values.push(v)}}}
                let outputs=publish_roots::<F>(&values);let reachable=outputs.iter().any(|b|bits_outcome(*b)==r.raw["expected"].as_str().unwrap());counts[0]+=1;counts[1]+=(!reachable)as usize;counts[2]+=outputs.len();counts[3]+=values.len();
                serde_json::to_writer(&mut f,&json!({"order":order,"trace_index":r.trace_index,"id":r.raw["id"],"expected":r.raw["expected"],"reachable":reachable,"outputs":outputs.iter().map(|b|bits_outcome(*b)).collect::<Vec<_>>()})).unwrap();writeln!(f).unwrap();
                if !reachable {serde_json::to_writer(&mut u,&json!({"order":order,"trace_index":r.trace_index,"expected":r.raw["expected"],"raw_roots":values.iter().map(|v|v.exact_encoding()).collect::<Vec<_>>()})).unwrap();writeln!(u).unwrap();}
                let done=progress.fetch_add(1,Ordering::Relaxed)+1;if done%500==0 {println!("count-one reachability {done}/{} {:.1}s",rows.len(),start.elapsed().as_secs_f64());}
            }counts
        }));}tasks.into_iter().map(|h|h.join().unwrap()).collect()});
    let mut totals=[0usize;4];for c in counts {for i in 0..4 {totals[i]+=c[i]}}assert_eq!(totals[0],rows.len());
    let summary=json!({"grammar_hash":grammar,"programs":total,"templates":ts.len(),"models_per_template":if level>0{json!(null)}else{json!(81)},"rows":rows.len(),"survivors":alive.len(),"evaluations":evaluations,"kills":kill_summary,"whole_family_unreachable":totals[1],"whole_family_reachable":rows.len()-totals[1],"published_outcomes":totals[2],"tagged_roots":totals[3],"complete":true,"seconds":start.elapsed().as_secs_f64(),"backend":std::any::type_name::<F>(),"workers":workers,"fresh_validation":0});
    json_write(out.join("summary.json"),&summary);println!("{summary}");
}
pub fn self_test(all:&[MassRow],const_path:&str) {test_inner(all,const_path,0)}
pub fn test_expanded(all:&[MassRow],const_path:&str) {test_inner(all,const_path,1)}
pub fn test_refined(all:&[MassRow],const_path:&str) {test_inner(all,const_path,2)}
fn test_inner(all:&[MassRow],const_path:&str,level:usize) {
    let ts=match level{0=>templates(),1=>expanded_templates(),_=>refined_templates()};let constants=self::constants::<Hard>(const_path);let software_constants=self::constants::<Soft>(const_path);let all_rows:Vec<_>=all.iter().filter(|r|r.count==1).collect();
    let rows:Vec<_>=[0,1,5,22,56,89,90,178,500,1000,2000,4231].into_iter().map(|i|all_rows[i]).collect();let mut tested=0;let mut backend_differences=0;
    for r in rows {for t in &ts {let expected=roots::<Hard>(&t.graph,r.mean,&constants);let mut explicit=vec![];
        let soft=roots::<Soft>(&t.graph,r.mean,&software_constants);backend_differences+=(expected.len()!=soft.len()||expected.iter().any(|v|!soft.contains(v)))as usize;
        for code in 0..variants(t) {let v=evaluate::<Hard>(&t.graph,r.mean,&constants,&models(code));if !explicit.contains(&v){explicit.push(v)}tested+=1;}
        assert_eq!(expected.len(),explicit.len());for v in expected {assert!(explicit.contains(&v));}
    }}println!("{}",json!({"explicit_template_model_checks":tested,"backend_template_set_differences":backend_differences,"note":"Synthetic arithmetic enumeration audit, not fresh Excel validation"}));
}
