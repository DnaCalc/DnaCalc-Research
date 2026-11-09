//! Fixed Loader binomial arithmetic programs, using the current cleared
//! Stirling table and near-mean series. No fitted constants or row selectors.
use crate::{dsl::{EvalModel, Op}, eval::{Val, binop, CW}, ported::{discrete_dist as current, excel_numeric::*}};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::{fs, io::{BufRead, BufReader, BufWriter, Write}, path::Path};

/// Fixed tiny-x beta prefactor trees. The original TOMS708 log-beta is an
/// explicit conditional leaf; public GAMMALN is not substituted for it.
pub fn beta_prefactors(input: &str, output: &str) {
    let rows = crate::lines::<Value>(input);
    let mut records = Vec::new();
    for row in rows {
        let args: Vec<_> = row["args"].as_array().unwrap().iter().map(|v|
            f64::from_bits(u64::from_str_radix(v.as_str().unwrap().trim_start_matches("0x"),16).unwrap())).collect();
        assert_eq!(args.len(), 3);
        let (x,a,b) = (args[0],args[1],args[2]);
        assert!(x > 0.0 && x < 1.0 && a >= 1.0 && b >= 1.0);
        let current_normalizer = crate::ported::special_math_common::bratio_betaln(a,b);
        let normalizer = row["normalizer"].as_str().map_or(current_normalizer, |raw|
            f64::from_bits(u64::from_str_radix(raw,16).unwrap()));
        assert!(normalizer.is_finite());
        let start = excel_x87_mul(a,excel_log(x));
        let difference = Hard::sub(&research_load(start), &research_load(normalizer), 0x133f);
        let resident = Hard::store(&research_exp::<Hard>(&difference), 0x133f);
        let variants = [
            ("stored", excel_exp(excel_x87_sub(start,normalizer))),
            ("compensated", excel_exp_compensated(start,&[normalizer])),
            ("running", excel_exp_running(start,&[normalizer])),
            ("resident", resident),
        ];
        let predictions: serde_json::Map<String,Value> = variants.into_iter().flat_map(|(name,value)| [
            (name.into(),json!(format!("{:016x}",excel_x87_div(value,a).to_bits()))),
            (format!("{name}_reciprocal"),json!(format!("{:016x}",excel_x87_mul(value,excel_x87_div(1.0,a)).to_bits()))),
        ]).collect();
        records.push(json!({"original":row,"normalizer":format!("{:016x}",normalizer.to_bits()),
            "current_normalizer":format!("{:016x}",current_normalizer.to_bits()),
            "start":format!("{:016x}",start.to_bits()),"predictions":predictions}));
    }
    crate::json_write(output,&json!({"complete":true,"records":records,
        "limits":"Conditional original TOMS708 beta normalizer. Fixed stored, compensated, running and resident exponent programs; no coefficient fitting, no application replacement or exactness claim."}));
}

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
    #[serde(default, skip_serializing_if="Option::is_none")]
    series: Option<SeriesProgram>,
}
#[derive(Clone,Copy,Debug,Serialize,Deserialize)]
struct SeriesProgram {model:EvalModel,local_store:bool}

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
                                        normalizer_association, normalizer_log, series:None});
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
        if let Some(program)=p.series{return Val::f64(series_deviance::<F>(x,mean,program));}
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

fn series_deviance<F:Fpu>(x:f64,mean:f64,p:SeriesProgram)->f64 {
    let m=p.model;let local=|v:Val|if p.local_store{Val::f64(v.to_f64::<F>())}else{v};
    let a=Val::f64(x);let b=Val::f64(mean);let d=local(sub::<F>(m,a,b));
    let mut v=local(div::<F>(m,d,add::<F>(m,a,b)));
    let mut sum=local(mul::<F>(m,d,v));
    let mut term=local(mul::<F>(m,mul::<F>(m,Val::f64(2.0),a),v));
    v=local(mul::<F>(m,v,v));
    for j in 1..10000 {
        term=local(mul::<F>(m,term,v));
        let denominator=add::<F>(m,mul::<F>(m,Val::f64(2.0),Val::f64(j as f64)),Val::f64(1.0));
        let next=local(add::<F>(m,sum,div::<F>(m,term,denominator)));
        if next.exact_encoding()==sum.exact_encoding(){return next.to_f64::<F>();}
        sum=next;
    }
    panic!("deviance series failed its geometric convergence bound"); // ubs:ignore — standalone research receipt must fail closed on a nonconverging candidate.
}

pub fn binomial_series<F:Fpu>(input:&str,output:&str) {
    let mut programs=vec![None];for model in crate::dsl::WIDE_MODELS{for local_store in [false,true]{programs.push(Some(SeriesProgram{model,local_store}));}}
    binomial_study::<F,_>(input,output,&programs,false,|x,n,p,series|{
        let profile=Profile{id:0,direct_model:EvalModel::Strict,log:0,association:0,store_ratio:true,
            correction_model:EvalModel::Strict,store_correction:true,normalizer_association:0,normalizer_log:0,series:*series};
        mass::<F>(x,n,p,&profile)
    },"Loader published deviance near-mean power series; existing cleared Stirling, far-deviance, CRT log and exponential program unchanged.",
      "11 fixed policies: original control plus5 coherent arithmetic models with local stores on/off. The original near-mean branch boundary and all coefficients remain fixed. Interior integer PMFs n<=1e6 are a research resource scope. Discovery; no fitted constant or selector.");
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

/// Apache Commons Math3.6.1 SaddlePointExpansion, commit16abfe5de688cc52fb0396e0609cb33044b15653.
/// The whole interior log-PMF program is fixed; arithmetic and elementary
/// providers are declared conditional adapters. See notices/CommonsMath.txt.
#[derive(Clone, Copy, Debug, Serialize)]
struct CommonsProfile { model: EvalModel, local_store: bool, log: usize, exponential: usize }

fn commons_local<F:Fpu>(v:Val,p:CommonsProfile)->Val {
    if p.local_store {Val::f64(v.to_f64::<F>())} else {v}
}
fn commons_log<F:Fpu>(v:Val,p:CommonsProfile)->Val {
    match p.log {
        0=>Val::f64(excel_log(v.to_f64::<F>())),
        1=>Val::f64(crt_log(v.to_f64::<F>())),
        2=>Val::Ext(research_log::<F>(&v.to_ext())),
        _=>unreachable!(),
    }
}
fn commons_stirling<F:Fpu>(z:f64,p:CommonsProfile)->Val {
    // Source constants, not values inferred from Excel. Half-integer inputs
    // below15 use the source table; the admitted count domain is integral.
    const HALF:[f64;31]=[0.0,0.1534264097200273452913848,0.0810614667953272582196702,
        0.0548141210519176538961390,0.0413406959554092940938221,0.03316287351993628748511048,
        0.02767792568499833914878929,0.02374616365629749597132920,0.02079067210376509311152277,
        0.01848845053267318523077934,0.01664469118982119216319487,0.01513497322191737887351255,
        0.01387612882307074799874573,0.01281046524292022692424986,0.01189670994589177009505572,
        0.01110455975820691732662991,0.010411265261972096497478567,0.009799416126158803298389475,
        0.009255462182712732917728637,0.008768700134139385462952823,0.008330563433362871256469318,
        0.007934114564314020547248100,0.007573675487951840794972024,0.007244554301320383179543912,
        0.006942840107209529865664152,0.006665247032707682442354394,0.006408994188004207068439631,
        0.006171712263039457647532867,0.005951370112758847735624416,0.005746216513010115682023589,
        0.005554733551962801371038690];
    if z<15.0 {return Val::f64(HALF[(2.0*z)as usize]);}
    let m=p.model;let z=Val::f64(z);let z2=commons_local::<F>(mul::<F>(m,z,z),p);
    let mut value=Val::f64(0.0008417508417508417508417508);
    for c in [0.000595238095238095238095238,0.00079365079365079365079365,
              0.00277777777777777777778,0.083333333333333333333] {
        value=sub::<F>(m,Val::f64(c),div::<F>(m,value,z2));
    }
    commons_local::<F>(div::<F>(m,value,z),p)
}
fn commons_deviance<F:Fpu>(x:f64,mu:Val,p:CommonsProfile)->Val {
    let m=p.model;let a=Val::f64(x);
    let d=sub::<F>(m,a,mu);let sum=add::<F>(m,a,mu);
    if d.to_f64::<F>().abs()<mul::<F>(m,Val::f64(0.1),sum).to_f64::<F>() {
        let d=commons_local::<F>(d,p);
        let mut v=commons_local::<F>(div::<F>(m,d,sum),p);
        let mut s1=commons_local::<F>(mul::<F>(m,v,d),p);
        let mut ej=commons_local::<F>(mul::<F>(m,mul::<F>(m,Val::f64(2.0),a),v),p);
        v=commons_local::<F>(mul::<F>(m,v,v),p);
        for j in 1..10000 {
            let s=s1;
            ej=commons_local::<F>(mul::<F>(m,ej,v),p);
            s1=commons_local::<F>(add::<F>(m,s,div::<F>(m,ej,Val::f64((2*j+1)as f64))),p);
            if s1==s {return s1;}
        }
        // ubs:ignore — abort the standalone research receipt on a resource failure; never omit a row or claim completeness.
        panic!("published deviance series did not converge within research bound");
    }
    commons_local::<F>(sub::<F>(m,add::<F>(m,mul::<F>(m,a,commons_log::<F>(div::<F>(m,a,mu),p)),mu),a),p)
}
fn commons_mass<F:Fpu>(x:f64,n:f64,prob:f64,p:CommonsProfile)->String {
    let m=p.model;let a=Val::f64(x);let total=Val::f64(n);let rest=Val::f64(n-x);
    // Public count/probability inputs and complement are stored binary64,
    // matching the Java method's supplied p/q argument boundary.
    let mut value=sub::<F>(m,sub::<F>(m,commons_stirling::<F>(n,p),commons_stirling::<F>(x,p)),commons_stirling::<F>(n-x,p));
    value=sub::<F>(m,value,commons_deviance::<F>(x,mul::<F>(m,total,Val::f64(prob)),p));
    value=commons_local::<F>(sub::<F>(m,value,commons_deviance::<F>(n-x,mul::<F>(m,total,Val::f64(1.0-prob)),p)),p);
    let factor=commons_local::<F>(div::<F>(m,mul::<F>(m,mul::<F>(m,Val::f64(2.0*std::f64::consts::PI),a),rest),total),p);
    value=commons_local::<F>(add::<F>(m,mul::<F>(m,Val::f64(-0.5),commons_log::<F>(factor,p)),value),p);
    let answer=match p.exponential {
        0=>excel_exp(value.to_f64::<F>()),
        1=>F::store(&research_exp::<F>(&value.to_ext()),CW),
        2=>F::store(&research_exp::<F>(&research_load(value.to_f64::<F>())),CW),
        _=>unreachable!(),
    };
    if !answer.is_finite() {"error:Num".into()}
    else {format!("0x{:016x}",if answer.is_subnormal()||answer==0.0{0}else{answer.to_bits()})}
}

pub fn commons_binomial<F:Fpu>(input:&str,output:&str) {
    let mut ps=Vec::new();
    for model in crate::dsl::WIDE_MODELS {for local_store in [true,false] {for log in 0..3 {for exponential in 0..3 {
        ps.push(CommonsProfile{model,local_store,log,exponential});
    }}}}
    assert_eq!(ps.len(),90);
    binomial_study::<F,_>(input,output,&ps,false,|x,n,p,profile|commons_mass::<F>(x,n,p,*profile),
        "Apache Commons Math3.6.1 commit16abfe5de688cc52fb0396e0609cb33044b15653 SaddlePointExpansion, Apache2; source fixed table and all asymptotic terms.",
        "Interior integer PMF only, n<=1e6 research resource bound. Fixed complete source log-PMF with5 arithmetic models,2local-store policies,3conditional log and3EXP adapters. Java FastMath itself is not identified. No endpoint/CDF/other producer exclusion; no fitted parameters or product change.");
}

#[derive(Clone,Copy,Serialize)]
struct BetaMassProgram {model:EvalModel,shape:usize}
pub fn beta_binomial<F:Fpu>(input:&str,output:&str) {
    let mut ps=Vec::new();for model in crate::dsl::WIDE_MODELS {for shape in 0..4 {ps.push(BetaMassProgram{model,shape});}}
    binomial_study::<F,_>(input,output,&ps,false,|x,n,p,profile| {
        let m=profile.model;let rest=n-x;let q=1.0-p;
        let leaf=Val::f64(if profile.shape<2 {crate::ported::special_math_common::bratio_brcomp(x,rest,p,q)}
            else{crate::ported::special_math_common::bratio_brcomp(x+1.0,rest+1.0,p,q)});
        let v=match profile.shape {
            0=>div::<F>(m,mul::<F>(m,leaf,Val::f64(n)),mul::<F>(m,Val::f64(x),Val::f64(rest))),
            1=>div::<F>(m,mul::<F>(m,leaf,div::<F>(m,Val::f64(n),Val::f64(x))),Val::f64(rest)),
            2=>div::<F>(m,leaf,mul::<F>(m,mul::<F>(m,Val::f64(n+1.0),Val::f64(p)),Val::f64(q))),
            3=>div::<F>(m,div::<F>(m,div::<F>(m,leaf,Val::f64(n+1.0)),Val::f64(p)),Val::f64(q)),
            _=>unreachable!(),
        }.to_f64::<F>();
        if !v.is_finite(){"error:Num".into()}else{format!("0x{:016x}",if v==0.0||v.is_subnormal(){0}else{v.to_bits()})}
    },"Beta factorial identity; current kept TOMS708 brcomp conditional leaf, original source unchanged except research visibility.",
      "Interior integer PMFs n<=1e6 only.20labels:4fixed beta-derivative identities/groupings,5outer arithmetic models. Existing brcomp and complement binary64 leaves; no complete binomial/CDF kernel identified or fitted.");
}

#[derive(Clone,Copy,Serialize)]
struct GammaDensityProgram {model:EvalModel,local_store:bool,log:usize,exponential:usize,normalizer:usize,shape:usize}
pub fn gamma_density<F:Fpu>(input:&str,output:&str) {
    let mut ps=Vec::new();for model in crate::dsl::WIDE_MODELS {for local_store in [true,false] {for log in 0..3 {for exponential in 0..3 {for normalizer in 0..2 {for shape in 0..2 {
        ps.push(GammaDensityProgram{model,local_store,log,exponential,normalizer,shape});
    }}}}}}
    assert_eq!(ps.len(),360);
    binomial_study::<F,_>(input,output,&ps,true,|x,a,b,p| {
        let m=p.model;let v=Val::f64;
        let local=|value:Val|if p.local_store{Val::f64(value.to_f64::<F>())}else{value};
        let log=|value:Val|match p.log {0=>v(crt_log(value.to_f64::<F>())),1=>v(excel_log(value.to_f64::<F>())),2=>Val::Ext(research_log::<F>(&value.to_ext())),_=>unreachable!()};
        let normalizer=if a==1.0{0.0}else if p.normalizer==0{crate::ported::special_math_common::ln_gamma(a)}else{gammaln_excel(a)};
        let reduced=local(div::<F>(m,v(x),v(b)));
        let alpha=sub::<F>(m,v(a),v(1.0));
        let value=if p.shape==0 {
            sub::<F>(m,sub::<F>(m,sub::<F>(m,mul::<F>(m,alpha,log(v(x))),reduced),v(normalizer)),mul::<F>(m,v(a),log(v(b))))
        }else{
            sub::<F>(m,sub::<F>(m,sub::<F>(m,mul::<F>(m,alpha,log(reduced)),reduced),v(normalizer)),log(v(b)))
        };
        let value=local(value);
        let output=match p.exponential {0=>crt_exp(value.to_f64::<F>()),1=>excel_exp(value.to_f64::<F>()),2=>F::store(&research_exp::<F>(&value.to_ext()),CW),_=>unreachable!()};
        crate::eval::outcome::<F>(v(output))
    },"Gamma density definition and equivalent scale normalization; current cleared MIT raw-log PDF, Lanczos normalizer or current worksheet GAMMALN as conditional leaves.",
      "Normal positive scalar x/shape/scale, PDF only.360 fixed labels:5coherent arithmetic modes,local stores on/off,3log and3EXP adapters,2normalizer leaves,2algebraic scale forms. Original strict label0 reproduces registered/frozen actuals. No CDF/admission/complete gamma producer identified; no fitted constant or selector.");
}

fn binomial_study<F:Fpu,P:Serialize>(input:&str,output:&str,ps:&[P],gamma:bool,prediction:impl Fn(f64,f64,f64,&P)->String,source:&str,limits:&str) {
    use std::collections::{BTreeMap,HashMap};
    let destination=Path::new(output);fs::create_dir(destination).unwrap();
    let rows=crate::lines::<Value>(input);
    let mut counts=vec![[0usize;5];ps.len()];let mut first=vec![None;ps.len()];
    let mut excluded=BTreeMap::<String,usize>::new();let mut cache=HashMap::new();let mut controls=0;
    for (index,row) in rows.iter().enumerate() {
        let surface=if gamma{matches!(row["function"].as_str(),Some("GAMMA.DIST"|"GAMMADIST"))}else{matches!(row["function"].as_str(),Some("BINOM.DIST"|"BINOMDIST"))};
        if !surface || !row["result_index"].is_null() {
            *excluded.entry("surface_or_selection".into()).or_default()+=1;continue;
        }
        // ubs:ignore — malformed frozen input must abort this standalone receipt, never silently leave its counted scope.
        let Some(a)=row["args"].as_array() else{panic!("missing original arguments");};
        if a.len()!=4 || !(a[3]==false || numeric(&a[3])==Some(0.0)) {
            *excluded.entry("cdf_or_shape".into()).or_default()+=1;continue;
        }
        let Some((x,n,p))=numeric(&a[0]).zip(numeric(&a[1])).zip(numeric(&a[2])).map(|((x,n),p)|(x,n,p)) else {
            *excluded.entry("non_numeric".into()).or_default()+=1;continue;
        };
        let admitted=if gamma{[x,n,p].iter().all(|v|v.is_normal()&&*v>0.0)}
            else{x.is_finite() && n.is_finite() && x==x.trunc() && n==n.trunc() && 0.0<x && x<n && n<=1e6 && 0.0<p && p<1.0};
        if !admitted {
            *excluded.entry(if gamma{"outside_normal_positive_domain"}else{"outside_interior_integer_domain_or_resource_bound"}.into()).or_default()+=1;continue;
        }
        let (baseline,predictions)=cache.entry((x.to_bits(),n.to_bits(),p.to_bits())).or_insert_with(|| {
            let args=a.iter().enumerate().map(|(i,v)|calc_verdict::decode_argument(v,&format!("args[{i}]")).unwrap()).collect();
            let (registered,_)=crate::invoke_control(row["function"].as_str().unwrap(),args);
            if !gamma {
                let mass=current::binom_mass(x,n,p,1.0-p);
                let direct=if !mass.is_finite(){"error:Num".into()}else{format!("0x{:016x}",if mass==0.0||mass.is_subnormal(){0}else{mass.to_bits()})};
                assert_eq!(direct,registered,"registered mass route control {}",row["id"]);
            }
            controls+=1;
            let predictions=ps.iter().map(|profile|prediction(x,n,p,profile)).collect::<Vec<_>>();
            if gamma{assert_eq!(predictions[0],registered,"strict original gamma density {}",row["id"]);}
            (registered,predictions)
        });
        assert_eq!(*baseline,row["actual"].as_str().unwrap(),"frozen registered numeric return {}",row["id"]);
        let expected=row["expected"].as_str().unwrap();let was=baseline==expected;
        for (id,got) in predictions.iter().enumerate() {
            let now=got==expected;
            for (slot,b) in counts[id].iter_mut().zip([true,now,was,now&&!was,was&&!now]){*slot+=usize::from(b);}
            if !now && first[id].is_none(){first[id]=Some(json!({"trace_index":index,"id":row["id"],"got":got,"expected":expected,"baseline":baseline}));}
        }
    }
    assert!(counts[0][0]>0);assert_eq!(counts[0][0]+excluded.values().sum::<usize>(),rows.len());
    let safe:Vec<_>=counts.iter().enumerate().filter(|(_,c)|c[3]>0&&c[4]==0).map(|(i,_)|i).collect();
    let result=json!({"complete":true,"input_rows":rows.len(),"scoped":counts[0][0],"unique_registered_controls":controls,
        "profiles":ps,"counts_all_exact_baseline_gains_losses":counts,"first_contradictions":first,"excluded":excluded,"zero_loss_improvements":safe,
        "source":source,"limits":limits});
    crate::json_write(destination.join("summary.json"),&result);
    println!("rows {} scoped {} controls {} zero-loss {:?}",rows.len(),counts[0][0],controls,safe);
}

#[cfg(test)]
mod commons_controls {
    use super::*;
    #[test]
    fn published_binomial_masses_follow_exact_dyadic_distribution() {
        let profile=CommonsProfile{model:EvalModel::Strict,local_store:true,log:0,exponential:0};
        for (n,x,numerator) in [(2.,1.,2u64),(4.,2.,6),(8.,3.,56),(10.,6.,210),(16.,8.,12870)] {
            let bits=commons_mass::<Hard>(x,n,0.5,profile);let answer=f64::from_bits(u64::from_str_radix(&bits[2..],16).unwrap());
            let exact=numerator as f64/2f64.powi(n as i32);
            assert!((answer-exact).abs()<=8.*f64::EPSILON*exact,"n={n} x={x}");
        }
    }
}

/// Only the final EXP return is retained as PC64; all other cleared mass
/// arithmetic stays fixed and its binary64 publication is checked below.
fn binomial_exp_return<F:Fpu>(x:f64,n:f64,p:f64,q:f64)->Val {
    let lc=if x==0.0 {
        if n==0.0{return Val::f64(1.0);}
        if p<0.1{-current::deviance(n,n*q)-n*p}else{n*excel_log(q)}
    }else if x==n{crate::ported::excel_numeric::excel_x87_mul(n,excel_log(p))}else{
        let lc=current::stirling_error(n)-current::stirling_error(x)-current::stirling_error(n-x)
            -current::deviance(x,n*p)-current::deviance(n-x,n*q);
        let lf=crt_log(2.0*core::f64::consts::PI*x*(n-x)/n);
        let half=F::div(&Val::f64(lf).to_ext(),&Val::f64(2.0).to_ext(),CW);
        let exponent=F::sub(&Val::f64(lc).to_ext(),&half,CW);
        return Val::Ext(research_exp::<F>(&exponent));
    };
    Val::Ext(research_exp::<F>(&Val::f64(lc).to_ext()))
}

/// Conditional private EXP return, not a reconstruction of Excel's full mass.
/// Named call-expression spill programs precede the all-stored control.
pub fn hypergeom_exp_returns<F:Fpu>(input:&str,output:&str){
    use calc_func::value::CalcValue;
    let destination=Path::new(output);fs::create_dir(destination).unwrap();let rows=crate::lines::<Value>(input);
    let mut counts=[[0_usize;5];5];let mut exclusions=std::collections::BTreeMap::<String,usize>::new();
    let mut predictions=Vec::new();let mut boundary=Vec::new();
    let store=|v:Val|Val::f64(v.to_f64::<F>());
    for row in &rows{
        let name=row["function"].as_str().unwrap();
        if !matches!(name,"HYPGEOM.DIST"|"HYPGEOMDIST")||!row["result_index"].is_null(){*exclusions.entry("surface_or_selection".into()).or_default()+=1;continue;}
        let args=row["args"].as_array().unwrap();
        let values:Option<Vec<_>>=args.iter().enumerate().map(|(i,v)|if i==4{v.as_bool().map(|x|if x{1.0}else{0.0}).or_else(||numeric(v))}else{numeric(v)}).collect();
        let Some(a)=values else{*exclusions.entry("non_numeric".into()).or_default()+=1;continue;};
        if a.len()!=if name=="HYPGEOMDIST"{4}else{5}||!a[..4].iter().all(|v|v.is_finite()&&*v>=0.0&&*v==v.trunc()){
            *exclusions.entry("shape_or_integer_domain".into()).or_default()+=1;continue;
        }
        let(s,n,k,total)=(a[0],a[1],a[2],a[3]);
        if n>total||k>total||n==0.0||n==total||s>n||s>k||n-s>total-k||total>10000.0||a.get(4).is_some_and(|x|*x!=0.0)
            ||!row["expected"].as_str().unwrap().starts_with("0x"){
            *exclusions.entry("outside_nontrivial_numeric_pmf_scope".into()).or_default()+=1;continue;
        }
        let p=n/total;let q=(total-n)/total;
        let leaves=[(s,k),(n-s,total-k),(n,total)].map(|(x,trials)|{
            let v=binomial_exp_return::<F>(x,trials,p,q);
            assert_eq!(v.to_f64::<F>().to_bits(),current::binom_mass(x,trials,p,q).to_bits(),"private EXP store control {}",row["id"]);
            v
        });
        let[x,y,z]=leaves;let m=EvalModel::ResidentPc64;
        let values=[
            div::<F>(m,store(mul::<F>(m,store(x),y)),z),
            div::<F>(m,mul::<F>(m,store(x),store(y)),z),
            mul::<F>(m,store(x),div::<F>(m,store(y),z)),
            mul::<F>(m,store(div::<F>(m,store(x),z)),y),
            div::<F>(EvalModel::Strict,mul::<F>(EvalModel::Strict,store(x),store(y)),store(z)),
        ];
        let baseline=crate::invoke_control(name,args.iter().map(|v|calc_verdict::decode_argument(v,"hypergeom return").unwrap()).collect()).0;
        assert_eq!(baseline,row["actual"],"registered frozen HYPGEOM baseline");
        let mut strings=Vec::new();for(i,v)in values.into_iter().enumerate(){
            let value=v.to_f64::<F>();let prediction=if value.is_finite(){calc_verdict::encode_result(&CalcValue::number(if value.is_subnormal(){0.0}else{value})).unwrap()}else{"refusal:nonfinite_return_program".into()};
            let old=baseline==row["expected"];let now=prediction==row["expected"];
            for(j,b)in[true,now,old,now&&!old,old&&!now].into_iter().enumerate(){counts[i][j]+=usize::from(b);}
            if i==4{assert_eq!(prediction,baseline,"all-stored HYPGEOM control");}
            strings.push(prediction);
        }
        let record=json!({"original":row,"predictions":strings,"private_exp_returns":leaves.map(|v|v.exact_encoding())});
        if k==0.0{boundary.push(record.clone());}predictions.push(record);
    }
    crate::json_write(destination.join("predictions.json"),&predictions);crate::json_write(destination.join("zero_population_success_controls.json"),&boundary);
    let summary=json!({"complete":true,"input_rows":rows.len(),"excluded":exclusions,"counts":counts,"count_order":["all","exact","baseline_exact","gains","losses"],
        "programs":["a stored before b; numerator stored before resident c; divide","stored a,b locals; resident numerator and c","a stored; b stored; b/c resident then times a","a stored; a/c quotient stored before resident b; multiply","all leaves and outer operations binary64 control"],
        "limits":"Conditional existing cleared Stirling/deviance/CRT-log bodies; only the final private EXP return remains PC64 instead of stored binary64. Every private return is asserted to store to the unchanged cleared binomial mass. The registered HYPGEOM baseline and all-stored control are separately asserted. These are call-expression/ABI hypotheses, not Excel private intermediates. Numeric integer PMFs only, total bounded10000 for research. All exclusions remain unscored obligations. No coefficient, cutoff, per-row mask or match-count ranking."});
    crate::json_write(destination.join("summary.json"),&summary);println!("{summary}");
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

#[cfg(gamma_exponent_audit)]
std::thread_local! {
    static GAMMA_EXPONENT_PROGRAM:std::cell::Cell<Option<(EvalModel,bool)>>=const{std::cell::Cell::new(None)};
    static GAMMA_TEMME_PROGRAM:std::cell::Cell<Option<(u8,bool)>>=const{std::cell::Cell::new(None)};
    static GAMMA_TEMME_VISITS:std::cell::Cell<usize>=const{std::cell::Cell::new(0)};
    static GAMMA_TEMME_LEAF:std::cell::Cell<Option<[u64;6]>>=const{std::cell::Cell::new(None)};
}
/// The unchanged Temme correction, with fixed binary64 or stored PC64 operations.
#[cfg(gamma_exponent_audit)]
pub fn gamma_temme_finish(l:f64,c:f64,w:f64,t:f64,a:f64,rt2pin:f64)->Option<(f64,f64)> {
    let (mask,software)=GAMMA_TEMME_PROGRAM.with(|p|p.get())?;
    GAMMA_TEMME_VISITS.with(|v|v.set(v.get()+1));
    GAMMA_TEMME_LEAF.with(|v|v.set(Some([l,c,w,t,a,rt2pin].map(f64::to_bits))));
    fn calculate<F:Fpu>(l:f64,c:f64,w:f64,t:f64,a:f64,rt2pin:f64,mask:u8)->(f64,f64) {
        let model=|bit:u8|if mask>=32{crate::dsl::WIDE_MODELS[(mask-32)as usize]}else if mask&(1u8<<bit)==0{EvalModel::Strict}else{EvalModel::StoredPc64};
        let v=Val::f64;
        let correction=div::<F>(model(1),mul::<F>(model(0),v(rt2pin),v(t)),v(a.sqrt()));
        let inner=if l>=1.0{add::<F>(model(2),v(w),correction)}else{sub::<F>(model(2),v(w),correction)};
        let answer=mul::<F>(model(3),v(c),inner).to_f64::<F>();
        let complement=add::<F>(model(4),v(0.5),sub::<F>(model(4),v(0.5),v(answer))).to_f64::<F>();
        if l>=1.0{(complement,answer)}else{(answer,complement)}
    }
    Some(if software{calculate::<Soft>(l,c,w,t,a,rt2pin,mask)}else{calculate::<Hard>(l,c,w,t,a,rt2pin,mask)})
}
/// Fixed evaluation of a*log(x)-x before its original EXP helper boundary.
#[cfg(gamma_exponent_audit)]
pub fn gamma_exponent(a:f64,x:f64)->f64 {
    let Some((model,software))=GAMMA_EXPONENT_PROGRAM.with(|p|p.get())else{return a*crt_log(x)-x;};
    fn calculate<F:Fpu>(a:f64,x:f64,model:EvalModel)->f64 {
        sub::<F>(model,mul::<F>(model,Val::f64(a),Val::f64(crt_log(x))),Val::f64(x)).to_f64::<F>()
    }
    if software{calculate::<Soft>(a,x,model)}else{calculate::<Hard>(a,x,model)}
}

#[cfg(gamma_exponent_audit)]
pub fn gamma_exponents(input:&str,output:&str,software:bool) {
    gamma_cdf_arithmetic(input,output,software,false);
}
#[cfg(gamma_exponent_audit)]
pub fn gamma_temme(input:&str,output:&str,software:bool) {
    gamma_cdf_arithmetic(input,output,software,true);
}
#[cfg(gamma_exponent_audit)]
fn gamma_cdf_arithmetic(input:&str,output:&str,software:bool,temme:bool) {
    let dest=Path::new(output);fs::create_dir(dest).unwrap();let rows=crate::lines::<Value>(input);
    let programs=if temme{(0u8..37).map(|mask|if mask<32{json!({"mask":mask})}else{json!({"coherent_model":crate::dsl::WIDE_MODELS[(mask-32)as usize]})}).collect::<Vec<_>>()}else{crate::dsl::WIDE_MODELS.iter().map(|p|json!(p)).collect()};
    let mut visited_rows=0;
    let mut counts=vec![[0usize;5];programs.len()];let mut first=vec![None;programs.len()];let mut excluded=std::collections::BTreeMap::<String,usize>::new();
    let mut file=BufWriter::new(fs::File::create_new(dest.join("changes.jsonl")).unwrap());
    for(index,row)in rows.iter().enumerate(){
        let name=row["function"].as_str().unwrap();let Some(a)=row["args"].as_array()else{*excluded.entry("shape".into()).or_default()+=1;continue;};
        let flag=|v:&Value|v==true||numeric(v).is_some_and(|v|v!=0.0);
        let args:Option<(f64,f64,bool)>=match name {
            "GAMMA.DIST"|"GAMMADIST" if a.len()==4&&flag(&a[3])=>{
                numeric(&a[0]).zip(numeric(&a[1])).zip(numeric(&a[2])).and_then(|((x,shape),scale)|
                    (x.is_finite()&&shape.is_finite()&&scale.is_finite()&&x>=0.0&&shape>0.0&&shape!=0.5&&shape!=1.0&&scale>0.0).then_some((shape,x/scale,false)))
            },
            "POISSON.DIST"|"POISSON" if a.len()==3&&flag(&a[2])=>{
                numeric(&a[0]).zip(numeric(&a[1])).and_then(|(k,mean)|
                    (k.is_finite()&&k>=1.0&&mean.is_finite()&&mean>0.0).then_some((k.trunc()+1.0,mean,true)))
            },
            _=>None,
        };
        let Some((shape,x,upper))=args.filter(|(shape,x,_)|shape.is_finite()&&x.is_finite())else{*excluded.entry("outside_gamma_prefactor_cdf_route".into()).or_default()+=1;continue;};
        let args=a.iter().enumerate().map(|(i,v)|calc_verdict::decode_argument(v,&format!("args[{i}]")).unwrap()).collect();let baseline=crate::invoke_control(name,args).0;assert_eq!(baseline,row["actual"],"frozen gamma exponent baseline {}",row["id"]);
        let evaluate=||{
            let value=if upper{crate::ported::special_math_common::regularized_gamma_q(shape,x)}else{crate::ported::special_math_common::regularized_gamma_p(shape,x)};
            if upper&&value.is_subnormal(){"0x0000000000000000".into()}else{crate::eval::outcome::<Hard>(Val::f64(value))}
        };
        GAMMA_EXPONENT_PROGRAM.with(|p|p.set(None));GAMMA_TEMME_PROGRAM.with(|p|p.set(None));assert_eq!(evaluate(),baseline,"unchanged gamma copy {}",row["id"]);
        let expected=row["expected"].as_str().unwrap();let was=baseline==expected;
        GAMMA_TEMME_VISITS.with(|v|v.set(0));
        for id in 0..programs.len(){
            if temme{GAMMA_TEMME_PROGRAM.with(|p|p.set(Some((id as u8,software))));}
            else{GAMMA_EXPONENT_PROGRAM.with(|p|p.set(Some((crate::dsl::WIDE_MODELS[id],software))));}
            let got=evaluate();let now=got==expected;
            for(slot,b)in counts[id].iter_mut().zip([true,now,was,now&&!was,was&&!now]){*slot+=usize::from(b);}
            if was&&!now&&first[id].is_none(){first[id]=Some(json!({"id":row["id"],"expected":expected,"prediction":got}));}
            if got!=baseline{writeln!(file,"{}",json!({"trace_index":index,"id":row["id"],"program":id,"expected":expected,"baseline":baseline,"prediction":got,"temme_leaf_l_c_w_t_a_rt2pin":if temme{GAMMA_TEMME_LEAF.with(|v|v.get()).map(|leaf|leaf.map(|bits|format!("0x{bits:016x}")))}else{None}})).unwrap();}
        }
        visited_rows+=usize::from(GAMMA_TEMME_VISITS.with(|v|v.get())>0);
    }
    GAMMA_EXPONENT_PROGRAM.with(|p|p.set(None));GAMMA_TEMME_PROGRAM.with(|p|p.set(None));assert!(counts[0][0]>0);assert_eq!(counts[0][0]+excluded.values().sum::<usize>(),rows.len());
    let result=json!({"complete":true,"input_rows":rows.len(),"programs":programs,"temme_visited_rows":visited_rows,"counts_all_exact_baseline_gains_losses":counts,"first_losses":first,"excluded":excluded,
        "source":"Existing retained TOMS654 gamma algorithms; IEEE754 and Intel SDM precision/store models. The CRT logarithm, EXP helpers, normalizer, branches and loops remain original conditional producers.","limits":if temme{"32 fixed Temme finish programs: binary64 or storedPC64 for correction multiply, divide, add/subtract, scaling and complement groups, plus five coherent arithmetic models. Original square-root and polynomial producers retained. No fitted coefficient/threshold/row selector or complete gamma identification."}else{"Five coherent evaluations of a*log(x)-x at the existing three helper sites; output crosses the original binary64 exponential boundary. No fitted coefficient/threshold/row selector or complete gamma identification."}});
    crate::json_write(dest.join("summary.json"),&result);println!("{}",result);
}
