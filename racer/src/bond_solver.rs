//! ExcelFinancialFunctions common.fs findRoot/Newton/bounds/bisection,
//! commit 97756a375769e45059c7d58fc51fc7e391bc4984, Apache-2.0.
//! Copyright 2010-2013 Luca Bolognese. Notice: notices/ExcelFinancialFunctions.txt.
//! Changed: explicit finite failures and coherent arithmetic models; the current
//! product's schedule/price/discount helpers are conditional numerical leaves.
use crate::{dsl::{EvalModel, Op, WIDE_MODELS}, eval::{Val, binop}, ported::excel_numeric::Fpu};
use serde_json::{Value, json};
use std::{fs, path::Path};

fn neg(a: Val) -> Val {
    match a { Val::F64(b) => Val::F64(b ^ (1 << 63)), Val::Ext(mut e) => { e.0[9] ^= 0x80; Val::Ext(e) } }
}
fn add<F:Fpu>(m:EvalModel,a:Val,b:Val)->Val { binop::<F>(Op::Add(0,1),m,a,b) }
fn sub<F:Fpu>(m:EvalModel,a:Val,b:Val)->Val { add::<F>(m,a,neg(b)) }
fn mul<F:Fpu>(m:EvalModel,a:Val,b:Val)->Val { binop::<F>(Op::Mul(0,1),m,a,b) }
fn div<F:Fpu>(m:EvalModel,a:Val,b:Val)->Val { binop::<F>(Op::Div(0,1),m,a,b) }
fn stored<F:Fpu>(a:Val)->Option<f64> { let v=a.to_f64::<F>(); v.is_finite().then_some(v) }

/// RATE's measured secant and stopping rule, tested on the original bond price
/// leaf. Named seeds, money normalization and signs are separate hypotheses.
pub fn rate_rule<F:Fpu>(input:&str,output:&str) { rate_rule_in::<F>(input,output,"rate"); }

pub fn rate_rule_in<F:Fpu>(input:&str, output:&str,variable:&str) {
    assert!(["rate","base","discount","discount_polynomial"].contains(&variable),"unknown secant variable");
    use crate::ported::bond_core_family::ResearchPrice;
    let rows=crate::lines::<Value>(input); assert!(!rows.is_empty(),"empty YIELD RATE-rule scope"); let v=Val::f64;
    let mut programs=Vec::new(); let mut survivors=Vec::new();
    let stops=if variable=="rate"{vec![("RATE",1e-7,1e-7)]}else{vec![("RATE",1e-7,1e-7),("step10",1e-10,f64::INFINITY),("both10",1e-10,1e-10),("value10",f64::INFINITY,1e-10)]};
    for (stop,xtol,ftol) in stops {for seed in ["coupon","zero","0.05","0.1","one"] { for delta in [0.001,0.0001,0.01] {
      for sign in [-1i32,0,1] {for normalization in [1.0,0.01] {for staged in [false,true] {for model in WIDE_MODELS {
        let id=programs.len();let mut first=None; let mut checked=0;
        for row in &rows {
            let a:Vec<f64>=row["args"].as_array().unwrap().iter().map(|x|
                f64::from_bits(u64::from_str_radix(x.as_str().unwrap().trim_start_matches("0x"),16).unwrap())).collect();
            let price=ResearchPrice::new(&a).unwrap(); assert!(price.count()>1);
            let to_yield=|x:f64|match variable {
                "rate"=>Some(x),"base"=>stored::<F>(mul::<F>(model,sub::<F>(model,v(x),v(1.0)),v(a[5]))),
                _=>stored::<F>(mul::<F>(model,sub::<F>(model,div::<F>(model,v(1.0),v(x)),v(1.0)),v(a[5]))),
            };
            let to_variable=|y:f64|match variable {
                "rate"=>Some(y),"base"=>stored::<F>(add::<F>(model,v(1.0),div::<F>(model,v(y),v(a[5])))),
                _=>stored::<F>(div::<F>(model,v(1.0),add::<F>(model,v(1.0),div::<F>(model,v(y),v(a[5]))))),
            };
            let f=|x:f64|{
                let pv=if variable=="discount_polynomial" {
                    let coupon=div::<F>(model,mul::<F>(model,v(100.0),v(a[2])),v(a[5]));
                    let mut sum=v(a[4]);
                    for _ in 0..price.count(){sum=mul::<F>(model,add::<F>(model,sum,coupon),v(x));}
                    stored::<F>(sum)
                }else{to_yield(x).and_then(|y|price.value(y,staged).ok())};
                pv.and_then(|p|stored::<F>(mul::<F>(model,sub::<F>(model,v(p),v(a[3])),v(normalization))))
            };
            let mut trace=Vec::new();
            let answer=(||->Option<f64>{
                let mut before=to_variable(match seed {"coupon"=>a[2],"zero"=>0.0,"0.05"=>0.05,"0.1"=>0.1,_=>1.0})?;
                let mut fb=f(before)?;if fb==0.0{return to_yield(before);}
                let direction=if sign==0 {(-fb).signum()} else {f64::from(sign)};
                let mut current=stored::<F>(add::<F>(model,v(before),v(delta*direction)))?;
                let mut fc=f(current)?;
                for iteration in 1..=200 {
                    let change=stored::<F>(sub::<F>(model,v(fc),v(fb)))?;
                    if change==0.0{return None;}
                    let next=stored::<F>(sub::<F>(model,v(current),div::<F>(model,
                        mul::<F>(model,v(fc),sub::<F>(model,v(current),v(before))),v(change))))?;
                    let step=stored::<F>(sub::<F>(model,v(next),v(current)))?.abs();
                    (before,fb)=(current,fc);current=next;fc=f(current)?;
                    trace.push(json!({"iteration":iteration,"x":format!("{:016x}",current.to_bits()),
                        "f":format!("{:016x}",fc.to_bits()),"step":format!("{:016x}",step.to_bits())}));
                    if trace.len()>4{trace.remove(0);}
                    if step<xtol && fc.abs()<ftol{return to_yield(current);}
                }
                None
            })();
            let actual=answer.map(|x|format!("0x{:016x}",x.to_bits())).unwrap_or_else(||"error:Num".into());
            checked+=1;
            if actual!=row["expected"].as_str().unwrap(){first=Some(json!({"row":row,"prediction":actual,"trace":trace}));break;}
        }
        if first.is_none(){survivors.push(id);}
        programs.push(json!({"id":id,"variable":variable,"stop":stop,"seed":seed,"delta":delta,"sign":sign,"normalization":normalization,
            "staged_price":staged,"model":model,"checked":checked,"first_contradiction":first}));
      }}}}
    }}}
    crate::json_write(output,&json!({"complete":true,"rows":rows.len(),"programs":programs,"survivors":survivors,
        "variable":variable,"limits":"RATE's simultaneous1e-7 tests and200-step limit; transformed-variable runs also test step-only, value-only and both1e-10 stops. Rate/base/discount or a direct discounted cash-flow Horner polynomial. Five named seeds, three deltas/signs, two money normalizations, original/staged price and five coherent arithmetic models. Every program refutes at its first exact contradiction. Product price and schedule remain conditional leaves; no ranking by counts or claim to exhaust solver families."}));
    println!("{} RATE-rule programs; {} survivors",programs.len(),survivors.len());
}

/// Bracket-width hypotheses suggested by the observed one-sided1e-10 cap.
/// The interpolation is ordinary linear or inverse quadratic; no coefficients
/// are fit. Keep every iterate of the first contradiction to test termination.
pub fn bracket_rule<F:Fpu>(input:&str, output:&str) {
    use crate::ported::bond_core_family::ResearchPrice;
    let rows=crate::lines::<Value>(input);assert!(!rows.is_empty(),"empty YIELD bracket scope");let v=Val::f64;
    let mut programs=Vec::new();let mut survivors=Vec::new();
    for (seed_low,seed_high) in [(0.0,1.0),(0.0,0.1),(0.0,0.05),(-0.5,1.0)] {
      for method in 0..6 {for tol in [1e-10,1e-8,1e-7] {for staged in [false,true] {for model in WIDE_MODELS {
        let id=programs.len();let mut first=None;let mut checked=0;
        for row in &rows {
            let a:Vec<f64>=row["args"].as_array().unwrap().iter().map(|x|f64::from_bits(
                u64::from_str_radix(x.as_str().unwrap().trim_start_matches("0x"),16).unwrap())).collect();
            let price=ResearchPrice::new(&a).unwrap();assert!(price.count()>1);
            let f=|x:f64|price.value(x,staged).ok().and_then(|p|stored::<F>(sub::<F>(model,v(p),v(a[3]))));
            let op=|kind:u8,a:f64,b:f64|stored::<F>(match kind {
                0=>add::<F>(model,v(a),v(b)),1=>sub::<F>(model,v(a),v(b)),
                2=>mul::<F>(model,v(a),v(b)),_=>div::<F>(model,v(a),v(b)),});
            let mut trace=Vec::new();
            let answer=(||->Option<f64>{
                let(mut low,mut high)=(seed_low,seed_high);let(mut fl,mut fh)=(f(low)?,f(high)?);
                if fl==0.0{return Some(low);}if fh==0.0{return Some(high);}
                for _ in 0..20 {
                    if fl>0.0 && fh<0.0 {break;}
                    if fl<0.0 {let next=op(1,low,op(1,high,low)?)?;
                        low=if next<=-a[5] {(-a[5]+low)*0.5}else{next};fl=f(low)?;}
                    else {high=op(0,high,op(1,high,low)?)?;fh=f(high)?;}
                    if fl==0.0{return Some(low);}if fh==0.0{return Some(high);}
                }
                if fl*fh>0.0{return None;}
                let(mut previous,mut fp)=(low,fl);
                for iteration in 0..100 {
                    if high-low<tol {return Some(low);}
                    let midpoint=||op(0,low,op(2,0.5,op(1,high,low)?)?);
                    let secant=||op(1,high,op(3,op(2,fh,op(1,high,low)?)?,op(1,fh,fl)?)?);
                    let mut x=if method==0 || (method==1 && iteration%2==0)
                        || (method==2 && iteration%2==1) {midpoint()?}
                    else if method>=4 && previous!=low && previous!=high && fl!=fp && fh!=fp {
                        let mut sum=0.0;
                        for(x,fx,fy,fz)in [(low,fl,fh,fp),(high,fh,fl,fp),(previous,fp,fl,fh)] {
                            let numerator=op(2,op(2,x,fy)?,fz)?;
                            let denominator=op(2,op(1,fx,fy)?,op(1,fx,fz)?)?;
                            sum=op(0,sum,op(3,numerator,denominator)?)?;
                        }
                        sum
                    } else {secant()?};
                    if x<=low || x>=high || (method==5 && (x-low).min(high-x)<0.01*(high-low)) {x=midpoint()?;}
                    let fx=f(x)?;trace.push(json!({"iteration":iteration,"low":format!("{:016x}",low.to_bits()),
                        "high":format!("{:016x}",high.to_bits()),"x":format!("{:016x}",x.to_bits()),"f":format!("{:016x}",fx.to_bits())}));
                    if fx==0.0{return Some(x);}
                    if fx>0.0 {(previous,fp)=(low,fl);(low,fl)=(x,fx);}
                    else {(previous,fp)=(high,fh);(high,fh)=(x,fx);}
                }
                None
            })();
            let actual=answer.map(|x|format!("0x{:016x}",x.to_bits())).unwrap_or_else(||"error:Num".into());checked+=1;
            if actual!=row["expected"].as_str().unwrap(){first=Some(json!({"id":row["id"],"args":row["args"],
                "expected":row["expected"],"prediction":actual,"trace":trace}));break;}
        }
        if first.is_none(){survivors.push(id);}
        programs.push(json!({"id":id,"low":seed_low,"high":seed_high,"method":method,"tolerance":tol,
            "staged_price":staged,"model":model,"checked":checked,"first_contradiction":first}));
      }}}}
    }
    crate::json_write(output,&json!({"complete":true,"rows":rows.len(),"programs":programs,"survivors":survivors,
        "methods":["bisection","bisection/secant","secant/bisection","secant","inverse quadratic","safeguarded inverse quadratic"],
        "limits":"Exploratory finite bracket programs suggested by a one-sided error cap; no claim these are Excel's routine. Original price/schedule conditional; first exact contradiction refutes. No fit by match count."}));
    println!("{} bracket programs; {} survivors",programs.len(),survivors.len());
}

/// Read fixed iteration counts on zero-coupon bonds, including the sixth-order
/// Newton/reused-derivative/Newton and Newton/Halley shapes.
pub fn fixed_iterations<F:Fpu>(input:&str,output:&str) {
    use crate::ported::bond_core_family::ResearchPrice;
    let rows:Vec<_>=crate::lines::<Value>(input).into_iter().filter(|r|
        r["args"][2]=="0x0000000000000000").collect();
    assert!(!rows.is_empty(),"empty zero-coupon fixed-iteration scope");
    let v=Val::f64;let mut programs=Vec::new();let mut survivors=Vec::new();
    for seed in [0.0,0.05,0.1,1.0] {for algorithm in 0..4 {for staged in [false,true] {for model in WIDE_MODELS {
        let mut first:Vec<Option<Value>>=vec![None;17];let mut checked=[0usize;17];
        for row in &rows {
            let a:Vec<f64>=row["args"].as_array().unwrap().iter().map(|x|f64::from_bits(
                u64::from_str_radix(x.as_str().unwrap().trim_start_matches("0x"),16).unwrap())).collect();
            let price=ResearchPrice::new(&a).unwrap();let n=price.count() as f64;
            let mut x:Option<f64>=Some(seed);let mut derivative=None;
            for iteration in 0..17 {
                let prediction=x.map(|x|format!("0x{:016x}",x.to_bits())).unwrap_or_else(||"error:Num".into());
                if first[iteration].is_none(){checked[iteration]+=1;
                    if prediction!=row["expected"].as_str().unwrap(){first[iteration]=Some(json!({"id":row["id"],
                        "args":row["args"],"expected":row["expected"],"prediction":prediction}));}}
                x=(||->Option<f64>{
                    let xx=x?;let px=price.value(xx,staged).ok()?;
                    let residual=stored::<F>(sub::<F>(model,v(px),v(a[3])))?;
                    let denominator=add::<F>(model,v(a[5]),v(xx));
                    let d=stored::<F>(neg(div::<F>(model,mul::<F>(model,v(n),v(px)),denominator)))?;
                    if algorithm!=1 || iteration%2==0 {derivative=Some(d);}
                    let d=derivative?;
                    let correction=if algorithm==2 || algorithm==3 && iteration%2==1 {
                        let second=div::<F>(model,mul::<F>(model,mul::<F>(model,v(n),v(n+1.0)),v(px)),
                            mul::<F>(model,denominator,denominator));
                        div::<F>(model,mul::<F>(model,mul::<F>(model,v(2.0),v(residual)),v(d)),
                            sub::<F>(model,mul::<F>(model,v(2.0),mul::<F>(model,v(d),v(d))),mul::<F>(model,v(residual),second)))
                    } else {div::<F>(model,v(residual),v(d))};
                    stored::<F>(sub::<F>(model,v(xx),correction))
                })();
            }
            if first.iter().all(Option::is_some){break;}
        }
        for(iteration,(first,checked))in first.into_iter().zip(checked).enumerate(){
            let id=programs.len();if first.is_none(){survivors.push(id);}
            programs.push(json!({"id":id,"seed":seed,"algorithm":algorithm,"iterations":iteration,
                "staged_price":staged,"model":model,"checked":checked,"first_contradiction":first}));
        }
    }}}}
    crate::json_write(output,&json!({"complete":true,"zero_coupon_rows":rows.len(),"programs":programs,
        "survivors":survivors,"algorithms":["Newton","derivative recomputed every second update","Halley","alternating Newton and Halley"],
        "limits":"Zero-coupon conditional price leaves and analytic derivatives only;0through16 updates, no convergence test or count selected by agreement percentage. Includes the low-order update combinations with local sixth-order error. These exact contradictions do not identify Excel's missing routine."}));
    println!("{} fixed-iteration programs; {} survivors",programs.len(),survivors.len());
}

/// All source constants and evaluation order are fixed before looking at gains.
fn root<F:Fpu>(f:impl Fn(f64)->Option<f64>,m:EvalModel)->Option<f64> {
    let precision=1e-6;
    let guess=0.05;
    let v=Val::f64;
    let mut x=guess;
    for count in 0..=21 {
        let fx=f(x)?;
        let above=f(stored::<F>(add::<F>(m,v(x),v(precision)))?)?;
        let below=f(stored::<F>(sub::<F>(m,v(x),v(precision)))?)?;
        let derivative=stored::<F>(div::<F>(m,sub::<F>(m,v(above),v(below)),mul::<F>(m,v(2.0),v(precision))))?;
        let next=stored::<F>(sub::<F>(m,v(x),div::<F>(m,v(fx),v(derivative))))?;
        let difference=stored::<F>(sub::<F>(m,v(next),v(x)))?;
        if difference.abs()<precision {
            if next>0.0 {return Some(next);}
            break;
        }
        if count>20 {break;}
        x=next;
    }
    let clamp=|value:f64| if value<=-1.0 {-1.0+precision}
        else if value>=f64::MAX {f64::MAX-precision} else {value};
    let mut low=clamp(stored::<F>(sub::<F>(m,v(guess),v(0.01)))?);
    let mut high=clamp(stored::<F>(add::<F>(m,v(guess),v(0.01)))?);
    let mut bounded=false;
    for _ in 1..60 {
        let lower=clamp(low); let upper=clamp(high);
        let product=stored::<F>(mul::<F>(m,v(f(lower)?),v(f(upper)?)))?;
        if product<=0.0 {low=lower; high=upper; bounded=true; break;}
        low=stored::<F>(add::<F>(m,v(lower),mul::<F>(m,v(1.6),sub::<F>(m,v(lower),v(upper)))))?;
        high=stored::<F>(add::<F>(m,v(upper),mul::<F>(m,v(1.6),sub::<F>(m,v(upper),v(lower)))))?;
    }
    if !bounded {return None;}
    for count in 0..=200 {
        if low==high {return None;}
        let fa=f(low)?;
        if fa.abs()<precision {return Some(low);}
        let fb=f(high)?;
        if fb.abs()<precision {return Some(high);}
        if count+1>200 || stored::<F>(mul::<F>(m,v(fa),v(fb)))?>0.0 {return None;}
        let mid=stored::<F>(add::<F>(m,v(low),mul::<F>(m,v(0.5),sub::<F>(m,v(high),v(low)))))?;
        let fm=f(mid)?;
        if fm.abs()<precision {return Some(mid);}
        let product=stored::<F>(mul::<F>(m,v(fa),v(fm)))?;
        if product<0.0 {high=mid;} else if product>0.0 {low=mid;} else {return None;}
    }
    None
}

/// Finite Newton hypotheses for independently observed on-coupon zero-coupon bonds.
/// Seeds, differences and stopping tolerances are written constants, never fitted.
pub fn zero_coupon_newton<F:Fpu>(input:&str, output:&str) {
    use crate::ported::bond_core_family::ResearchPrice;
    let rows=crate::lines::<Value>(input);
    let destination=Path::new(output);fs::create_dir(destination).unwrap();
    let mut programs=Vec::new();
    let mut derivatives=vec![(0,0.0)];
    for method in 1..=3 {for step in [1e-6,1e-8,1e-4] {derivatives.push((method,step));}}
    let mut survivors=Vec::new();let mut executions=0usize;
    for seed in [0.1,0.05,1.0,0.0] {for staged in [true,false] {
        for model in WIDE_MODELS {for &(method,step) in &derivatives {
            for tolerance in [1e-8,1e-6,1e-10,1e-4] {
                let id=programs.len();let mut first=None;let mut checked=0usize;
                for row in &rows {
                    assert_eq!(row["function"],"YIELD");
                    let a:Vec<f64>=row["args"].as_array().unwrap().iter().map(|a|
                        f64::from_bits(u64::from_str_radix(a.as_str().unwrap().trim_start_matches("0x"),16).unwrap())).collect();
                    assert_eq!(a.len(),7);assert_eq!(a[2],0.0);
                    let price=ResearchPrice::new(&a).unwrap();assert!(price.count()>1);
                    let mut trace=Vec::new();
                    let answer=(||->Option<f64>{
                        if price.is_par(a[3]) {return Some(0.0);}
                        let v=Val::f64;let mut x=seed;
                        for iteration in 0..100 {
                            let px=price.value(x,staged).ok()?;
                            let residual=stored::<F>(sub::<F>(model,v(px),v(a[3])))?;
                            if residual==0.0 {return Some(x);}
                            let derivative=if method==0 {
                                stored::<F>(neg(div::<F>(model,
                                    mul::<F>(model,v(price.count()as f64),v(px)),
                                    add::<F>(model,v(a[5]),v(x)))))?
                            } else {
                                let above=if method!=3 {
                                    price.value(stored::<F>(add::<F>(model,v(x),v(step)))?,staged).ok()?
                                }else{px};
                                let below=if method!=1 {
                                    price.value(stored::<F>(sub::<F>(model,v(x),v(step)))?,staged).ok()?
                                }else{px};
                                stored::<F>(div::<F>(model,
                                    sub::<F>(model,v(above),v(below)),
                                    v(if method==2 {2.0*step}else{step})))?
                            };
                            if derivative==0.0 {return None;}
                            let next=stored::<F>(sub::<F>(model,v(x),div::<F>(model,v(residual),v(derivative))))?;
                            trace.push(json!({"iteration":iteration,"x":format!("{:016x}",x.to_bits()),
                                "price":format!("{:016x}",px.to_bits()),"residual":format!("{:016x}",residual.to_bits()),
                                "derivative":format!("{:016x}",derivative.to_bits()),"next":format!("{:016x}",next.to_bits())}));
                            if trace.len()>3 {trace.remove(0);}
                            let difference=stored::<F>(sub::<F>(model,v(next),v(x)))?;
                            if difference.abs()<tolerance {return Some(next);}
                            x=next;
                        }
                        None
                    })();
                    let prediction=answer.map(|x|format!("0x{:016x}",x.to_bits())).unwrap_or_else(||"error:Num".into());
                    checked+=1;executions+=1;
                    if prediction!=row["expected"].as_str().unwrap() {
                        first=Some(json!({"original":row,"prediction":prediction,"last_iterations":trace}));break;
                    }
                }
                if first.is_none(){survivors.push(id);}
                programs.push(json!({"id":id,"seed":seed,"staged_price":staged,"model":model,
                    "derivative":method,"difference_step":step,"tolerance":tolerance,
                    "checked_until_first_contradiction":checked,"first_contradiction":first}));
            }
        }}
    }}
    crate::json_write(destination.join("summary.json"),&json!({"complete":true,"rows":rows.len(),
        "programs":programs,"survivors":survivors,"executions":executions,
        "method_order":["analytic zero-coupon derivative","forward difference","central difference","backward difference"],
        "limits":"Finite Newton programs only. Existing product date/price helpers are conditional leaves. Each named iterate and price result publishes binary64; expression operations use the stated coherent precision. Existing par shortcut retained. No fallback, derivative step, tolerance or seed fitted to counts; first contradiction refutes a complete hypothesis."}));
    println!("{} programs; {} survivors; {} executions",programs.len(),survivors.len(),executions);
}

pub fn run<F:Fpu>(input:&str, output:&str) {
    use crate::ported::bond_core_family::ResearchPrice;
    use calc_func::locale::{CalcEnvironment,DateSystem};
    let destination=Path::new(output);fs::create_dir(destination).unwrap();
    let rows=crate::lines::<Value>(input);
    let mut ranking=Vec::new(); let mut first_losses=Vec::new();
    for (model_index,&model) in WIDE_MODELS.iter().enumerate() {
        for staged in [false,true] {
            let id=model_index*2+usize::from(staged);
            let mut counts=[0_usize;6]; // all,exact,baseline,gains,losses,finite-program failure
            let mut contradictions=Vec::new();
            for row in &rows {
                assert_eq!(row["function"],"YIELD");
                let args:Vec<_>=row["args"].as_array().unwrap().iter()
                    .map(|a|calc_verdict::decode_argument(a,"YIELD solver study").unwrap()).collect();
                let Some(numbers)=args.iter().map(|a|a.as_number()).collect::<Option<Vec<_>>>() else {continue;};
                let environment=CalcEnvironment::dates(DateSystem::System1900);
                let baseline=crate::invoke_selected_control_in("YIELD",args,None,Some(&environment)).0;
                assert_eq!(baseline,row["actual"],"original registered baseline");
                let (prediction,failed)=if let Ok(price)=ResearchPrice::new(&numbers) {
                    if price.count()<=1 {(baseline.clone(),false)} else {
                        let f=|x|price.value(x,staged).ok().and_then(|a|stored::<F>(sub::<F>(model,Val::f64(a),Val::f64(numbers[3]))));
                        let answer=root::<F>(f,model);
                        (answer.map(|a|format!("0x{:016x}",a.to_bits())).unwrap_or_else(||"error:Num".into()),answer.is_none())
                    }
                } else {(baseline.clone(),false)};
                let now=prediction==row["expected"].as_str().unwrap();
                let old=baseline==row["expected"].as_str().unwrap();
                for (j,b) in [true,now,old,now&&!old,old&&!now,failed].into_iter().enumerate(){counts[j]+=usize::from(b);}
                if old&&!now&&contradictions.is_empty(){contradictions.push(json!({"original":row,"prediction":prediction}));}
            }
            ranking.push(json!({"id":id,"model":model,"staged_discount":staged,"counts":counts,
                "inclusive_failures_all_original_rows":rows.len()-counts[1]}));
            first_losses.push(json!({"id":id,"first_losses":contradictions}));
        }
    }
    crate::json_write(destination.join("first-losses.json").to_str().unwrap(),&first_losses);
    let scope=ranking[0]["counts"][0].as_u64().unwrap() as usize;
    assert!(scope>0);assert!(ranking.iter().all(|r|r["counts"][0]==scope));
    let result=json!({"complete":true,"program_labels":10,"rows":rows.len(),"scoped_numeric_rows":scope,"unscoped_rows_still_fail":rows.len()-scope,"count_order":["all","exact","baseline_exact","gains","losses","finite_program_failed"],"programs":ranking,
        "source":"ExcelFinancialFunctions97756a375769e45059c7d58fc51fc7e391bc4984 common.fs findRoot/newton/findBounds/bisection; Apache2.0, notice retained.",
        "limits":"Current original product schedule and price/discount helpers are conditional. Ten coherent arithmetic/discount labels; every local/helper boundary stores binary64. Published constants and stopping order remain fixed. Nonfinite intermediate/leaf errors yield a typed NUM hypothesis, not a claim to reproduce all F# exceptional behavior. Nonnumeric rows are outside this finite diagnostic scope. Single-coupon/domain baseline routes are preserved. No fitted constant, threshold or row-specific program selector; all answers are discovery."});
    crate::json_write(destination.join("summary.json").to_str().unwrap(),&result);println!("{result}");
}

/// Price-factor controls with the original product solver and caller branches.
/// The alternate A factor is the published day-count numerator already used
/// by DURATION; this does not alter or identify the solve/discount arithmetic.
pub fn components(input:&str,output:&str) {
    use crate::ported::bond_core_family::ResearchPrice;
    use calc_func::locale::{CalcEnvironment,DateSystem};
    let destination=Path::new(output);fs::create_dir(destination).unwrap();
    let rows=crate::lines::<Value>(input);
    let mut counts=[[0_usize;5];4];let mut first_losses=vec![None;4];let mut changes=Vec::new();
    let mut excluded=0;
    for row in &rows {
        let args:Vec<_>=row["args"].as_array().unwrap().iter()
            .map(|a|calc_verdict::decode_argument(a,"YIELD factor study").unwrap()).collect();
        let Some(numbers)=args.iter().map(|a|a.as_number()).collect::<Option<Vec<_>>>() else {excluded+=1;continue;};
        let environment=CalcEnvironment::dates(DateSystem::System1900);
        let baseline=crate::invoke_selected_control_in("YIELD",args,None,Some(&environment)).0;
        assert_eq!(baseline,row["actual"],"original registered factor baseline");
        for numerator_days in [false,true] {for staged in [false,true] {
            let id=usize::from(numerator_days)*2+usize::from(staged);
            let prediction=if let Ok(price)=ResearchPrice::new(&numbers) {
                if price.count()<=1 || price.is_par(numbers[3]) || !numbers[3].is_finite() || numbers[3]<=0.0 {baseline.clone()}
                else {price.original_root(numbers[3],staged,numerator_days)
                    .map(|x|format!("0x{:016x}",x.to_bits())).unwrap_or_else(|_|"error:Num".into())}
            } else {baseline.clone()};
            if id==0 {assert_eq!(prediction,baseline,"unchanged solver and date/discount leaves");}
            let now=prediction==row["expected"].as_str().unwrap();let old=baseline==row["expected"].as_str().unwrap();
            for (j,v) in [true,now,old,now&&!old,old&&!now].into_iter().enumerate(){counts[id][j]+=usize::from(v);}
            if old&&!now&&first_losses[id].is_none(){first_losses[id]=Some(json!({"original":row,"prediction":prediction}));}
            if baseline!=prediction {changes.push(json!({"program":id,"original":row,"prediction":prediction}));}
        }}
    }
    crate::json_write(destination.join("changed-predictions.json").to_str().unwrap(),&changes);
    let summary=json!({"complete":true,"program_labels":4,"original_rows":rows.len(),"excluded_still_fail":excluded,
        "count_order":["all","exact","baseline_exact","gains","losses"],"counts":counts,"first_losses":first_losses,
        "source":"ExcelFinancialFunctions97756a375769e45059c7d58fc51fc7e391bc4984 bonds.fs getPriceYieldFactors uses DaysBetween(...,Numerator); Apache2 notice retained. Current cleared days_between_num implements that public date rule.",
        "limits":"Four complete declared date-factor/discount combinations with the original product solve and unchanged caller domain, single-coupon and par identities. Original solver, schedule and discount remain conditional numerical leaves. Baseline program0 equals every original registered numeric return. No fitted constant/tolerance/selector or complete YIELD identification."});
    crate::json_write(destination.join("summary.json").to_str().unwrap(),&summary);println!("{summary}");
}

/// EFF getMatFactors uses the existing numerator-position day-count helper;
/// the priceMat association is the unchanged registered equation.
pub fn mat<F:Fpu>(input:&str,output:&str) {
    let destination=Path::new(output);fs::create_dir(destination).unwrap();
    let rows=crate::lines::<Value>(input);let mut counts=[[0_usize;5];10];
    let mut changes=Vec::new();let mut excluded=0;
    for row in &rows {
        assert_eq!(row["function"],"PRICEMAT");
        let args:Vec<_>=row["args"].as_array().unwrap().iter().map(|a|calc_verdict::decode_argument(a,"PRICEMAT study").unwrap()).collect();
        let Some(a)=args.iter().map(|a|a.as_number()).collect::<Option<Vec<_>>>() else {excluded+=1;continue;};
        let environment=calc_func::locale::CalcEnvironment::dates(calc_func::locale::DateSystem::System1900);
        let (baseline,numeric)=crate::invoke_selected_control_in("PRICEMAT",args,None,Some(&environment));
        assert_eq!(baseline,row["actual"],"original registered PRICEMAT return");
        for numerator in [false,true] {for (i,m)in WIDE_MODELS.into_iter().enumerate() {
            let id=usize::from(numerator)*5+i;let v=Val::f64;
            let prediction=if numeric.is_some() {
                let (b,dim,accrued_days)=crate::ported::bond_core_family::research_mat_factors(&a,numerator).unwrap();
                let compute=|| {
                    let dsm=stored::<F>(sub::<F>(m,v(dim),v(accrued_days)))?;
                    let future=stored::<F>(add::<F>(m,v(100.0),mul::<F>(m,mul::<F>(m,div::<F>(m,v(dim),v(b)),v(a[3])),v(100.0))))?;
                    let accrued=stored::<F>(mul::<F>(m,mul::<F>(m,div::<F>(m,v(accrued_days),v(b)),v(a[3])),v(100.0)))?;
                    let den=stored::<F>(add::<F>(m,v(1.0),mul::<F>(m,div::<F>(m,v(dsm),v(b)),v(a[4]))))?;
                    if den<=0.0 {return None;}
                    stored::<F>(sub::<F>(m,div::<F>(m,v(future),v(den)),v(accrued)))
                };
                compute().map(|x|format!("0x{:016x}",x.to_bits())).unwrap_or_else(||"error:Num".into())
            } else {baseline.clone()};
            if id==0 {assert_eq!(prediction,baseline,"original day factors and strict association");}
            let now=prediction==row["expected"].as_str().unwrap();let old=baseline==row["expected"].as_str().unwrap();
            for (j,v)in [true,now,old,now&&!old,old&&!now].into_iter().enumerate(){counts[id][j]+=usize::from(v);}
            if prediction!=baseline {changes.push(json!({"program":id,"original":row,"prediction":prediction}));}
        }}
    }
    crate::json_write(destination.join("changed-predictions.json").to_str().unwrap(),&changes);
    let result=json!({"complete":true,"program_labels":10,"original_rows":rows.len(),"excluded_still_fail":excluded,
        "count_order":["all","exact","baseline_exact","gains","losses"],"counts":counts,
        "source":"ExcelFinancialFunctions97756a375769e45059c7d58fc51fc7e391bc4984 bonds.fs getMatFactors/priceMat; Apache2, notice retained. Numerator-position helper is the existing cleared product date rule.",
        "programs":"Two fixed date-factor producers times five coherent arithmetic precisions. Existing equation association; named local results store binary64. No fitted coefficient, tolerance or row-specific selector.",
        "limits":"Registered domain/error branches are conditional and retained. Date formulas/elementary arithmetic outside this declared finite scope are not identified. Complete changes versus each checked original baseline are retained."});
    crate::json_write(destination.join("summary.json").to_str().unwrap(),&result);println!("{result}");
}

/// Apache POI 3.17 Irr.java, ab700dc2ae3a880f66a10075cbf1a502da928afb.
/// Copyright 2003-2017 The Apache Software Foundation; Apache-2.0.
/// Changed: coherent precision models and explicit finite failure, without
/// changing the published twenty iterations, tolerance or recurrence order.
/// Notice and license: notices/Apache-POI.txt.
fn poi_irr<F:Fpu>(values:&[f64],guess:f64,m:EvalModel)->Option<f64> {
    let v=Val::f64;
    let mut x=guess;
    for _ in 0..20 {
        let factor=stored::<F>(add::<F>(m,v(1.0),v(x)))?;
        let mut denominator=factor;
        let mut value=values[0];
        let mut derivative=0.0;
        for (k,&cash) in values.iter().enumerate().skip(1) {
            value=stored::<F>(add::<F>(m,v(value),div::<F>(m,v(cash),v(denominator))))?;
            denominator=stored::<F>(mul::<F>(m,v(denominator),v(factor)))?;
            derivative=stored::<F>(sub::<F>(m,v(derivative),div::<F>(m,mul::<F>(m,v(k as f64),v(cash)),v(denominator))))?;
        }
        let next=stored::<F>(sub::<F>(m,v(x),div::<F>(m,v(value),v(derivative))))?;
        if stored::<F>(sub::<F>(m,v(next),v(x)))?.abs()<=1e-7 {return Some(next);}
        x=next;
    }
    None
}

pub fn periodic<F:Fpu>(input:&str,output:&str) {
    let destination=Path::new(output);fs::create_dir(destination).unwrap();
    let rows=crate::lines::<Value>(input);
    let mut counts=[[0_usize;5];10];let mut losses=vec![None;10];let mut changes=Vec::new();
    let mut excluded=0;
    for row in &rows {
        assert_eq!(row["function"],"IRR");
        let sources=row["args"].as_array().unwrap();
        let mut values=Vec::new();let mut numeric=true;
        fn collect(argument:&Value,values:&mut Vec<f64>)->bool {
            if let Some(array)=argument.as_array() {array.iter().all(|v|collect(v,values))}
            else if let Some(number)=calc_verdict::decode_argument(argument,"IRR numeric study").unwrap().as_number() {
                values.push(number);true
            } else {false}
        }
        if !sources[0].is_array() || !collect(&sources[0],&mut values) {numeric=false;}
        let guess=if sources.len()==1 {Some(0.1)} else {
            calc_verdict::decode_argument(&sources[1],"IRR guess").unwrap().as_number()
        };
        let Some(guess)=guess else {excluded+=1;continue;};
        if !numeric || values.is_empty() || values.iter().any(|v|!v.is_finite()) || !guess.is_finite() {
            excluded+=1;continue;
        }
        let args=sources.iter().map(|a|calc_verdict::decode_argument(a,"IRR baseline").unwrap()).collect();
        let baseline=crate::invoke_control("IRR",args).0;
        assert_eq!(baseline,row["actual"],"original registered IRR baseline");
        for (i,m) in WIDE_MODELS.into_iter().enumerate() {for keep_two in [false,true] {
            let id=i*2+usize::from(keep_two);
            let prediction=if keep_two && values.len()==2 {baseline.clone()} else {
                poi_irr::<F>(&values,guess,m).map(|x|format!("0x{:016x}",x.to_bits())).unwrap_or_else(||"error:Num".into())
            };
            let now=prediction==row["expected"].as_str().unwrap();let old=baseline==row["expected"].as_str().unwrap();
            for (j,v) in [true,now,old,now&&!old,old&&!now].into_iter().enumerate(){counts[id][j]+=usize::from(v);}
            if old&&!now&&losses[id].is_none(){losses[id]=Some(json!({"original":row,"prediction":prediction}));}
            if now&&!old {changes.push(json!({"program":id,"original":row,"prediction":prediction}));}
        }}
    }
    crate::json_write(destination.join("gains.json").to_str().unwrap(),&changes);
    let result=json!({"complete":true,"program_labels":10,"original_rows":rows.len(),"excluded_still_fail":excluded,
        "count_order":["all","exact","baseline_exact","gains","losses"],"counts":counts,"first_losses":losses,
        "source":"Apache POI3.17 Irr.java ab700dc2ae3a880f66a10075cbf1a502da928afb, Apache2; full source/license/notice retained.",
        "programs":"Five coherent arithmetic models times two explicit domains: published complete recurrence or current two-cashflow route retained as a conditional leaf. Each declared local boundary stores binary64.",
        "limits":"Original numeric array readings and optional numeric guesses only; skipped text/error/reference obligations remain failures. Nonfinite intermediates map to NUM, without claiming all Java exceptional semantics. No fitted stopping constant, coefficient or row selector; discovery answers; no complete IRR identification."});
    crate::json_write(destination.join("summary.json").to_str().unwrap(),&result);println!("{result}");
}

/// Fixed NPV operation trees inside the measured discount-factor secant.
/// This is original numerical research; the older POI trial above stays separate.
fn irr_discount_program<F:Fpu>(flows:&[f64],guess:f64,kind:usize,reverse:bool,m:EvalModel,normalization:&str)->Option<f64> {
    let v=Val::f64;
    let scale=match normalization {"first"=>flows.iter().find(|x|**x!=0.0)?.abs(),"maximum"=>flows.iter().map(|x|x.abs()).fold(0.0,f64::max),_=>1.0};
    let normalized=flows.iter().map(|x|stored::<F>(div::<F>(m,v(*x),v(scale)))).collect::<Option<Vec<_>>>()?;
    let flows=&normalized;
    let npv=|q:f64|->Option<f64>{
        let base=stored::<F>(div::<F>(m,v(1.0),v(q)))?;
        if kind==6 {
            let mut sum=v(0.0);
            for &flow in flows.iter().rev() {sum=add::<F>(m,mul::<F>(m,sum,v(q)),v(flow));}
            return stored::<F>(sum);
        }
        let mut powers=vec![1.0;flows.len()];
        if kind==4 || kind==5 {
            for k in 1..flows.len(){powers[k]=stored::<F>(mul::<F>(m,v(powers[k-1]),v(if kind==4{base}else{q})))?;}
        }
        let mut sum=v(0.0);
        for index in 0..flows.len() {
            let k=if reverse{flows.len()-index-1}else{index};
            let power=match kind {
                0|2=>crate::ported::power_fn::power_kernel(if kind==0{base}else{q},k as f64).ok()?,
                1|3=>crate::ported::excel_numeric::crt_pow(if kind==1{base}else{q},k as f64),
                4|5=>powers[k], _=>unreachable!(),
            };
            if !power.is_finite(){return None;}
            let term=if kind==0||kind==1||kind==4{div::<F>(m,v(flows[k]),v(power))}else{mul::<F>(m,v(flows[k]),v(power))};
            sum=add::<F>(m,sum,term);
        }
        stored::<F>(sum)
    };
    let mut previous=stored::<F>(div::<F>(m,v(1.0),add::<F>(m,v(1.0),v(guess))))?;
    let mut fp=npv(previous)?;
    if fp==0.0{return Some(guess);}
    let mut current=stored::<F>(add::<F>(m,v(previous),v(0.001_f64.copysign(fp))))?;
    let mut fc=npv(current)?;
    for _ in 0..200 {
        let change=stored::<F>(sub::<F>(m,v(fc),v(fp)))?;
        if change==0.0{return None;}
        let next=stored::<F>(sub::<F>(m,v(current),div::<F>(m,mul::<F>(m,v(fc),sub::<F>(m,v(current),v(previous))),v(change))))?;
        let step=stored::<F>(sub::<F>(m,v(next),v(current)))?.abs();
        (previous,fp)=(current,fc);current=next;fc=npv(current)?;
        if step<1e-7 && fc.abs()<1e-7{return stored::<F>(sub::<F>(m,div::<F>(m,v(1.0),v(current)),v(1.0)));}
    }
    None
}

pub fn irr_discount_trees<F:Fpu>(input:&str,output:&str) {
    let rows=crate::lines::<Value>(input);let mut scope=Vec::new();let mut excluded=0;
    fn collect(a:&Value,out:&mut Vec<f64>)->bool {
        if let Some(a)=a.as_array(){a.iter().all(|x|collect(x,out))}
        else if let Some(x)=calc_verdict::decode_argument(a,"IRR tree").unwrap().as_number(){out.push(x);x.is_finite()}
        else{false}
    }
    for row in &rows {
        assert_eq!(row["function"],"IRR");let args=row["args"].as_array().unwrap();let mut flows=Vec::new();
        let guess=if args.len()==1{Some(0.1)}else{calc_verdict::decode_argument(&args[1],"guess").unwrap().as_number()};
        if !args[0].is_array()||!collect(&args[0],&mut flows)||flows.is_empty()||guess.is_none_or(|x|!x.is_finite())||!flows.iter().any(|x|*x<0.0)||!flows.iter().any(|x|*x>0.0){excluded+=1;continue;}
        scope.push((row,flows,guess.unwrap()));
    }
    assert!(!scope.is_empty(),"empty IRR tree scope");
    let mut controlled=0;
    for (row,flows,guess) in &scope {
        let args=row["args"].as_array().unwrap().iter().map(|x|calc_verdict::decode_argument(x,"IRR control").unwrap()).collect();
        let registered=crate::invoke_control("IRR",args).0;
        assert_eq!(registered,row["actual"],"frozen original IRR registered route");
        if registered!="error:Ref" {
            let got=irr_discount_program::<F>(flows,*guess,0,false,EvalModel::StoredPc64,"none")
                .map(|x|format!("0x{:016x}",x.to_bits())).unwrap_or_else(||"error:Num".into());
            assert_eq!(got,registered,"original whole stored discount-secant control");controlled+=1;
        }
    }
    let mut programs=Vec::new();
    for normalization in ["none","first","maximum"] {for kind in 0..7 {for reverse in [false,true] {if kind==6&&reverse{continue;}for model in WIDE_MODELS {
        let mut first=None;let mut checked=0;
        for (row,flows,guess)in &scope {
            let prediction=irr_discount_program::<F>(flows,*guess,kind,reverse,model,normalization)
                .map(|x|format!("0x{:016x}",x.to_bits())).unwrap_or_else(||"error:Num".into());
            checked+=1;
            if prediction!=row["expected"].as_str().unwrap(){first=Some(json!({"original":row,"prediction":prediction}));break;}
        }
        programs.push(json!({"kind":kind,"reverse":reverse,"normalization":normalization,"model":model,"checked":checked,"first_contradiction":first}));
    }}}}

    let result=json!({"complete":true,"rows":rows.len(),"scoped":scope.len(),"excluded_still_fail":excluded,"registered_whole_controls":controlled,
        "programs":programs,"survivors":programs.iter().filter(|p|p["first_contradiction"].is_null()).collect::<Vec<_>>(),
        "order":"Worksheet POWER quotient, CRT-power quotient, direct-q worksheet POWER product, direct-q CRT product, base-power recurrence, q-power recurrence, reverse Horner; forward/reverse sums; five coherent word policies; original flows or all flows divided by the first nonzero absolute value or maximum absolute value. Secant seeds and 1e-7 step-and-value stop fixed to the existing measured program.",
        "limits":"All original finite numeric arrays with both signs; rectangular arrays are a diagnostic flattened input, not an admitted application change. Complete trees are rejected by first exact contradiction, without match-percentage selection. No new Excel answers or uniquely recovered IRR callback."});
    crate::json_write(output,&result);println!("{result}");
}

/// The published discount-security equations in the registered bodies. These
/// variants change coherent arithmetic precision only; date rules are fixed.
#[cfg(closed_security_audit)]
fn security_value<F:Fpu>(function:&str,a:&[f64],days:f64,year:f64,m:EvalModel)->Option<f64> {
    let v=Val::f64;
    let scale=div::<F>(m,v(year),v(days));
    let result=match function {
        "DISC"=>mul::<F>(m,sub::<F>(m,v(1.0),div::<F>(m,v(a[2]),v(a[3]))),scale),
        "INTRATE"=>mul::<F>(m,div::<F>(m,sub::<F>(m,v(a[3]),v(a[2])),v(a[2])),scale),
        "RECEIVED"=>{
            let den=sub::<F>(m,v(1.0),div::<F>(m,mul::<F>(m,v(a[3]),v(days)),v(year)));
            if den.to_f64::<F>()<=0.0 {return None;}
            div::<F>(m,v(a[2]),den)
        },
        "PRICEDISC"=>sub::<F>(m,v(a[3]),div::<F>(m,mul::<F>(m,mul::<F>(m,v(a[2]),v(a[3])),v(days)),v(year))),
        "TBILLYIELD"=>mul::<F>(m,div::<F>(m,sub::<F>(m,v(100.0),v(a[2])),v(a[2])),scale),
        _=>panic!("undeclared security function"), // ubs:ignore -- the caller asserts the fixed study family before invoking this private helper.
    };
    stored::<F>(result)
}

#[cfg(closed_security_audit)]
fn bill_value<F:Fpu>(days:f64,discount:f64,m:EvalModel,staged_sqrt:bool)->Option<f64> {
    let v=Val::f64;
    if days<=182.0 {
        let den=stored::<F>(sub::<F>(m,v(360.0),mul::<F>(m,v(discount),v(days))))?;
        if den<=0.0 {return None;}
        return stored::<F>(div::<F>(m,mul::<F>(m,v(365.0),v(discount)),v(den)));
    }
    let price=stored::<F>(sub::<F>(m,v(100.0),div::<F>(m,mul::<F>(m,mul::<F>(m,v(discount),v(100.0)),v(days)),v(360.0))))?;
    if price<=0.0 {return None;}
    let year=if days==366.0 {366.0} else {365.0};
    let t=stored::<F>(div::<F>(m,v(days),v(year)))?;
    let term=stored::<F>(sub::<F>(m,v(1.0),div::<F>(m,v(1.0),div::<F>(m,v(price),v(100.0)))))?;
    let disc=stored::<F>(sub::<F>(m,mul::<F>(m,v(t),v(t)),mul::<F>(m,sub::<F>(m,mul::<F>(m,v(2.0),v(t)),v(1.0)),v(term))))?;
    if disc<0.0 {return None;}
    let root=if staged_sqrt {
        F::store(&F::sqrt(&crate::ported::excel_numeric::research_load(disc),crate::eval::CW),crate::eval::CW)
    } else {disc.sqrt()};
    stored::<F>(div::<F>(m,add::<F>(m,v(-t),v(root)),sub::<F>(m,v(t),v(0.5))))
}

#[cfg(closed_security_audit)]
pub fn bill_equivalent<F:Fpu>(input:&str,output:&str) {
    let destination=Path::new(output);fs::create_dir(destination).unwrap();
    let rows=crate::lines::<Value>(input);let mut counts=[[0_usize;5];10];
    let mut changes=Vec::new();let mut excluded=0;
    for row in &rows {
        assert_eq!(row["function"],"TBILLEQ");
        let args:Vec<_>=row["args"].as_array().unwrap().iter().map(|a|calc_verdict::decode_argument(a,"TBILLEQ study").unwrap()).collect();
        let Some(a)=args.iter().map(|a|a.as_number()).collect::<Option<Vec<_>>>() else {excluded+=1;continue;};
        let environment=calc_func::locale::CalcEnvironment::dates(calc_func::locale::DateSystem::System1900);
        let (baseline,numeric)=crate::invoke_selected_control_in("TBILLEQ",args,None,Some(&environment));
        assert_eq!(baseline,row["actual"],"original registered TBILLEQ return");
        for (i,m)in WIDE_MODELS.into_iter().enumerate() {for staged in [false,true] {
            let id=2*i+usize::from(staged);
            let prediction=if numeric.is_some() {
                let (days,_)=crate::ported::discount_bill_yearfrac_family::research_security_days(&a,true).unwrap();
                bill_value::<F>(days,a[2],m,staged).map(|x|format!("0x{:016x}",x.to_bits())).unwrap_or_else(||"error:Num".into())
            } else {baseline.clone()};
            if id==0 {assert_eq!(prediction,baseline,"strict original bill equation and SQRT");}
            let now=prediction==row["expected"].as_str().unwrap();let old=baseline==row["expected"].as_str().unwrap();
            for (j,v)in [true,now,old,now&&!old,old&&!now].into_iter().enumerate(){counts[id][j]+=usize::from(v);}
            if prediction!=baseline {changes.push(json!({"program":id,"original":row,"prediction":prediction}));}
        }}
    }
    crate::json_write(destination.join("changed-predictions.json").to_str().unwrap(),&changes);
    let result=json!({"complete":true,"program_labels":10,"original_rows":rows.len(),"excluded_still_fail":excluded,
        "count_order":["all","exact","baseline_exact","gains","losses"],"counts":counts,
        "source":"The existing published short-bill discount equation and long-bill quadratic identity; IEEE754 and shared x87 PC64 SQRT. No coefficient or cutoff changed.",
        "programs":"Five coherent arithmetic modes times two fixed SQRT leaves. Named local values store binary64; all original equation associations preserved.",
        "limits":"Date helper and registered error branches are conditional. Complete baseline-change predictions retained. No fitted tolerance, coefficient or per-input selector, and no whole-family exclusion claim."});
    crate::json_write(destination.join("summary.json").to_str().unwrap(),&result);println!("{result}");
}

#[cfg(closed_security_audit)]
pub fn securities<F:Fpu>(input:&str,output:&str) {
    use std::collections::BTreeMap;
    let destination=Path::new(output);fs::create_dir(destination).unwrap();
    let rows=crate::lines::<Value>(input);
    let mut counts:BTreeMap<String,[[usize;5];5]>=BTreeMap::new();
    let mut losses=BTreeMap::new();let mut changes=Vec::new();let mut excluded=0;
    for row in &rows {
        let name=row["function"].as_str().unwrap();
        assert!(["DISC","INTRATE","RECEIVED","PRICEDISC","TBILLYIELD"].contains(&name));
        let args:Vec<_>=row["args"].as_array().unwrap().iter().map(|a|calc_verdict::decode_argument(a,"security study").unwrap()).collect();
        let Some(numbers)=args.iter().map(|a|a.as_number()).collect::<Option<Vec<_>>>() else {excluded+=1;continue;};
        let environment=calc_func::locale::CalcEnvironment::dates(calc_func::locale::DateSystem::System1900);
        let (baseline,numeric)=crate::invoke_selected_control_in(name,args,None,Some(&environment));
        assert_eq!(baseline,row["actual"],"original registered security return");
        for (id,m) in WIDE_MODELS.into_iter().enumerate() {
            let prediction=if numeric.is_some() {
                let (days,year)=crate::ported::discount_bill_yearfrac_family::research_security_days(&numbers,name=="TBILLYIELD").unwrap();
                security_value::<F>(name,&numbers,days,year,m)
                    .map(|x|format!("0x{:016x}",x.to_bits())).unwrap_or_else(||"error:Num".into())
            } else {baseline.clone()};
            if id==0 {assert_eq!(prediction,baseline,"strict published equation equals registered return");}
            let now=prediction==row["expected"].as_str().unwrap();let old=baseline==row["expected"].as_str().unwrap();
            for (j,v) in [true,now,old,now&&!old,old&&!now].into_iter().enumerate(){counts.entry(name.into()).or_insert([[0;5];5])[id][j]+=usize::from(v);}
            if old&&!now {losses.entry(format!("{name}:{id}")).or_insert(json!({"original":row,"prediction":prediction}));}
            if prediction!=baseline {changes.push(json!({"function":name,"program":id,"original":row,"prediction":prediction}));}
        }
    }
    crate::json_write(destination.join("changed-predictions.json").to_str().unwrap(),&changes);
    let result=json!({"complete":true,"programs_per_function":5,"original_rows":rows.len(),"excluded_still_fail":excluded,
        "count_order":["all","exact","baseline_exact","gains","losses"],"models":WIDE_MODELS,"counts":counts,"first_losses":losses,
        "source":"Published discount-security equations, as in cleared registered bodies; IEEE754/Intel arithmetic precision. Original MIT implementation; no coefficients or thresholds fitted.",
        "limits":"Five coherent arithmetic models retain the equation association and original date helper. Each program applies uniformly within each named function. Registered error branches unchanged; nonnumeric rows remain excluded failures. Original recorded answers are discovery, not proof about unmeasured inputs."});
    crate::json_write(destination.join("summary.json").to_str().unwrap(),&result);println!("{result}");
}

#[cfg(test)]
mod tests {
    use super::*;
    /// Native-only: the research bridge uses native instruction controls.
    #[test]
    fn fixed_newton_and_fallback_control_roots() {
        use crate::ported::excel_numeric::Hard;
        for m in WIDE_MODELS {
            let result=root::<Hard>(|x|Some(x-0.125),m).unwrap();
            assert!((result-0.125).abs()<1e-6);
            // A negative root forces the source's opposite-sign fallback.
            let result=root::<Hard>(|x|Some(x+0.125),m).unwrap();
            assert!((result+0.125).abs()<1e-6);
        }
    }

    /// Native-only: the research precision models use native instructions.
    #[test]
    fn published_irr_recurrence_control() {
        for model in WIDE_MODELS {
            let value=poi_irr::<crate::ported::excel_numeric::Hard>(&[-100.0,121.0],0.1,model).unwrap();
            assert!((value-0.21).abs()<1e-7);
            assert!(poi_irr::<crate::ported::excel_numeric::Hard>(&[100.0,100.0],0.1,model).is_none());
        }
    }
}
