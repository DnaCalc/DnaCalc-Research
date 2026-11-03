//! Trace the existing registered and workbook routes, without selecting an alternative kernel.
use crate::ported::complex_family::{ParsedComplex, format_complex_text, parse_complex_text};
use crate::{Given, json_write, ported::financial_time_value_family as tv};
use calc_func::locale::{CalcEnvironment, DateSystem};
use calc_func::{EvalServices, FunctionCallScratch, Notation, catalogue, invoke};
use calc_host::value::{CalcValue, CellKey, CoreRef, NumericMode};
use calc_host::{CellEdit, CellInput, CellText, FiniteNumber, ProjectedValue, SkinIntent};
use dna_calc_worker::workbook::Workbook;
use serde_json::{Value, json};
use std::{fs, path::Path};

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
        let ProjectedValue::Calculated { value, .. } = cell.value() else {
            panic!("question prediction was not calculated")
        };
        writeln!(output, "{}", json!({"id":question["id"],"formula":question["formula"],
            "row":question["row"],"column":question["column"],"calc_func_prediction":encoded(value),
            "not_excel_evidence":true,"route":"ordinary Workbook::enter on the whole question sheet"})).unwrap();
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
