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

/// Check literal operands and whole calls through the ordinary workbook route.
/// This extends the existing audit; it never supplies replacement Excel answers.
pub fn function_routes(input: &str, out: &str) {
    use std::io::Write;
    let mut output = std::io::BufWriter::new(fs::OpenOptions::new().write(true).create_new(true).open(out).unwrap());
    let rows = crate::lines::<Value>(input);
    for raw in &rows {
        let mut book = Workbook::new_blank();
        let sheet = book.sheets()[0].id;
        let at = |row, column| CellKey::new(sheet, row, column).unwrap();
        let mut echoes = vec![];
        for (index, source) in raw["argument_sources"].as_array().unwrap().iter().enumerate() {
            let argument = calc_verdict::decode_argument(&raw["args"][index], "original argument").unwrap();
            let top = 2 + 40 * index as u32;
            assert!(book.enter(at(top, 1), &format!("={}", source.as_str().unwrap())).accepted());
            if let Some(array) = argument.as_array() {
                for r in 0..array.shape().rows {
                    for c in 0..array.shape().cols {
                        let cell = book.cell(at(top+r as u32, 1+c as u32)).unwrap();
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
                assert_eq!(encoded(value), encoded(&argument), "operand changed: {}", raw["id"]);
                echoes.push(encoded(value));
            }
        }
        let args = raw["args"].as_array().unwrap().iter().map(|a| calc_verdict::decode_argument(a, "original argument").unwrap()).collect();
        let selection = raw["result_index"].as_array().map(|v| (v[0].as_u64().unwrap() as usize, v[1].as_u64().unwrap() as usize));
        let (registered, _) = crate::invoke_selected_control(raw["function"].as_str().unwrap(), args, selection);
        let result = at(1, 30);
        assert!(book.enter(result, raw["formula"].as_str().unwrap()).accepted());
        let cell = book.cell(result).unwrap();
        let ProjectedValue::Calculated { value, .. } = cell.value() else {
            panic!("function call was not calculated") // ubs:ignore — failed route cannot count as a semantic observation
        };
        let workbook = encoded(value);
        assert_eq!(workbook, registered, "workbook/registered discrepancy: {}", raw["id"]);
        assert_eq!(workbook, raw["actual"], "current/frozen product discrepancy: {}", raw["id"]);
        writeln!(output, "{}", json!({"raw":raw,"workbook":workbook,"registered":registered,
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
