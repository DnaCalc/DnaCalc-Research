//! One-period annuity identity, a separate branch hypothesis from the helper.
use crate::{
    Row,
    boundary::{neg, nper},
    dsl::*,
    eval::*,
    json_write,
    ported::excel_numeric::*,
};
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::{
    fs,
    io::{BufWriter, Write},
    path::Path,
};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Identity {
    pub id: usize,
    pub form: usize,
    pub models: [EvalModel; 5],
}
fn programs() -> Vec<Identity> {
    let mut out = vec![];
    for form in 0..3 {
        for mut code in 0..243 {
            let mut models = [EvalModel::Strict; 5];
            for m in &mut models {
                *m = MODELS[code % 3];
                code /= 3;
            }
            out.push(Identity {
                id: out.len(),
                form,
                models,
            });
        }
    }
    out
}
pub fn evaluate<F: Fpu>(p: &Identity, r: f64, pv: f64, fv: f64, timing: f64) -> Val {
    let [a, b, c, d, e] = p.models;
    let [r, pv, fv, one] = [r, pv, fv, 1.0].map(Val::f64);
    let add = |model, x, y| binop::<F>(Op::Add(0, 1), model, x, y);
    let mul = |model, x, y| binop::<F>(Op::Mul(0, 1), model, x, y);
    let div = |model, x, y| binop::<F>(Op::Div(0, 1), model, x, y);
    let value = match p.form {
        0 => {
            let s = add(a, one, r);
            if timing == 1.0 {
                add(e, pv, div(d, fv, s))
            } else {
                add(c, mul(b, pv, s), fv)
            }
        }
        1 => {
            let s = add(a, one, r);
            let n = add(c, mul(b, pv, s), fv);
            let tf = add(d, one, if timing == 1.0 { r } else { Val::f64(0.0) });
            div(e, n, tf)
        }
        2 => {
            let n = add(c, add(b, pv, mul(a, pv, r)), fv);
            let tf = add(d, one, if timing == 1.0 { r } else { Val::f64(0.0) });
            div(e, n, tf)
        }
        _ => unreachable!(),
    };
    neg(value)
}
pub fn run<F: Fpu>(all: &[Row], output: &str) {
    fs::create_dir_all(output).unwrap();
    let out = Path::new(output);
    let ps = programs();
    let rows: Vec<_> = all.iter().filter(|r| nper(r) == 1.0).collect();
    let mut alive: Vec<_> = (0..ps.len()).collect();
    json_write(out.join("programs.json"), &ps);
    let mut f = BufWriter::new(fs::File::create(out.join("first-kills.jsonl")).unwrap());
    let mut stats = vec![];
    let mut evaluations = 0;
    for (order, r) in rows.iter().enumerate() {
        let before = alive.len();
        let mut next = vec![];
        let [rate, pv, fv, timing, _] = r.args();
        for id in alive {
            let got = outcome::<F>(evaluate::<F>(&ps[id], rate, pv, fv, timing));
            evaluations += 1;
            if got == r.expected() {
                next.push(id)
            } else {
                serde_json::to_writer(&mut f,&json!({"candidate_id":id,"order":order,"trace_index":r.trace_index,"witness":r.raw,"expected":r.expected(),"got":got})).unwrap();
                writeln!(f).unwrap();
            }
        }
        if before != next.len() {
            let stat = json!({"order":order,"trace_index":r.trace_index,"id":r.raw["id"],"killed":before-next.len(),"remaining":next.len()});
            println!("identity {stat}");
            stats.push(stat)
        }
        alive = next;
    }
    f.flush().unwrap();
    let mut unreachable = 0;
    let mut union = BufWriter::new(fs::File::create(out.join("reachability.jsonl")).unwrap());
    for (order, r) in rows.iter().enumerate() {
        let [rate, pv, fv, timing, _] = r.args();
        let mut roots = vec![];
        for p in &ps {
            let v = evaluate::<F>(p, rate, pv, fv, timing);
            if !roots.contains(&v) {
                roots.push(v)
            }
        }
        let outputs = publish_roots::<F>(&roots);
        let reachable = outputs.iter().any(|b| bits_outcome(*b) == r.expected());
        unreachable += (!reachable) as usize;
        serde_json::to_writer(&mut union,&json!({"order":order,"trace_index":r.trace_index,"id":r.raw["id"],"expected":r.expected(),"reachable":reachable,"outputs":outputs.iter().map(|b|bits_outcome(*b)).collect::<Vec<_>>()})).unwrap();
        writeln!(union).unwrap();
    }
    let survivors: Vec<_> = alive.iter().map(|i| &ps[*i]).collect();
    json_write(out.join("survivors.json"), &survivors);
    let summary = json!({"programs":ps.len(),"grammar_hash":hash_bytes(&serde_json::to_vec(&ps).unwrap()),"rows":rows.len(),"survivors":survivors.len(),"evaluations":evaluations,"kills":stats,"whole_family_unreachable":unreachable,"whole_family_reachable":rows.len()-unreachable,"complete":true,"backend":std::any::type_name::<F>()});
    json_write(out.join("summary.json"), &summary);
    println!("{summary}");
}
