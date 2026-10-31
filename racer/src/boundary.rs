//! R7: explicit expm1 precision/store programs and FV=0 rational associations.
//! Own annuity algebra and IEEE 754; current cleared Fpu/elementary chains.
use crate::{Row, dsl::*, eval::*, json_write, ported::excel_numeric::*};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::{
    collections::{HashMap, HashSet},
    fs,
    io::{BufWriter, Write},
    path::Path,
    sync::atomic::{AtomicUsize, Ordering},
    time::Instant,
};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Helper {
    pub id: usize,
    pub t: EvalModel,
    pub exp_store: bool,
    pub difference: EvalModel,
    pub log_store: bool,
    pub association: usize,
    pub arithmetic: [EvalModel; 2],
    pub em_store: bool,
    pub correction: bool,
}

pub fn helpers() -> Vec<Helper> {
    let mut out = vec![];
    for t in MODELS {
        for exp_store in [true, false] {
            for difference in MODELS {
                for correction in [true, false] {
                    for log_store in if correction {
                        vec![true, false]
                    } else {
                        vec![true]
                    } {
                        for association in 0..if correction { 3 } else { 1 } {
                            for a in if correction {
                                MODELS.to_vec()
                            } else {
                                vec![EvalModel::Strict]
                            } {
                                for b in if correction {
                                    MODELS.to_vec()
                                } else {
                                    vec![EvalModel::Strict]
                                } {
                                    for em_store in [true, false] {
                                        out.push(Helper {
                                            id: out.len(),
                                            t,
                                            exp_store,
                                            difference,
                                            log_store,
                                            association,
                                            arithmetic: [a, b],
                                            em_store,
                                            correction,
                                        });
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    assert_eq!(out.len(), 1980);
    out
}

pub fn neg(v: Val) -> Val {
    match v {
        Val::F64(b) => Val::F64(b ^ (1 << 63)),
        Val::Ext(mut e) => {
            e.0[9] ^= 0x80;
            Val::Ext(e)
        }
    }
}
pub fn store<F: Fpu>(v: Val, yes: bool) -> Val {
    if yes { Val::f64(v.to_f64::<F>()) } else { v }
}
fn less_one<F: Fpu>(v: Val) -> bool {
    // Compare the absolute extended value without reducing it to binary64.
    let mut e = v.to_ext();
    e.0[9] &= 0x7f;
    let one = research_load(1.0);
    let d = F::sub(&e, &one, CW);
    d.0[9] & 0x80 != 0 && d.0[..8].iter().any(|b| *b != 0)
}
pub fn helper<F: Fpu>(h: &Helper, rate: f64, n: f64, log: f64) -> Val {
    let _ = rate; // Rate's frozen log1p is an explicit leaf, never refitted.
    let t = neg(binop::<F>(Op::Mul(0, 1), h.t, Val::f64(n), Val::f64(log)));
    let u = store::<F>(Val::Ext(research_exp::<F>(&t.to_ext())), h.exp_store);
    let one = Val::f64(1.0);
    if u.to_ext().0 == one.to_ext().0 {
        return store::<F>(t, h.em_store);
    }
    let d = binop::<F>(Op::Add(0, 1), h.difference, u, neg(one));
    if !h.correction || !less_one::<F>(t) {
        return store::<F>(d, h.em_store);
    }
    let logu = store::<F>(Val::Ext(research_log::<F>(&u.to_ext())), h.log_store);
    let [a, b] = h.arithmetic;
    let e = match h.association {
        0 => binop::<F>(Op::Div(0, 1), b, binop::<F>(Op::Mul(0, 1), a, d, t), logu),
        1 => binop::<F>(Op::Mul(0, 1), b, d, binop::<F>(Op::Div(0, 1), a, t, logu)),
        2 => binop::<F>(Op::Mul(0, 1), b, binop::<F>(Op::Div(0, 1), a, d, logu), t),
        _ => unreachable!(),
    };
    store::<F>(e, h.em_store)
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Expr {
    Leaf(usize),
    Mul(Box<Expr>, Box<Expr>),
    Div(Box<Expr>, Box<Expr>),
}
fn l(i: usize) -> Expr {
    Expr::Leaf(i)
}
fn m(a: Expr, b: Expr) -> Expr {
    Expr::Mul(Box::new(a), Box::new(b))
}
fn d(a: Expr, b: Expr) -> Expr {
    Expr::Div(Box::new(a), Box::new(b))
}
pub fn expressions() -> Vec<Expr> {
    // Leaves P=0,R=1,E=2,T=3. Each appears exactly once; all ten programs
    // are algebraically P*R/(E*T), but floating-point associations stay distinct.
    vec![
        m(d(d(l(0), l(2)), l(3)), l(1)),
        d(m(l(0), l(1)), m(l(2), l(3))),
        m(d(l(0), l(3)), d(l(1), l(2))),
        m(l(0), d(d(l(1), l(2)), l(3))),
        m(l(0), d(l(1), m(l(2), l(3)))),
        d(d(m(l(0), l(1)), l(2)), l(3)),
        m(d(l(0), l(2)), d(l(1), l(3))),
        d(l(0), d(m(l(2), l(3)), l(1))),
        d(l(0), m(l(2), d(l(3), l(1)))),
        d(m(l(0), d(l(1), l(3))), l(2)),
    ]
}
pub fn all_expressions() -> Vec<Expr> {
    fn generate(mask: u8) -> Vec<(Expr, [i8; 4])> {
        if mask.count_ones() == 1 {
            let i = mask.trailing_zeros() as usize;
            let mut powers = [0i8; 4];
            powers[i] = 1;
            return vec![(l(i), powers)];
        }
        let mut out = vec![];
        let mut a = (mask - 1) & mask;
        while a != 0 {
            let b = mask ^ a;
            for (left, lp) in generate(a) {
                for (right, rp) in generate(b) {
                    out.push((
                        m(left.clone(), right.clone()),
                        std::array::from_fn(|i| lp[i] + rp[i]),
                    ));
                    out.push((
                        d(left.clone(), right),
                        std::array::from_fn(|i| lp[i] - rp[i]),
                    ));
                }
            }
            a = (a - 1) & mask;
        }
        out
    }
    let mut out: Vec<_> = generate(15)
        .into_iter()
        .filter_map(|(e, p)| (p == [1, 1, -1, -1]).then_some(e))
        .collect();
    out.sort_by_cached_key(|e| serde_json::to_string(e).unwrap());
    assert_eq!(out.len(), 60);
    for (i, e) in out.iter().enumerate() {
        assert!(!out[..i].contains(e))
    }
    for e in expressions() {
        assert!(out.contains(&e))
    }
    out
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Combine {
    pub id: usize,
    pub expression: usize,
    pub models: [EvalModel; 4],
}
pub fn combines() -> Vec<Combine> {
    combines_for(expressions().len())
}
fn combines_for(count: usize) -> Vec<Combine> {
    let mut out = vec![];
    for expression in 0..count {
        for mut code in 0..81 {
            let mut models = [EvalModel::Strict; 4];
            for v in &mut models {
                *v = MODELS[code % 3];
                code /= 3;
            }
            out.push(Combine {
                id: out.len(),
                expression,
                models,
            });
        }
    }
    out
}
fn eval_expr<F: Fpu>(
    expr: &Expr,
    args: &[Val; 4],
    models: &[EvalModel; 4],
    site: &mut usize,
) -> Val {
    let (a, b, op) = match expr {
        Expr::Leaf(i) => return args[*i],
        Expr::Mul(a, b) => (a, b, Op::Mul(0, 1)),
        Expr::Div(a, b) => (a, b, Op::Div(0, 1)),
    };
    let av = eval_expr::<F>(a, args, models, site);
    let bv = eval_expr::<F>(b, args, models, site);
    let model = models[*site];
    *site += 1;
    binop::<F>(op, model, av, bv)
}
pub fn combine<F: Fpu>(c: &Combine, exprs: &[Expr], r: &Row, em: Val) -> Val {
    let [rate, pv, _, timing, _] = r.args();
    // rate * canonical timing is exact on the scoped finite inputs.
    let tf = binop::<F>(
        Op::Add(0, 1),
        c.models[0],
        Val::f64(1.0),
        Val::f64(rate * timing),
    );
    eval_expr::<F>(
        &exprs[c.expression],
        &[Val::f64(pv), Val::f64(rate), em, tf],
        &c.models,
        &mut 1,
    )
}
fn value_set<F: Fpu>(expr: &Expr, args: &[Vec<Val>; 4]) -> Vec<Val> {
    let (a, b, op) = match expr {
        Expr::Leaf(i) => return args[*i].clone(),
        Expr::Mul(a, b) => (a, b, Op::Mul(0, 1)),
        Expr::Div(a, b) => (a, b, Op::Div(0, 1)),
    };
    let av = value_set::<F>(a, args);
    let bv = value_set::<F>(b, args);
    let mut out = vec![];
    for a in av {
        for b in &bv {
            for model in MODELS {
                let v = binop::<F>(op, model, a, *b);
                if !out.contains(&v) {
                    out.push(v)
                }
            }
        }
    }
    out
}
pub fn roots<F: Fpu>(exprs: &[Expr], r: &Row, ems: &[Val]) -> Vec<Val> {
    let [rate, pv, _, timing, _] = r.args();
    let mut tfs = vec![];
    for model in MODELS {
        let tf = binop::<F>(Op::Add(0, 1), model, Val::f64(1.0), Val::f64(rate * timing));
        if !tfs.contains(&tf) {
            tfs.push(tf)
        }
    }
    let args = [vec![Val::f64(pv)], vec![Val::f64(rate)], ems.to_vec(), tfs];
    let mut out = vec![];
    for e in exprs {
        for v in value_set::<F>(e, &args) {
            if !out.contains(&v) {
                out.push(v)
            }
        }
    }
    out
}
fn group(r: &Row) -> (u64, u64, u64) {
    (r.leaves[0], nper(r).to_bits(), r.log1p)
}
pub fn nper(r: &Row) -> f64 {
    let a = &r.raw["args"][1];
    let v = calc_verdict::decode_argument(a, "args[1]").unwrap();
    match v.core() {
        crate::value::CoreValue::Number(n) => *n,
        _ => panic!("not numeric nper"),
    }
}
fn fv0(r: &Row) -> bool {
    r.args()[2] == 0.0
}

pub fn run<F: Fpu + 'static>(all: &[Row], output: &str) {
    run_inner::<F>(all, output, false)
}
pub fn run_wide<F: Fpu + 'static>(all: &[Row], output: &str) {
    run_inner::<F>(all, output, true)
}
fn run_inner<F: Fpu + 'static>(all: &[Row], output: &str, wide: bool) {
    fs::create_dir_all(output).unwrap();
    let out = Path::new(output);
    let start = Instant::now();
    let rows: Vec<_> = all.iter().filter(|r| fv0(r)).cloned().collect();
    let hs = helpers();
    let es = if wide {
        all_expressions()
    } else {
        expressions()
    };
    let cs = combines_for(es.len());
    json_write(out.join("helpers.json"), &hs);
    json_write(out.join("combines.json"), &cs);
    json_write(out.join("expressions.json"), &es);
    let grammar = hash_bytes(&serde_json::to_vec(&(&hs, &cs, &es)).unwrap());
    let mut index = HashMap::new();
    let mut inputs = vec![];
    let mut indices = vec![];
    for r in &rows {
        let g = group(r);
        let next = inputs.len();
        let id = *index.entry(g).or_insert_with(|| {
            inputs.push(g);
            next
        });
        indices.push(id);
    }
    let worker_limit = std::env::var("DNA_RESEARCH_WORKERS")
        .ok()
        .map(|v| v.parse::<usize>().unwrap())
        .unwrap_or(8);
    assert!(worker_limit >= 1 && worker_limit <= 8);
    let cursor = AtomicUsize::new(0);
    let workers = worker_limit.min(inputs.len());
    let values: Vec<(usize, Vec<Val>, Vec<Val>)> = std::thread::scope(|scope| {
        let mut tasks = vec![];
        for _ in 0..workers {
            let cursor = &cursor;
            let inputs = &inputs;
            let hs = &hs;
            tasks.push(scope.spawn(move || {
                let mut found = vec![];
                loop {
                    let i = cursor.fetch_add(1, Ordering::Relaxed);
                    if i >= inputs.len() {
                        break;
                    }
                    let (r, n, l) = inputs[i];
                    let vals: Vec<_> = hs
                        .iter()
                        .map(|h| {
                            helper::<F>(h, f64::from_bits(r), f64::from_bits(n), f64::from_bits(l))
                        })
                        .collect();
                    let mut distinct = vec![];
                    for &v in &vals {
                        if !distinct.contains(&v) {
                            distinct.push(v)
                        }
                    }
                    found.push((i, vals, distinct));
                }
                found
            }));
        }
        tasks.into_iter().flat_map(|h| h.join().unwrap()).collect()
    });
    let mut values = values;
    values.sort_by_key(|x| x.0);
    println!(
        "prepared {} input groups, {} helpers, {} combines in {:.2}s",
        inputs.len(),
        hs.len(),
        cs.len(),
        start.elapsed().as_secs_f64()
    );
    let mut emfile = BufWriter::new(fs::File::create(out.join("em-values.jsonl")).unwrap());
    for (i, _, distinct) in &values {
        serde_json::to_writer(&mut emfile,&json!({"group":i,"rate":format!("0x{:016x}",inputs[*i].0),"nper":format!("0x{:016x}",inputs[*i].1),"log1p":format!("0x{:016x}",inputs[*i].2),"distinct":distinct.iter().map(|v|v.exact_encoding()).collect::<Vec<_>>()})).unwrap();
        writeln!(emfile).unwrap();
    }
    emfile.flush().unwrap();
    let total = hs.len() * cs.len();
    for (name, without_one) in [("all-fv0", false), ("nonunit-fv0", true)] {
        let mut alive: Vec<_> = (0..total).collect();
        let mut evals = 0usize;
        let mut kill_summary = vec![];
        let mut kills = BufWriter::new(
            fs::File::create(out.join(format!("{name}-first-kills.jsonl"))).unwrap(),
        );
        let mut stats = BufWriter::new(
            fs::File::create(out.join(format!("{name}-witness-kills.jsonl"))).unwrap(),
        );
        let mut orders = 0;
        for (order, r) in rows.iter().enumerate() {
            if without_one && nper(r) == 1.0 {
                continue;
            }
            orders += 1;
            let before = alive.len();
            let mut next = vec![];
            let ems = &values[indices[order]].1;
            for id in alive {
                let h = id / cs.len();
                let c = id % cs.len();
                let got = outcome::<F>(combine::<F>(&cs[c], &es, r, ems[h]));
                evals += 1;
                if got == r.expected() {
                    next.push(id)
                } else {
                    serde_json::to_writer(&mut kills,&json!({"candidate_id":id,"helper":h,"combine":c,"order":order,"trace_index":r.trace_index,"expected":r.expected(),"got":got})).unwrap();
                    writeln!(kills).unwrap();
                }
            }
            let record = json!({"order":order,"trace_index":r.trace_index,"id":r.raw["id"],"killed":before-next.len(),"remaining":next.len()});
            if before != next.len() {
                println!("{name} {record}");
                kill_summary.push(record.clone());
            }
            serde_json::to_writer(&mut stats, &record).unwrap();
            writeln!(stats).unwrap();
            alive = next;
        }
        kills.flush().unwrap();
        stats.flush().unwrap();
        json_write(out.join(format!("{name}-survivors.json")), &alive);
        json_write(
            out.join(format!("{name}-summary.json")),
            &json!({"grammar_hash":grammar,"programs":total,"survivors":alive.len(),"rows":orders,"evaluations":evals,"nonzero_kills":kill_summary,"complete":true}),
        );
    }
    let cursor = AtomicUsize::new(0);
    let progress = AtomicUsize::new(0);
    let counts: Vec<_> = std::thread::scope(|scope| {
        let mut tasks = vec![];
        for worker in 0..workers {
            let cursor = &cursor;
            let progress = &progress;
            let rows = &rows;
            let indices = &indices;
            let values = &values;
            let es = &es;
            tasks.push(scope.spawn(move||{
            let mut f=BufWriter::new(fs::File::create(out.join(format!("reachability-{worker}.jsonl"))).unwrap());
            let mut u=BufWriter::new(fs::File::create(out.join(format!("unreachable-roots-{worker}.jsonl"))).unwrap());
            let mut counts=[0usize;5];
            loop {let order=cursor.fetch_add(1,Ordering::Relaxed);if order>=rows.len(){break}let r=&rows[order];
                let roots=roots::<F>(es,r,&values[indices[order]].2);let published=publish_roots::<F>(&roots);
                let reachable=published.iter().any(|b|bits_outcome(*b)==r.expected());
                counts[0]+=1;counts[1]+=(!reachable)as usize;counts[2]+=published.len();counts[3]+=roots.len();counts[4]+=(!reachable&&nper(r)!=1.0)as usize;
                let record=json!({"order":order,"trace_index":r.trace_index,"id":r.raw["id"],"expected":r.expected(),"reachable":reachable,"outputs":published.iter().map(|b|bits_outcome(*b)).collect::<Vec<_>>()});
                serde_json::to_writer(&mut f,&record).unwrap();writeln!(f).unwrap();
                if !reachable {serde_json::to_writer(&mut u,&json!({"order":order,"trace_index":r.trace_index,"expected":r.expected(),"raw_roots":roots.iter().map(|v|v.exact_encoding()).collect::<Vec<_>>()})).unwrap();writeln!(u).unwrap();}
                let done=progress.fetch_add(1,Ordering::Relaxed)+1;if done%25000==0 {println!("reachability {done}/{} {:.1}s",rows.len(),start.elapsed().as_secs_f64());}
            }counts
        }));
        }
        tasks.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let mut totals = [0usize; 5];
    for c in counts {
        for i in 0..5 {
            totals[i] += c[i]
        }
    }
    assert_eq!(totals[0], rows.len());
    let summary = json!({"grammar_hash":grammar,"programs":total,"helpers":hs.len(),"combines":cs.len(),"scoped_recorded_rows":rows.len(),"input_groups":inputs.len(),"whole_family_unreachable":totals[1],"whole_family_reachable":rows.len()-totals[1],"published_outcomes":totals[2],"tagged_roots":totals[3],"nonunit_unreachable":totals[4],"complete":true,"seconds":start.elapsed().as_secs_f64(),"fresh_validation":0,"backend":std::any::type_name::<F>(),"workers":workers});
    json_write(out.join("summary.json"), &summary);
    println!("{summary}");
}

pub fn self_test(all: &[Row]) {
    test_inner(all, false)
}
pub fn test_wide(all: &[Row]) {
    test_inner(all, true)
}
fn test_inner(all: &[Row], wide: bool) {
    let hs = helpers();
    let es = if wide {
        all_expressions()
    } else {
        expressions()
    };
    let cs = combines_for(es.len());
    let mut tests = 0usize;
    for r in all.iter().filter(|r| fv0(r)).take(8) {
        let n = nper(r);
        let log = f64::from_bits(r.log1p);
        let mut ems = vec![];
        for h in &hs {
            let hard = helper::<Hard>(h, r.args()[0], n, log);
            let soft = helper::<Soft>(h, r.args()[0], n, log);
            assert_eq!(
                hard.exact_encoding(),
                soft.exact_encoding(),
                "helper {} {}",
                h.id,
                r.raw["id"]
            );
            if !ems.contains(&hard) {
                ems.push(hard)
            }
        }
        let mut hard = roots::<Hard>(&es, r, &ems);
        let soft = roots::<Soft>(&es, r, &ems);
        assert_eq!(hard.len(), soft.len());
        for v in &hard {
            assert!(soft.contains(v));
        }
        let mut explicit = vec![];
        for h in &hs {
            let em = helper::<Hard>(h, r.args()[0], n, log);
            for c in &cs {
                let v = combine::<Hard>(c, &es, r, em);
                if !explicit.contains(&v) {
                    explicit.push(v)
                }
                tests += 1;
            }
        }
        assert_eq!(explicit.len(), hard.len());
        for v in hard.drain(..) {
            assert!(explicit.contains(&v));
        }
    }
    println!(
        "R7 self-test: {tests} explicit graph evaluations; all helper/root native-software controls agree"
    );
}

pub fn inverse_controls(all: &[Row], output: &str) {
    let rows: Vec<_> = all.iter().filter(|r| fv0(r)).collect();
    let cs = combines();
    let es = expressions();
    let mut f = BufWriter::new(fs::File::create(output).unwrap());
    let mut controls = 0;
    let mut selected: Vec<_> = (0..8).collect();
    for i in 1..=8 {
        selected.push(i * (rows.len() - 1) / 8)
    }
    for row_index in selected {
        let r = rows[row_index];
        let base = Val::F64(r.leaves[4]);
        let ems = [
            base,
            Val::Ext(base.to_ext()),
            binop::<Hard>(
                Op::Div(0, 1),
                EvalModel::ResidentPc64,
                binop::<Hard>(
                    Op::Add(0, 1),
                    EvalModel::ResidentPc64,
                    base,
                    Val::F64(r.leaves[4] + 1),
                ),
                Val::f64(2.0),
            ),
        ];
        for c in &cs {
            for em in ems {
                let got = outcome::<Hard>(combine::<Hard>(c, &es, r, em));
                let soft = outcome::<Soft>(combine::<Soft>(c, &es, r, em));
                serde_json::to_writer(&mut f,&json!({"row_index":row_index,"combine_id":c.id,"em":em.exact_encoding(),"control_target":got,"software_control_target":soft,"kind":"synthetic forward/inverse arithmetic control, never an Excel answer"})).unwrap();
                writeln!(f).unwrap();
                controls += 1;
            }
        }
    }
    f.flush().unwrap();
    println!("{controls} independent native forward/inverse controls exported");
}
