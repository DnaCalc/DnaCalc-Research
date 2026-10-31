//! Count-one Poisson normalization identities, with independently derived constants.
use crate::poisson::{outcome, publish_roots};
use crate::{
    boundary::{neg, store},
    dsl::*,
    eval::*,
    json_write,
    poisson::{self, MassRow},
    ported::excel_numeric::*,
};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::{
    collections::HashMap,
    fs,
    io::{BufWriter, Write},
    path::Path,
    sync::atomic::{AtomicUsize, Ordering},
    time::Instant,
};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Node {
    Mean,
    One,
    Number {
        bits: u64,
    },
    Constant {
        index: usize,
        extended: bool,
    },
    Neg(Box<Node>),
    ReciprocalLog {
        provider: usize,
        reciprocal_model: EvalModel,
    },
    DerivedG {
        tau_extended: bool,
        model: EvalModel,
    },
    Add {
        site: usize,
        a: Box<Node>,
        b: Box<Node>,
    },
    Mul {
        site: usize,
        a: Box<Node>,
        b: Box<Node>,
    },
    Div {
        site: usize,
        a: Box<Node>,
        b: Box<Node>,
    },
    Log {
        provider: usize,
    },
    Exp {
        provider: usize,
        entry_store: bool,
        input: Box<Node>,
    },
    Power {
        provider: usize,
        exponent: i32,
    },
    ScaledUnderflow {
        scaled: Box<Node>,
        denominator: Box<Node>,
        fallback: Box<Node>,
        site: usize,
        reciprocal_site: Option<usize>,
    },
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Template {
    pub id: usize,
    pub deviance: usize,
    pub normalization: usize,
    pub log_provider: usize,
    pub exp_provider: usize,
    pub exp_entry_store: bool,
    pub s_extended: bool,
    pub other_extended: bool,
    pub graph: Node,
}
fn add(site: usize, a: Node, b: Node) -> Node {
    Node::Add {
        site,
        a: Box::new(a),
        b: Box::new(b),
    }
}
fn mul(site: usize, a: Node, b: Node) -> Node {
    Node::Mul {
        site,
        a: Box::new(a),
        b: Box::new(b),
    }
}
fn div(site: usize, a: Node, b: Node) -> Node {
    Node::Div {
        site,
        a: Box::new(a),
        b: Box::new(b),
    }
}
fn minus(a: Node) -> Node {
    Node::Neg(Box::new(a))
}
fn sub(site: usize, a: Node, b: Node) -> Node {
    add(site, a, minus(b))
}
fn count_sites(node: &Node, counts: &mut [usize; 8]) {
    match node {
        Node::Add { site, a, b } | Node::Mul { site, a, b } | Node::Div { site, a, b } => {
            counts[*site] += 1;
            count_sites(a, counts);
            count_sites(b, counts)
        }
        Node::Neg(a) | Node::Exp { input: a, .. } => count_sites(a, counts),
        Node::ScaledUnderflow {
            scaled,
            denominator,
            fallback,
            site,
            reciprocal_site,
        } => {
            counts[*site] += 1;
            if let Some(r) = reciprocal_site {
                counts[*r] += 1
            }
            count_sites(scaled, counts);
            count_sites(denominator, counts);
            count_sites(fallback, counts);
        }
        _ => {}
    }
}
pub fn templates() -> Vec<Template> {
    let mut out = vec![];
    for deviance in 0..3 {
        for normalization in 0..5 {
            for log_provider in 0..4 {
                for exp_provider in 0..5 {
                    for exp_entry_store in if exp_provider < 2 {
                        vec![true, false]
                    } else {
                        vec![false]
                    } {
                        for s_extended in [false, true] {
                            for other_extended in [false, true] {
                                let log = Node::Log {
                                    provider: log_provider,
                                };
                                let d = match deviance {
                                    0 => sub(1, sub(0, Node::Mean, Node::One), log),
                                    1 => sub(1, Node::Mean, add(0, Node::One, log)),
                                    2 => sub(1, sub(0, Node::Mean, log), Node::One),
                                    _ => unreachable!(),
                                };
                                let s = Node::Constant {
                                    index: 1,
                                    extended: s_extended,
                                };
                                let other = Node::Constant {
                                    index: if normalization == 1 || normalization == 4 {
                                        0
                                    } else {
                                        2
                                    },
                                    extended: other_extended,
                                };
                                let exp = |input| Node::Exp {
                                    provider: exp_provider,
                                    entry_store: exp_entry_store,
                                    input: Box::new(input),
                                };
                                let graph = match normalization {
                                    0 => div(3, exp(minus(add(2, s, d))), other),
                                    1 => exp(sub(3, minus(add(2, s, d)), other)),
                                    2 => div(3, mul(2, exp(minus(d)), exp(minus(s))), other),
                                    3 => div(3, exp(minus(d)), mul(2, other, exp(s))),
                                    4 => mul(3, exp(minus(d)), exp(minus(add(2, s, other)))),
                                    _ => unreachable!(),
                                };
                                let mut sites = [0; 8];
                                count_sites(&graph, &mut sites);
                                assert_eq!(sites, [1, 1, 1, 1, 0, 0, 0, 0]);
                                out.push(Template {
                                    id: out.len(),
                                    deviance,
                                    normalization,
                                    log_provider,
                                    exp_provider,
                                    exp_entry_store,
                                    s_extended,
                                    other_extended,
                                    graph,
                                });
                            }
                        }
                    }
                }
            }
        }
    }
    assert_eq!(out.len(), 1680);
    let mut hashes = HashMap::<String, Vec<&Node>>::new();
    for t in &out {
        let hash = hash_bytes(&serde_json::to_vec(&t.graph).unwrap());
        let bucket = hashes.entry(hash).or_default();
        assert!(!bucket.contains(&&t.graph));
        bucket.push(&t.graph);
    }
    out
}
pub fn expanded_templates() -> Vec<Template> {
    let mut out = templates();
    for variant in 0..10 {
        for log_provider in 0..4 {
            for exp_provider in 0..5 {
                for entry in if exp_provider < 2 {
                    vec![true, false]
                } else {
                    vec![false]
                } {
                    let l = Node::Log {
                        provider: log_provider,
                    };
                    let two = || Node::Number {
                        bits: 2.0f64.to_bits(),
                    };
                    let exp = |input| Node::Exp {
                        provider: exp_provider,
                        entry_store: entry,
                        input: Box::new(input),
                    };
                    let scaled = || exp(sub(1, mul(0, two(), l.clone()), Node::Mean));
                    let graph = match variant {
                        0 => exp(sub(0, l, Node::Mean)),
                        1 => mul(0, exp(l), exp(minus(Node::Mean))),
                        2 => exp(sub(2, sub(1, mul(0, two(), l.clone()), Node::Mean), l)),
                        3 => exp(sub(2, mul(0, two(), l.clone()), add(1, Node::Mean, l))),
                        4 => exp(sub(2, sub(1, mul(0, two(), l.clone()), l), Node::Mean)),
                        5 => div(2, scaled(), Node::Mean),
                        6 => mul(3, scaled(), div(2, Node::One, Node::Mean)),
                        7 => mul(
                            2,
                            div(1, exp(mul(0, two(), l)), Node::Mean),
                            exp(minus(Node::Mean)),
                        ),
                        8 => Node::ScaledUnderflow {
                            scaled: Box::new(scaled()),
                            denominator: Box::new(Node::Mean),
                            fallback: Box::new(exp(sub(3, l, Node::Mean))),
                            site: 2,
                            reciprocal_site: None,
                        },
                        9 => Node::ScaledUnderflow {
                            scaled: Box::new(scaled()),
                            denominator: Box::new(Node::Mean),
                            fallback: Box::new(exp(sub(4, l, Node::Mean))),
                            site: 3,
                            reciprocal_site: Some(2),
                        },
                        _ => unreachable!(),
                    };
                    out.push(Template {
                        id: out.len(),
                        deviance: 0,
                        normalization: 5 + variant,
                        log_provider,
                        exp_provider,
                        exp_entry_store: entry,
                        s_extended: false,
                        other_extended: false,
                        graph,
                    });
                }
            }
        }
    }
    for variant in 0..4 {
        for power_provider in 0..4 {
            for exp_provider in 0..5 {
                for entry in if exp_provider < 2 {
                    vec![true, false]
                } else {
                    vec![false]
                } {
                    let exp = Node::Exp {
                        provider: exp_provider,
                        entry_store: entry,
                        input: Box::new(minus(Node::Mean)),
                    };
                    let pow = Node::Power {
                        provider: power_provider,
                        exponent: match variant {
                            0 => 1,
                            1 => -1,
                            _ => 2,
                        },
                    };
                    let graph = match variant {
                        0 => mul(0, pow, exp),
                        1 => div(0, exp, pow),
                        2 => mul(1, exp, div(0, pow, Node::Mean)),
                        3 => div(1, mul(0, exp, pow), Node::Mean),
                        _ => unreachable!(),
                    };
                    out.push(Template {
                        id: out.len(),
                        deviance: power_provider,
                        normalization: 15 + variant,
                        log_provider: 0,
                        exp_provider,
                        exp_entry_store: entry,
                        s_extended: false,
                        other_extended: false,
                        graph,
                    });
                }
            }
        }
    }
    let mut hashes = HashMap::<String, Vec<&Node>>::new();
    for t in &out {
        let hash = hash_bytes(&serde_json::to_vec(&t.graph).unwrap());
        let bucket = hashes.entry(hash).or_default();
        assert!(!bucket.contains(&&t.graph), "duplicate graph {}", t.id);
        bucket.push(&t.graph);
        variants(t);
    }
    assert_eq!(out.len(), 2072);
    assert_eq!(out.iter().map(variants).sum::<usize>(), 152040);
    out
}
fn models(mut code: usize) -> [EvalModel; 8] {
    std::array::from_fn(|_| {
        let m = MODELS[code % 3];
        code /= 3;
        m
    })
}
fn models_in(mut code: usize, choices: &[EvalModel]) -> [EvalModel; 8] {
    std::array::from_fn(|_| {
        let m = choices[code % choices.len()];
        code /= choices.len();
        m
    })
}
fn model_index(m: EvalModel) -> usize {
    match m {
        EvalModel::Strict => 0,
        EvalModel::StoredPc64 => 1,
        EvalModel::ResidentPc64 => 2,
        EvalModel::StoredPc53 => 3,
        EvalModel::ResidentPc53 => 4,
    }
}
fn rewrite(node: &Node, mapper: &impl Fn(&Node) -> Option<Node>) -> Node {
    if let Some(n) = mapper(node) {
        return n;
    }
    match node {
        Node::Neg(a) => minus(rewrite(a, mapper)),
        Node::Add { site, a, b } => add(*site, rewrite(a, mapper), rewrite(b, mapper)),
        Node::Mul { site, a, b } => mul(*site, rewrite(a, mapper), rewrite(b, mapper)),
        Node::Div { site, a, b } => div(*site, rewrite(a, mapper), rewrite(b, mapper)),
        Node::Exp {
            provider,
            entry_store,
            input,
        } => Node::Exp {
            provider: *provider,
            entry_store: *entry_store,
            input: Box::new(rewrite(input, mapper)),
        },
        Node::ScaledUnderflow {
            scaled,
            denominator,
            fallback,
            site,
            reciprocal_site,
        } => Node::ScaledUnderflow {
            scaled: Box::new(rewrite(scaled, mapper)),
            denominator: Box::new(rewrite(denominator, mapper)),
            fallback: Box::new(rewrite(fallback, mapper)),
            site: *site,
            reciprocal_site: *reciprocal_site,
        },
        _ => node.clone(),
    }
}
fn has(node: &Node, p: &impl Fn(&Node) -> bool) -> bool {
    if p(node) {
        return true;
    }
    match node {
        Node::Neg(a) | Node::Exp { input: a, .. } => has(a, p),
        Node::Add { a, b, .. } | Node::Mul { a, b, .. } | Node::Div { a, b, .. } => {
            has(a, p) || has(b, p)
        }
        Node::ScaledUnderflow {
            scaled,
            denominator,
            fallback,
            ..
        } => has(scaled, p) || has(denominator, p) || has(fallback, p),
        _ => false,
    }
}
pub fn refined_templates() -> Vec<Template> {
    make_refined(false)
}
pub fn full_templates() -> Vec<Template> {
    make_refined(true)
}
fn make_refined(full: bool) -> Vec<Template> {
    let choices: &[EvalModel] = if full { &WIDE_MODELS } else { &MODELS };
    let base = expanded_templates();
    let mut out = base.clone();
    for t in base
        .iter()
        .filter(|t| has(&t.graph, &|n| matches!(n, Node::Log { .. })))
    {
        for &reciprocal_model in choices {
            let mut c = t.clone();
            c.id = out.len();
            c.graph = rewrite(&t.graph, &|n| {
                if let Node::Log { provider } = n {
                    Some(minus(Node::ReciprocalLog {
                        provider: *provider,
                        reciprocal_model,
                    }))
                } else {
                    None
                }
            });
            out.push(c);
        }
    }
    assert_eq!(out.len(), if full { 11872 } else { 7952 });
    let mut hashes = HashMap::<String, Vec<usize>>::new();
    for (i, t) in out.iter().enumerate() {
        hashes
            .entry(hash_bytes(&serde_json::to_vec(&t.graph).unwrap()))
            .or_default()
            .push(i);
    }
    let original = out.clone();
    for t in original
        .iter()
        .filter(|t| has(&t.graph, &|n| matches!(n, Node::Constant { index: 2, .. })))
    {
        for tau_extended in [false, true] {
            for &model in choices {
                let mut c = t.clone();
                c.graph = rewrite(&t.graph, &|n| {
                    if matches!(n, Node::Constant { index: 2, .. }) {
                        Some(Node::DerivedG {
                            tau_extended,
                            model,
                        })
                    } else {
                        None
                    }
                });
                let hash = hash_bytes(&serde_json::to_vec(&c.graph).unwrap());
                let bucket = hashes.entry(hash).or_default();
                if !bucket.iter().any(|i| out[*i].graph == c.graph) {
                    c.id = out.len();
                    bucket.push(c.id);
                    out.push(c);
                }
            }
        }
    }
    assert_eq!(out.len(), if full { 42112 } else { 20048 });
    assert_eq!(
        out.iter()
            .map(|t| variants_in(t, choices.len()))
            .sum::<usize>(),
        if full { 26043360 } else { 1585920 }
    );
    for t in &out {
        variants_in(t, choices.len());
    }
    out
}
fn variants(t: &Template) -> usize {
    variants_in(t, 3)
}
fn variants_in(t: &Template, base: usize) -> usize {
    let mut sites = [0; 8];
    count_sites(&t.graph, &mut sites);
    let arity = sites.iter().filter(|c| **c != 0).count();
    assert!(sites[..arity].iter().all(|c| *c == 1));
    assert!(sites[arity..].iter().all(|c| *c == 0));
    base.pow(arity as u32)
}
fn is_zero(v: Val) -> bool {
    let e = v.to_ext();
    e.0[..8].iter().all(|b| *b == 0) && (e.0[8] == 0) && (e.0[9] & 0x7f == 0)
}
fn power<F: Fpu>(mean: u64, provider: usize, exponent: i32) -> Val {
    let x = f64::from_bits(mean);
    let e = exponent as f64;
    match provider {
        0 | 1 => store::<F>(Val::Ext(research_pow::<F>(x, e)), provider == 1),
        2 => Val::f64(crt_pow(x, e)),
        3 => Val::f64(libm::pow(x, e)),
        _ => unreachable!(),
    }
}
type Constants = [[Val; 2]; 9];
pub fn constants<F: Fpu>(path: &str) -> Constants {
    let data: Value = serde_json::from_slice(&fs::read(path).unwrap()).unwrap();
    let mut result = [[Val::f64(0.0); 2]; 9];
    for (idx, name) in ["C", "S", "G", "TAU"].iter().enumerate() {
        let c = &data["constants"][*name];
        if c.is_null() {
            assert_eq!(idx, 3);
            continue;
        }
        let b = u64::from_str_radix(c["binary64"].as_str().unwrap(), 16).unwrap();
        let raw = c["ext80_le"].as_str().unwrap();
        let mut ext = [0u8; 10];
        for i in 0..10 {
            ext[i] = u8::from_str_radix(&raw[2 * i..2 * i + 2], 16).unwrap()
        }
        result[idx] = [Val::F64(b), Val::Ext(Extended(ext))];
    }
    if !data["constants"]["TAU"].is_null() {
        for t in 0..2 {
            for m in 0..5 {
                let tau = result[3][t];
                let v = if m == 0 {
                    Val::f64(tau.to_f64::<F>().sqrt())
                } else {
                    let cw = if m >= 3 { (CW & !0x300) | 0x200 } else { CW };
                    store::<F>(Val::Ext(F::sqrt(&tau.to_ext(), cw)), m == 1 || m == 3)
                };
                result[4 + m][t] = v;
            }
        }
    }
    result
}
// The research evaluator memoizes only pure requests with complete tagged
// operands, operation/provider and precision model, within one input/backend.
// Physical program descriptions remain separate; no finite prediction-vector
// equality is used as a global program equivalence. This grammar has no FP
// status-word observers, unmasked traps or other side effects.
struct RowCache<F: Fpu> {
    mean: u64,
    enabled: bool,
    log: [Option<Val>; 4],
    reciprocal: [[Option<Val>; 5]; 4],
    powers: [[Option<Val>; 3]; 4],
    arithmetic: HashMap<(u8, EvalModel, Val, Val), Val>,
    exponentials: HashMap<(usize, Val), Val>,
    backend: std::marker::PhantomData<F>,
}
impl<F: Fpu> RowCache<F> {
    fn new(mean: u64, enabled: bool) -> Self {
        Self {
            mean,
            enabled,
            log: [None; 4],
            reciprocal: [[None; 5]; 4],
            powers: [[None; 3]; 4],
            arithmetic: HashMap::new(),
            exponentials: HashMap::new(),
            backend: std::marker::PhantomData,
        }
    }
    fn logarithm(&mut self, p: usize) -> Val {
        if self.enabled {
            if let Some(v) = self.log[p] {
                return v;
            }
        }
        let v = poisson::log::<F>(f64::from_bits(self.mean), p);
        if self.enabled {
            self.log[p] = Some(v)
        }
        v
    }
    fn reciprocal_log(&mut self, p: usize, m: EvalModel) -> Val {
        let i = model_index(m);
        if self.enabled {
            if let Some(v) = self.reciprocal[p][i] {
                return v;
            }
        }
        let v = poisson::log_value::<F>(
            binop::<F>(Op::Div(0, 1), m, Val::f64(1.0), Val::F64(self.mean)),
            p,
        );
        if self.enabled {
            self.reciprocal[p][i] = Some(v)
        }
        v
    }
    fn power(&mut self, p: usize, e: i32) -> Val {
        let i = match e {
            -1 => 0,
            1 => 1,
            2 => 2,
            _ => unreachable!(),
        };
        if self.enabled {
            if let Some(v) = self.powers[p][i] {
                return v;
            }
        }
        let v = power::<F>(self.mean, p, e);
        if self.enabled {
            self.powers[p][i] = Some(v)
        }
        v
    }
    fn exp(&mut self, p: usize, x: Val) -> Val {
        let k = (p, x);
        if self.enabled {
            if let Some(v) = self.exponentials.get(&k) {
                return *v;
            }
        }
        let v = poisson::exp::<F>(x, p);
        if self.enabled {
            self.exponentials.insert(k, v);
        }
        v
    }
    fn op(&mut self, op: Op, m: EvalModel, a: Val, b: Val) -> Val {
        let kind = match op {
            Op::Add(..) => 0,
            Op::Mul(..) => 1,
            Op::Div(..) => 2,
            _ => unreachable!(),
        };
        let k = (kind, m, a, b);
        if self.enabled {
            if let Some(v) = self.arithmetic.get(&k) {
                return *v;
            }
        }
        let v = binop::<F>(op, m, a, b);
        if self.enabled {
            self.arithmetic.insert(k, v);
        }
        v
    }
}
fn evaluate<F: Fpu>(
    node: &Node,
    mean: u64,
    constants: &Constants,
    models: &[EvalModel; 8],
    cache: &mut RowCache<F>,
) -> Val {
    match node {
        Node::Mean => Val::F64(mean),
        Node::One => Val::f64(1.0),
        Node::Number { bits } => Val::F64(*bits),
        Node::Constant { index, extended } => constants[*index][*extended as usize],
        Node::Power { provider, exponent } => cache.power(*provider, *exponent),
        Node::ReciprocalLog {
            provider,
            reciprocal_model,
        } => cache.reciprocal_log(*provider, *reciprocal_model),
        Node::DerivedG {
            tau_extended,
            model,
        } => {
            assert!(
                !is_zero(constants[3][0]),
                "derived sqrt(2*pi) requires a TAU constant certificate"
            );
            constants[4 + model_index(*model)][*tau_extended as usize]
        }
        Node::ScaledUnderflow {
            scaled,
            denominator,
            fallback,
            site,
            reciprocal_site,
        } => {
            let s = evaluate::<F>(scaled, mean, constants, models, cache);
            if is_zero(s) {
                evaluate::<F>(fallback, mean, constants, models, cache)
            } else {
                let d = evaluate::<F>(denominator, mean, constants, models, cache);
                if let Some(r) = reciprocal_site {
                    let inv = cache.op(Op::Div(0, 1), models[*r], Val::f64(1.0), d);
                    cache.op(Op::Mul(0, 1), models[*site], s, inv)
                } else {
                    cache.op(Op::Div(0, 1), models[*site], s, d)
                }
            }
        }
        Node::Neg(a) => neg(evaluate::<F>(a, mean, constants, models, cache)),
        Node::Log { provider } => cache.logarithm(*provider),
        Node::Exp {
            provider,
            entry_store,
            input,
        } => {
            let x = store::<F>(
                evaluate::<F>(input, mean, constants, models, cache),
                *entry_store,
            );
            cache.exp(*provider, x)
        }
        Node::Add { site, a, b } | Node::Mul { site, a, b } | Node::Div { site, a, b } => {
            let op = match node {
                Node::Add { .. } => Op::Add(0, 1),
                Node::Mul { .. } => Op::Mul(0, 1),
                Node::Div { .. } => Op::Div(0, 1),
                _ => unreachable!(),
            };
            let av = evaluate::<F>(a, mean, constants, models, cache);
            let bv = evaluate::<F>(b, mean, constants, models, cache);
            cache.op(op, models[*site], av, bv)
        }
    }
}
fn roots<F: Fpu>(
    node: &Node,
    mean: u64,
    constants: &Constants,
    choices: &[EvalModel],
    cache: &mut RowCache<F>,
) -> Vec<Val> {
    match node {
        Node::Mean => vec![Val::F64(mean)],
        Node::One => vec![Val::f64(1.0)],
        Node::Number { bits } => vec![Val::F64(*bits)],
        Node::Constant { index, extended } => vec![constants[*index][*extended as usize]],
        Node::Power { provider, exponent } => vec![cache.power(*provider, *exponent)],
        Node::ReciprocalLog {
            provider,
            reciprocal_model,
        } => vec![cache.reciprocal_log(*provider, *reciprocal_model)],
        Node::DerivedG {
            tau_extended,
            model,
        } => {
            assert!(
                !is_zero(constants[3][0]),
                "derived sqrt(2*pi) requires a TAU constant certificate"
            );
            vec![constants[4 + model_index(*model)][*tau_extended as usize]]
        }
        Node::ScaledUnderflow {
            scaled,
            denominator,
            fallback,
            reciprocal_site,
            ..
        } => {
            let mut out = vec![];
            for s in roots::<F>(scaled, mean, constants, choices, cache) {
                if is_zero(s) {
                    for f in roots::<F>(fallback, mean, constants, choices, cache) {
                        if !out.contains(&f) {
                            out.push(f)
                        }
                    }
                } else {
                    for d in roots::<F>(denominator, mean, constants, choices, cache) {
                        for &m in choices {
                            if reciprocal_site.is_some() {
                                for &q in choices {
                                    let inv = cache.op(Op::Div(0, 1), q, Val::f64(1.0), d);
                                    let v = cache.op(Op::Mul(0, 1), m, s, inv);
                                    if !out.contains(&v) {
                                        out.push(v)
                                    }
                                }
                            } else {
                                let v = cache.op(Op::Div(0, 1), m, s, d);
                                if !out.contains(&v) {
                                    out.push(v)
                                }
                            }
                        }
                    }
                }
            }
            out
        }
        Node::Neg(a) => roots::<F>(a, mean, constants, choices, cache)
            .into_iter()
            .map(neg)
            .collect(),
        Node::Log { provider } => vec![cache.logarithm(*provider)],
        Node::Exp {
            provider,
            entry_store,
            input,
        } => {
            let mut out = vec![];
            for v in roots::<F>(input, mean, constants, choices, cache) {
                let x = cache.exp(*provider, store::<F>(v, *entry_store));
                if !out.contains(&x) {
                    out.push(x)
                }
            }
            out
        }
        Node::Add { a, b, .. } | Node::Mul { a, b, .. } | Node::Div { a, b, .. } => {
            let op = match node {
                Node::Add { .. } => Op::Add(0, 1),
                Node::Mul { .. } => Op::Mul(0, 1),
                Node::Div { .. } => Op::Div(0, 1),
                _ => unreachable!(),
            };
            let av = roots::<F>(a, mean, constants, choices, cache);
            let bv = roots::<F>(b, mean, constants, choices, cache);
            let mut out = vec![];
            for a in av {
                for b in &bv {
                    for &model in choices {
                        let x = cache.op(op, model, a, *b);
                        if !out.contains(&x) {
                            out.push(x)
                        }
                    }
                }
            }
            out
        }
    }
}
pub fn run<F: Fpu + 'static>(all: &[MassRow], const_path: &str, output: &str) {
    run_inner::<F>(all, const_path, output, 0)
}
pub fn run_expanded<F: Fpu + 'static>(all: &[MassRow], const_path: &str, output: &str) {
    run_inner::<F>(all, const_path, output, 1)
}
pub fn run_refined<F: Fpu + 'static>(all: &[MassRow], const_path: &str, output: &str) {
    run_inner::<F>(all, const_path, output, 2)
}
pub fn run_full<F: Fpu + 'static>(all: &[MassRow], const_path: &str, output: &str) {
    run_inner::<F>(all, const_path, output, 3)
}
fn run_inner<F: Fpu + 'static>(all: &[MassRow], const_path: &str, output: &str, level: usize) {
    fs::create_dir_all(output).unwrap();
    let out = Path::new(output);
    let ts = match level {
        0 => templates(),
        1 => expanded_templates(),
        2 => refined_templates(),
        _ => full_templates(),
    };
    let constants = constants::<F>(const_path);
    let start = Instant::now();
    let choices: &[EvalModel] = if level >= 3 { &WIDE_MODELS } else { &MODELS };
    let rows: Vec<_> = all.iter().filter(|r| r.count == 1).collect();
    let mut programs = vec![];
    let mut layout = vec![];
    for t in &ts {
        let count = variants_in(t, choices.len());
        layout.push(json!({"template":t.id,"first_candidate_id":programs.len(),"variants":count}));
        for code in 0..count {
            programs.push((t.id, code));
        }
    }
    let total = programs.len();
    let grammar = hash_bytes(&serde_json::to_vec(&ts).unwrap());
    json_write(out.join("program-layout.json"), &layout);
    json_write(out.join("templates.json"), &ts);
    let mut alive: Vec<_> = (0..total).collect();
    let mut evaluations = 0;
    let mut kill_summary = vec![];
    let mut kills = BufWriter::new(fs::File::create(out.join("first-kills.jsonl")).unwrap());
    let mut stats = BufWriter::new(fs::File::create(out.join("witness-kills.jsonl")).unwrap());
    for (order, r) in rows.iter().enumerate() {
        let before = alive.len();
        let mut next = vec![];
        let mut cache = RowCache::<F>::new(r.mean, true);
        for id in alive {
            let (template_id, model_code) = programs[id];
            let got = outcome::<F>(evaluate::<F>(
                &ts[template_id].graph,
                r.mean,
                &constants,
                &models_in(model_code, choices),
                &mut cache,
            ));
            evaluations += 1;
            if got == r.raw["expected"].as_str().unwrap() {
                next.push(id)
            } else {
                serde_json::to_writer(&mut kills,&json!({"candidate_id":id,"template_id":template_id,"model_code":model_code,"order":order,"trace_index":r.trace_index,"expected":r.raw["expected"],"got":got})).unwrap();
                writeln!(kills).unwrap();
            }
        }
        let record = json!({"order":order,"trace_index":r.trace_index,"id":r.raw["id"],"killed":before-next.len(),"remaining":next.len()});
        if before != next.len() {
            println!("count-one {record}");
            kill_summary.push(record.clone());
        }
        serde_json::to_writer(&mut stats, &record).unwrap();
        writeln!(stats).unwrap();
        alive = next;
    }
    kills.flush().unwrap();
    stats.flush().unwrap();
    json_write(out.join("survivors.json"), &alive);
    let cursor = AtomicUsize::new(0);
    let progress = AtomicUsize::new(0);
    let workers = std::env::var("DNA_RESEARCH_WORKERS")
        .ok()
        .map(|v| v.parse::<usize>().unwrap())
        .unwrap_or(4);
    assert!(workers >= 1 && workers <= 8);
    let counts: Vec<_> = std::thread::scope(|scope| {
        let mut tasks = vec![];
        for worker in 0..workers {
            let cursor = &cursor;
            let progress = &progress;
            let rows = &rows;
            let ts = &ts;
            let constants = &constants;
            tasks.push(scope.spawn(move||{
            let mut f=BufWriter::new(fs::File::create(out.join(format!("reachability-{worker}.jsonl"))).unwrap());let mut u=BufWriter::new(fs::File::create(out.join(format!("unreachable-roots-{worker}.jsonl"))).unwrap());let mut counts=[0usize;4];
            loop {let order=cursor.fetch_add(1,Ordering::Relaxed);if order>=rows.len(){break}let r=rows[order];let mut values=vec![];let mut cache=RowCache::<F>::new(r.mean,true);
                for t in ts {for v in roots::<F>(&t.graph,r.mean,constants,choices,&mut cache) {if !values.contains(&v){values.push(v)}}}
                let outputs=publish_roots::<F>(&values);let reachable=outputs.iter().any(|b|bits_outcome(*b)==r.raw["expected"].as_str().unwrap());counts[0]+=1;counts[1]+=(!reachable)as usize;counts[2]+=outputs.len();counts[3]+=values.len();
                serde_json::to_writer(&mut f,&json!({"order":order,"trace_index":r.trace_index,"id":r.raw["id"],"expected":r.raw["expected"],"reachable":reachable,"outputs":outputs.iter().map(|b|bits_outcome(*b)).collect::<Vec<_>>()})).unwrap();writeln!(f).unwrap();
                if !reachable {serde_json::to_writer(&mut u,&json!({"order":order,"trace_index":r.trace_index,"expected":r.raw["expected"],"raw_roots":values.iter().map(|v|v.exact_encoding()).collect::<Vec<_>>()})).unwrap();writeln!(u).unwrap();}
                let done=progress.fetch_add(1,Ordering::Relaxed)+1;if done%500==0 {println!("count-one reachability {done}/{} {:.1}s",rows.len(),start.elapsed().as_secs_f64());}
            }counts
        }));
        }
        tasks.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let mut totals = [0usize; 4];
    for c in counts {
        for i in 0..4 {
            totals[i] += c[i]
        }
    }
    assert_eq!(totals[0], rows.len());
    let summary = json!({"grammar_hash":grammar,"programs":total,"templates":ts.len(),"models_per_template":if level>0{json!(null)}else{json!(81)},"rows":rows.len(),"survivors":alive.len(),"evaluations":evaluations,"kills":kill_summary,"whole_family_unreachable":totals[1],"whole_family_reachable":rows.len()-totals[1],"published_outcomes":totals[2],"tagged_roots":totals[3],"complete":true,"seconds":start.elapsed().as_secs_f64(),"backend":std::any::type_name::<F>(),"workers":workers,"fresh_validation":0});
    json_write(out.join("summary.json"), &summary);
    println!("{summary}");
}
pub fn self_test(all: &[MassRow], const_path: &str) {
    test_inner(all, const_path, 0)
}
pub fn test_expanded(all: &[MassRow], const_path: &str) {
    test_inner(all, const_path, 1)
}
pub fn test_refined(all: &[MassRow], const_path: &str) {
    test_inner(all, const_path, 2)
}
pub fn test_full(all: &[MassRow], const_path: &str) {
    test_inner(all, const_path, 3)
}
fn test_inner(all: &[MassRow], const_path: &str, level: usize) {
    let ts = match level {
        0 => templates(),
        1 => expanded_templates(),
        2 => refined_templates(),
        _ => full_templates(),
    };
    let constants = self::constants::<Hard>(const_path);
    let software_constants = self::constants::<Soft>(const_path);
    let all_rows: Vec<_> = all.iter().filter(|r| r.count == 1).collect();
    let choices: &[EvalModel] = if level >= 3 { &WIDE_MODELS } else { &MODELS };
    let rows: Vec<_> = [0, 1, 5, 22, 56, 89, 90, 178, 500, 1000, 2000, 4231]
        .into_iter()
        .map(|i| all_rows[i])
        .collect();
    let mut tested = 0;
    let mut backend_differences = 0;
    let mut raw_checks = 0;
    for r in rows {
        let mut hc = RowCache::<Hard>::new(r.mean, true);
        let mut sc = RowCache::<Soft>::new(r.mean, true);
        let mut raw = RowCache::<Hard>::new(r.mean, false);
        for t in &ts {
            let expected = roots::<Hard>(&t.graph, r.mean, &constants, choices, &mut hc);
            let mut explicit = vec![];
            let soft = roots::<Soft>(&t.graph, r.mean, &software_constants, choices, &mut sc);
            backend_differences += (expected.len() != soft.len()
                || expected.iter().any(|v| !soft.contains(v)))
                as usize;
            let count = variants_in(t, choices.len());
            for code in 0..count {
                let ms = models_in(code, choices);
                let v = evaluate::<Hard>(&t.graph, r.mean, &constants, &ms, &mut hc);
                if !explicit.contains(&v) {
                    explicit.push(v)
                }
                tested += 1;
                if code % 257 == 0 || code == count - 1 {
                    assert_eq!(
                        v,
                        evaluate::<Hard>(&t.graph, r.mean, &constants, &ms, &mut raw)
                    );
                    raw_checks += 1;
                }
            }
            assert_eq!(expected.len(), explicit.len());
            for v in expected {
                assert!(explicit.contains(&v));
            }
        }
    }
    println!(
        "{}",
        json!({"explicit_template_model_checks":tested,"uncached_reference_checks":raw_checks,"backend_template_set_differences":backend_differences,"note":"Synthetic arithmetic enumeration audit, not fresh Excel validation"})
    );
}

/// Keep recorded builds as separate experimental cohorts; no build-dependent
/// branch is introduced into any candidate or application calculation.
pub fn build_cohorts<F: Fpu>(all: &[MassRow], const_path: &str, output: &str) {
    fs::create_dir_all(output).unwrap();
    let out = Path::new(output);
    let ts = refined_templates();
    let constants = constants::<F>(const_path);
    let mut programs = vec![];
    let mut layout = vec![];
    for t in &ts {
        let count = variants(t);
        layout.push(json!({"template":t.id,"first_candidate_id":programs.len(),"variants":count}));
        for code in 0..count {
            programs.push((t.id, code));
        }
    }
    json_write(out.join("templates.json"), &ts);
    json_write(out.join("program-layout.json"), &layout);
    let mut grouped = std::collections::BTreeMap::<String, Vec<&MassRow>>::new();
    for r in all.iter().filter(|r| r.count == 1) {
        let build = r.raw["reading"]["excel_build"]
            .as_str()
            .unwrap_or("missing")
            .to_string();
        grouped.entry(build).or_default().push(r);
    }
    let mut order: Vec<_> = grouped.keys().cloned().collect();
    order.sort_by_key(|b| (b != "20527", b.clone()));
    let mut summaries = vec![];
    for build in order {
        let rows = &grouped[&build];
        let dir = out.join(format!("build-{build}"));
        fs::create_dir_all(&dir).unwrap();
        let mut kills = BufWriter::new(fs::File::create(dir.join("first-kills.jsonl")).unwrap());
        let mut stats = BufWriter::new(fs::File::create(dir.join("witness-kills.jsonl")).unwrap());
        let mut alive: Vec<_> = (0..programs.len()).collect();
        let mut evaluations = 0;
        let mut kill_summary = vec![];
        for (order, r) in rows.iter().enumerate() {
            let before = alive.len();
            let mut next = vec![];
            let mut cache = RowCache::<F>::new(r.mean, true);
            for id in alive {
                let (tid, code) = programs[id];
                let got = outcome::<F>(evaluate::<F>(
                    &ts[tid].graph,
                    r.mean,
                    &constants,
                    &models(code),
                    &mut cache,
                ));
                evaluations += 1;
                if got == r.raw["expected"].as_str().unwrap() {
                    next.push(id)
                } else {
                    serde_json::to_writer(&mut kills,&json!({"candidate_id":id,"template_id":tid,"model_code":code,"order":order,"trace_index":r.trace_index,"expected":r.raw["expected"],"got":got})).unwrap();
                    writeln!(kills).unwrap();
                }
            }
            let stat = json!({"order":order,"trace_index":r.trace_index,"id":r.raw["id"],"killed":before-next.len(),"remaining":next.len()});
            if before != next.len() {
                println!("build-{build} {stat}");
                kill_summary.push(stat.clone());
            }
            serde_json::to_writer(&mut stats, &stat).unwrap();
            writeln!(stats).unwrap();
            alive = next;
        }
        kills.flush().unwrap();
        stats.flush().unwrap();
        json_write(dir.join("survivors.json"), &alive);
        let summary = json!({"build":build,"programs":programs.len(),"grammar_hash":hash_bytes(&serde_json::to_vec(&ts).unwrap()),"rows":rows.len(),"survivors":alive.len(),"evaluations":evaluations,"first_kills":kill_summary,"complete":true,"backend":std::any::type_name::<F>()});
        json_write(dir.join("summary.json"), &summary);
        println!("{summary}");
        summaries.push(summary);
    }
    json_write(
        out.join("summary.json"),
        &json!({"cohorts":summaries,"policy":"Exogenous recorded-build cohorts only; full mixed-build family stays rejected. No application build selector, newly read answer, or fitted input cutoff."}),
    );
}

fn prediction_codes<F: Fpu>(
    mean: u64,
    ids: &[usize],
    programs: &[(usize, usize)],
    ts: &[Template],
    constants: &Constants,
) -> (Vec<u16>, Vec<String>) {
    let mut cache = RowCache::<F>::new(mean, true);
    let mut names = vec![];
    let mut lookup = HashMap::new();
    let mut codes = Vec::with_capacity(ids.len());
    for id in ids {
        let (t, c) = programs[*id];
        let got = outcome::<F>(evaluate::<F>(
            &ts[t].graph,
            mean,
            constants,
            &models(c),
            &mut cache,
        ));
        let code = *lookup.entry(got.clone()).or_insert_with(|| {
            let n = names.len();
            assert!(n < 65536);
            names.push(got);
            n as u16
        });
        codes.push(code);
    }
    (codes, names)
}
fn refine_labels(labels: &[u32], codes: &[u16]) -> Vec<u32> {
    let mut map = HashMap::new();
    labels
        .iter()
        .zip(codes)
        .map(|(a, b)| {
            let next = map.len() as u32;
            *map.entry((*a, *b)).or_insert(next)
        })
        .collect()
}
fn pair_score(labels: &[u32], codes: &[u16]) -> u64 {
    let mut before = HashMap::<u32, u64>::new();
    let mut after = HashMap::<(u32, u16), u64>::new();
    for (l, c) in labels.iter().zip(codes) {
        *before.entry(*l).or_default() += 1;
        *after.entry((*l, *c)).or_default() += 1;
    }
    before.values().map(|n| n * n).sum::<u64>() - after.values().map(|n| n * n).sum::<u64>()
}
/// Questions cover a declared current-build cohort; the mixed-build rejection
/// remains unchanged. Requests never contain predictions or recorded answers.
pub fn questions(all: &[MassRow], const_path: &str, native: &str, software: &str, output: &str) {
    let out = Path::new(output);
    assert!(!out.exists(), "never overwrite a frozen question packet");
    fs::create_dir_all(out).unwrap();
    let native = Path::new(native);
    let software = Path::new(software);
    let n: Vec<usize> =
        serde_json::from_slice(&fs::read(native.join("build-20527/survivors.json")).unwrap())
            .unwrap();
    let s: std::collections::HashSet<usize> =
        serde_json::from_slice(&fs::read(software.join("build-20527/survivors.json")).unwrap())
            .unwrap();
    let ids: Vec<_> = n.into_iter().filter(|id| s.contains(id)).collect();
    assert!(ids.len() > 1);
    let ts = refined_templates();
    let hc = constants::<Hard>(const_path);
    let sc = constants::<Soft>(const_path);
    let mut programs = vec![];
    for t in &ts {
        for c in 0..variants(t) {
            programs.push((t.id, c));
        }
    }
    let p2 = |e: i32| f64::from_bits(((1023 + e) as u64) << 52);
    let anchors = [
        f64::MIN_POSITIVE,
        p2(-538),
        p2(-512),
        p2(-40),
        p2(-27),
        1e-6,
        0.1,
        0.5,
        1.0,
        2.0,
        5.0,
        10.0,
        100.0,
        700.0,
        720.0,
    ];
    let mut pool = std::collections::BTreeSet::new();
    for a in anchors {
        for delta in [-16i64, -2, -1, 0, 1, 2, 16] {
            let b = a.to_bits().wrapping_add_signed(delta);
            let x = f64::from_bits(b);
            if x.is_normal() && x > 0.0 && x < 1024.0 {
                pool.insert(b);
            }
        }
    }
    let mut state = 0x52375f504f495353u64;
    for _ in 0..128 {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        let e = (state % 1032) as i32 - 1022;
        pool.insert((((e + 1023) as u64) << 52) | (state & 0xfffffffffffff));
    }
    let means: Vec<_> = pool.into_iter().collect();
    let mut codes = vec![];
    let mut names = vec![];
    for mean in &means {
        let (c, n) = prediction_codes::<Hard>(*mean, &ids, &programs, &ts, &hc);
        codes.push(c);
        names.push(n);
    }
    let mut selected = vec![];
    let mut labels = vec![0u32; ids.len()];
    // Exact cross-context controls, then greedy disagreement questions. No
    // observed answer is used in scoring an unanswered question.
    for m in [f64::MIN_POSITIVE, p2(-40), 1e-6, 1.0, 720.0] {
        let i = means.binary_search(&m.to_bits()).unwrap();
        selected.push(i);
        labels = refine_labels(&labels, &codes[i]);
    }
    while selected.len() < 32 {
        let mut best = None;
        for i in 0..means.len() {
            if selected.contains(&i) {
                continue;
            }
            let score = pair_score(&labels, &codes[i]);
            let e = ((means[i] >> 52) & 0x7ff) as i32;
            let distance = selected
                .iter()
                .map(|j| (e - ((means[*j] >> 52) & 0x7ff) as i32).unsigned_abs())
                .min()
                .unwrap();
            let key = (score, distance, names[i].len());
            if best.as_ref().is_none_or(|(k, _)| key > *k) {
                best = Some((key, i));
            }
        }
        let Some(((score, _, _), i)) = best else {
            break;
        };
        if score == 0 {
            break;
        }
        selected.push(i);
        labels = refine_labels(&labels, &codes[i]);
    }
    let mut requests = vec![];
    let mut np = BufWriter::new(fs::File::create(out.join("predictions-native.jsonl")).unwrap());
    let mut sp = BufWriter::new(fs::File::create(out.join("predictions-software.jsonl")).unwrap());
    let mut controls = vec![];
    for (q, i) in selected.iter().enumerate() {
        let mean = means[*i];
        let (soft_codes, soft_names) = prediction_codes::<Soft>(mean, &ids, &programs, &ts, &sc);
        for logical in [true, false] {
            let id = format!(
                "r7-poisson-final-k1-{q:02}-{}",
                if logical { "bool" } else { "number" }
            );
            let flag = if logical {
                json!(false)
            } else {
                json!("0x0000000000000000")
            };
            requests
                .push(json!({"id":id,"args":["0x3ff0000000000000",format!("0x{mean:016x}"),flag]}));
            for (writer, cv, nm) in [
                (&mut np, &codes[*i], &names[*i]),
                (&mut sp, &soft_codes, &soft_names),
            ] {
                let mut groups = vec![Vec::<usize>::new(); nm.len()];
                for (id, code) in ids.iter().zip(cv) {
                    groups[*code as usize].push(*id)
                }
                let partitions: Vec<_> = groups
                    .into_iter()
                    .enumerate()
                    .map(|(j, ids)| json!({"outcome":nm[j],"candidate_ids":ids}))
                    .collect();
                serde_json::to_writer(&mut *writer,&json!({"id":id,"source_build_cohort":"20527","kind":"program predictions, not Excel answers","partitions":partitions})).unwrap();
                writeln!(*writer).unwrap();
            }
            if q < 5 {
                controls.push(id)
            }
        }
    }
    np.flush().unwrap();
    sp.flush().unwrap();
    assert!(requests.len() <= 64);
    json_write(
        out.join("requests.json"),
        &json!({"function":"POISSON.DIST","witnesses":requests}),
    );
    json_write(
        out.join("hypotheses.json"),
        &json!({"grammar_hash":hash_bytes(&serde_json::to_vec(&ts).unwrap()),"candidate_ids":ids,"models":MODELS,"templates_path":native.join("templates.json"),"layout_path":native.join("program-layout.json"),"scope":"Exact on all88 build20527 count-one observations on each backend. Rejected globally; older builds and argument encoding are confounded, not explained."}),
    );
    let summary = json!({"questions":requests.len(),"means":selected.len(),"candidate_programs":ids.len(),"pool_means":means.len(),"discovery_rows":all.iter().filter(|r|r.count==1&&r.raw["reading"]["excel_build"]=="20527").count(),"prediction_signatures":labels.iter().copied().collect::<std::collections::HashSet<_>>().len(),"signature_scope":"Finite signatures on the selected pool, not global program equivalence.","control_ids":controls,"selection":"Five cross-context/domain controls (including mean720) plus exact greedy pair disagreement; exponent distance breaks ties; no answer labels in selection.","answers":0,"reader":"cargo xtask witness read --func POISSON.DIST --inputs requests.json (Windows; includes rechecks, omit --new-only)"});
    json_write(out.join("selection.json"), &summary);
    println!("{summary}");
}

/// Expose computed numerator states to the existing exact interval engine.
/// This projects the original three-model normalization form only. An unknown
/// denominator is a constraint variable, never an inferred product constant.
pub fn normalizer_matrix<F: Fpu>(all: &[MassRow], const_path: &str, output: &str) {
    let out = Path::new(output);
    assert!(!out.exists());
    fs::create_dir_all(out).unwrap();
    let constants = constants::<F>(const_path);
    let mut graphs = vec![];
    for t in refined_templates().iter().filter(|t| t.normalization == 0) {
        let Node::Div { site: 3, a, .. } = &t.graph else {
            panic!("normalization form")
        };
        if !graphs.contains(a.as_ref()) {
            graphs.push(a.as_ref().clone())
        }
    }
    assert_eq!(graphs.len(), 672);
    json_write(out.join("numerator-graphs.json"), &graphs);
    let programs = graphs.len() * 27;
    assert_eq!(programs, 18144);
    let rows: Vec<_> = all
        .iter()
        .filter(|r| r.count == 1 && r.raw["reading"]["excel_build"] == "20131")
        .collect();
    let mut file = BufWriter::new(fs::File::create(out.join("matrix.jsonl")).unwrap());
    for (order, r) in rows.iter().take(64).enumerate() {
        let mut cache = RowCache::<F>::new(r.mean, true);
        let mut groups = vec![];
        let mut index = HashMap::<Val, usize>::new();
        for (g, graph) in graphs.iter().enumerate() {
            for code in 0..27 {
                let value = evaluate::<F>(graph, r.mean, &constants, &models(code), &mut cache);
                let next = groups.len();
                let i = *index.entry(value).or_insert_with(|| {
                    groups.push((value, Vec::<usize>::new()));
                    next
                });
                groups[i].1.push(g * 27 + code);
            }
        }
        let values: Vec<_> = groups
            .into_iter()
            .map(|(v, ids)| json!({"value":v.exact_encoding(),"numerator_ids":ids}))
            .collect();
        serde_json::to_writer(
            &mut file,
            &json!({"order":order,"trace_index":r.trace_index,"witness":r.raw,"values":values}),
        )
        .unwrap();
        writeln!(file).unwrap();
    }
    file.flush().unwrap();
    let summary = json!({"numerator_templates":graphs.len(),"numerator_programs":programs,"denominator_profiles":programs*3,"exported_rows":rows.len().min(64),"available_build20131_rows":rows.len(),"models":MODELS,"backend":std::any::type_name::<F>(),"kind":"Computed research numerator states plus original recorded witness references; no new Excel answers"});
    json_write(out.join("summary.json"), &summary);
    println!("{summary}");
}
