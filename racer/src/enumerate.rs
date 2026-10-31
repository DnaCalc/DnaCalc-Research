//! Bounded structural templates and the racer's Cartesian model assignments.
use crate::dsl::*;

pub fn template(numerator: usize, combine: usize) -> Graph {
    let mut b = GraphBuilder::default();
    // Leaves: rate, PV, FV, canonical timing (0 or 1), frozen expm1.
    let r = b.strict(Op::Arg(0));
    let pv = b.strict(Op::Arg(1));
    let fv = b.strict(Op::Arg(2));
    let timing = b.strict(Op::Arg(3));
    let em = b.strict(Op::Arg(4));
    let one = b.strict(Op::Const {
        bits: 1f64.to_bits(),
    });
    let n = if numerator == 0 {
        let v = b.strict(Op::Add(one, em));
        let f = b.strict(Op::Mul(fv, v));
        b.strict(Op::Add(pv, f))
    } else {
        let sum = b.strict(Op::Add(pv, fv));
        let f = b.strict(Op::Mul(fv, em));
        b.strict(Op::Add(sum, f))
    };
    let rt = b.strict(Op::Mul(r, timing));
    let tf = b.strict(Op::Add(one, rt));
    let output = match combine {
        0 => {
            let q = b.strict(Op::Div(n, em));
            let q = b.strict(Op::Div(q, tf));
            b.strict(Op::Mul(q, r))
        }
        1 => {
            let nr = b.strict(Op::Mul(n, r));
            let et = b.strict(Op::Mul(em, tf));
            b.strict(Op::Div(nr, et))
        }
        2 => {
            let nt = b.strict(Op::Div(n, tf));
            let re = b.strict(Op::Div(r, em));
            b.strict(Op::Mul(nt, re))
        }
        _ => unreachable!(),
    };
    b.finish(output)
}

pub fn model_variants(g: &Graph) -> Vec<([EvalModel; 8], Graph)> {
    let ids: Vec<_> = g
        .nodes
        .iter()
        .enumerate()
        .filter_map(|(i, n)| matches!(n.op, Op::Add(..) | Op::Mul(..) | Op::Div(..)).then_some(i))
        .collect();
    assert_eq!(ids.len(), 8);
    let space = 3usize.pow(ids.len() as u32);
    assert!(space <= 100_000);
    (0..space)
        .map(|mut assignment| {
            let mut graph = g.clone();
            let mut models = [EvalModel::Strict; 8];
            for (slot, id) in ids.iter().enumerate() {
                models[slot] = MODELS[assignment % 3];
                assignment /= 3;
                graph.nodes[*id].model = models[slot];
            }
            (models, graph)
        })
        .collect()
}

pub fn family() -> Vec<Candidate> {
    let mut out = Vec::new();
    for numerator in 0..2 {
        for combine in 0..3 {
            for (models, graph) in model_variants(&template(numerator, combine)) {
                out.push(Candidate {
                    id: out.len(),
                    numerator,
                    combine,
                    models,
                    graph,
                });
            }
        }
    }
    assert_eq!(out.len(), 39_366);
    // Hashes accelerate lookup; graph equality is checked before deduplication.
    let mut buckets = std::collections::HashMap::<String, Vec<&Graph>>::new();
    for c in &out {
        let bucket = buckets.entry(c.hash()).or_default();
        assert!(!bucket.contains(&&c.graph), "duplicate serialized graph");
        bucket.push(&c.graph);
    }
    out
}
