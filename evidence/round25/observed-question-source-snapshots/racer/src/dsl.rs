//! Arithmetic subset of OxFunc's serialized calculation graphs.
use serde::{Deserialize, Serialize};

pub type NodeId = usize;

#[derive(Serialize, Deserialize, Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[serde(rename_all = "snake_case")]
pub enum EvalModel {
    Strict,
    StoredPc64,
    ResidentPc64,
    StoredPc53,
    ResidentPc53,
}
pub const MODELS: [EvalModel; 3] = [
    EvalModel::Strict,
    EvalModel::StoredPc64,
    EvalModel::ResidentPc64,
];
pub const WIDE_MODELS: [EvalModel; 5] = [
    EvalModel::Strict,
    EvalModel::StoredPc64,
    EvalModel::ResidentPc64,
    EvalModel::StoredPc53,
    EvalModel::ResidentPc53,
];

#[derive(Serialize, Deserialize, Clone, Copy, Debug, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Op {
    Arg(usize),
    Const { bits: u64 },
    Add(NodeId, NodeId),
    Mul(NodeId, NodeId),
    Div(NodeId, NodeId),
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct Node {
    pub op: Op,
    pub model: EvalModel,
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct Graph {
    pub nodes: Vec<Node>,
    pub output: NodeId,
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct Candidate {
    pub id: usize,
    pub numerator: usize,
    pub combine: usize,
    pub models: [EvalModel; 8],
    pub graph: Graph,
}

pub fn hash_bytes(bytes: &[u8]) -> String {
    let mut h = 0xcbf29ce484222325u64;
    for b in bytes {
        h ^= *b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("fnv1a64:{h:016x}")
}

impl Candidate {
    pub fn hash(&self) -> String {
        hash_bytes(&serde_json::to_vec(&self.graph).unwrap())
    }
}

#[derive(Default)]
pub struct GraphBuilder {
    nodes: Vec<Node>,
}
impl GraphBuilder {
    pub fn push(&mut self, op: Op, model: EvalModel) -> NodeId {
        self.nodes.push(Node { op, model });
        self.nodes.len() - 1
    }
    pub fn strict(&mut self, op: Op) -> NodeId {
        self.push(op, EvalModel::Strict)
    }
    pub fn finish(self, output: NodeId) -> Graph {
        Graph {
            nodes: self.nodes,
            output,
        }
    }
}

/// Declared accumulation trees, shared by the cumulative research sketches.
pub fn sum_graph(
    g: &mut GraphBuilder,
    terms: &[NodeId],
    order: usize,
    model: EvalModel,
    zero: NodeId,
) -> NodeId {
    fn paired(g: &mut GraphBuilder, terms: &[NodeId], model: EvalModel, zero: NodeId) -> NodeId {
        match terms.len() {
            0 => zero,
            1 => terms[0],
            n => {
                let a = paired(g, &terms[..n / 2], model, zero);
                let b = paired(g, &terms[n / 2..], model, zero);
                g.push(Op::Add(a, b), model)
            }
        }
    }
    match order {
        0 => terms
            .iter()
            .fold(zero, |a, &b| g.push(Op::Add(a, b), model)),
        1 => terms
            .iter()
            .rev()
            .fold(zero, |a, &b| g.push(Op::Add(a, b), model)),
        2 => paired(g, terms, model, zero),
        _ => panic!("undeclared accumulation order"), // ubs:ignore — fail closed on an out-of-grammar research profile
    }
}
