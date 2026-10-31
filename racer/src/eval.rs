//! Full tagged values flow until a binary64 store; strict consumers store inputs.
use crate::{
    dsl::*,
    ported::excel_numeric::{Extended, Fpu, research_load},
};
pub const CW: u16 = 0x133f;

#[derive(Clone, Copy, Debug)]
pub enum Val {
    F64(u64),
    Ext(Extended),
}
impl PartialEq for Val {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::F64(a), Self::F64(b)) => a == b,
            (Self::Ext(a), Self::Ext(b)) => a.0 == b.0,
            _ => false,
        }
    }
}
impl Eq for Val {}
impl std::hash::Hash for Val {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        match self {
            Self::F64(b) => std::hash::Hash::hash(&(0u8, *b), state),
            Self::Ext(e) => std::hash::Hash::hash(&(1u8, e.0), state),
        }
    }
}
impl Val {
    pub fn f64(x: f64) -> Self {
        Self::F64(x.to_bits())
    }
    pub fn to_f64<F: Fpu>(self) -> f64 {
        match self {
            Self::F64(b) => f64::from_bits(b),
            Self::Ext(e) => F::store(&e, CW),
        }
    }
    pub fn to_ext(self) -> Extended {
        match self {
            Self::F64(b) => research_load(f64::from_bits(b)),
            Self::Ext(e) => e,
        }
    }
    pub fn exact_encoding(self) -> String {
        match self {
            Self::F64(bits) => format!("binary64:{bits:016x}"),
            Self::Ext(e) => format!(
                "ext80-le:{}",
                e.0.iter().map(|b| format!("{b:02x}")).collect::<String>()
            ),
        }
    }
}

#[inline]
pub fn binop<F: Fpu>(op: Op, model: EvalModel, a: Val, b: Val) -> Val {
    if model == EvalModel::Strict {
        let a = a.to_f64::<F>();
        let b = b.to_f64::<F>();
        return Val::f64(match op {
            Op::Add(..) => a + b,
            Op::Mul(..) => a * b,
            Op::Div(..) => a / b,
            _ => unreachable!(),
        });
    }
    let a = a.to_ext();
    let b = b.to_ext();
    let cw = if matches!(model, EvalModel::StoredPc53 | EvalModel::ResidentPc53) {
        (CW & !0x0300) | 0x0200
    } else {
        CW
    };
    let x = match op {
        Op::Add(..) => F::add(&a, &b, cw),
        Op::Mul(..) => F::mul(&a, &b, cw),
        Op::Div(..) => F::div(&a, &b, cw),
        _ => unreachable!(),
    };
    if matches!(model, EvalModel::StoredPc64 | EvalModel::StoredPc53) {
        Val::f64(F::store(&x, CW))
    } else {
        Val::Ext(x)
    }
}

pub fn precision_controls() {
    use crate::ported::excel_numeric::{Hard, Soft};
    let epsilon = Val::f64(f64::from_bits(((1023 - 63) as u64) << 52));
    let a = binop::<Hard>(
        Op::Add(0, 1),
        EvalModel::ResidentPc64,
        Val::f64(1.0),
        epsilon,
    );
    let minus_one = Val::f64(-1.0);
    assert_eq!(
        binop::<Hard>(Op::Add(0, 1), EvalModel::Strict, a, minus_one).to_f64::<Hard>(),
        0.0
    );
    for model in [EvalModel::StoredPc53, EvalModel::ResidentPc53] {
        let v = binop::<Hard>(Op::Add(0, 1), model, a, minus_one);
        assert_eq!(
            v.to_f64::<Hard>().to_bits(),
            epsilon.to_f64::<Hard>().to_bits()
        );
        assert_eq!(
            v.exact_encoding(),
            binop::<Soft>(Op::Add(0, 1), model, a, minus_one).exact_encoding()
        );
    }
    let mut state = 0x504335335f5237u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut checks = 0;
    for _ in 0..20000 {
        let ab = (next() & 0x800fffffffffffff) | (((next() % 2046) + 1) << 52);
        let bb = (next() & 0x800fffffffffffff) | (((next() % 2046) + 1) << 52);
        let mut a = Val::F64(ab).to_ext();
        let mut b = Val::F64(bb).to_ext();
        // Add precision below the binary64 significand without changing the
        // exponent or explicit integer bit. These are raw arithmetic controls.
        a.0[0] = next() as u8;
        b.0[0] = next() as u8;
        for op in [Op::Add(0, 1), Op::Mul(0, 1), Op::Div(0, 1)] {
            for model in [EvalModel::StoredPc53, EvalModel::ResidentPc53] {
                let h = binop::<Hard>(op, model, Val::Ext(a), Val::Ext(b));
                let s = binop::<Soft>(op, model, Val::Ext(a), Val::Ext(b));
                assert_eq!(h.exact_encoding(), s.exact_encoding());
                checks += 1;
                for continuation in MODELS {
                    assert_eq!(
                        binop::<Hard>(Op::Mul(0, 1), continuation, h, Val::f64(0.75))
                            .exact_encoding(),
                        binop::<Soft>(Op::Mul(0, 1), continuation, s, Val::f64(0.75))
                            .exact_encoding()
                    );
                    checks += 1;
                }
            }
        }
    }
    println!(
        "PC53 controls: {checks} native/software arithmetic/continuation checks; full-operand cancellation preserved"
    );
}

pub fn eval_graph<F: Fpu>(g: &Graph, args: &[f64; 5]) -> Val {
    let mut values = [Val::F64(0); 14];
    assert!(g.nodes.len() <= values.len());
    for (i, n) in g.nodes.iter().enumerate() {
        values[i] = match n.op {
            Op::Arg(a) => Val::f64(args[a]),
            Op::Const { bits } => Val::F64(bits),
            Op::Add(a, b) | Op::Mul(a, b) | Op::Div(a, b) => {
                assert!(a < i && b < i);
                binop::<F>(n.op, n.model, values[a], values[b])
            }
        };
    }
    values[g.output]
}

pub fn outcome<F: Fpu>(x: Val) -> String {
    let x = x.to_f64::<F>();
    if !x.is_finite() {
        "error:Num".into()
    } else if x.is_subnormal() {
        "fault:NotImplemented".into()
    } else {
        format!("0x{:016x}", if x == 0.0 { 0 } else { x.to_bits() })
    }
}

pub fn publication_bits<F: Fpu>(x: Val) -> u64 {
    let x = x.to_f64::<F>();
    if !x.is_finite() {
        f64::NAN.to_bits()
    } else if x.is_subnormal() {
        f64::NAN.to_bits() + 1
    } else if x == 0.0 {
        0
    } else {
        x.to_bits()
    }
}

pub fn bits_outcome(bits: u64) -> String {
    if bits == f64::NAN.to_bits() {
        "error:Num".into()
    } else if bits == f64::NAN.to_bits() + 1 {
        "fault:NotImplemented".into()
    } else {
        format!("0x{bits:016x}")
    }
}

/// Exact set propagation for the whole family, including eliminated candidates.
/// Each computed subtree is used once; only immutable leaves are shared. Thus
/// independent Cartesian joins cannot invent inconsistent model assignments.
pub fn root_values<F: Fpu>(args: &[f64; 5]) -> Vec<Val> {
    let mut roots = Vec::new();
    for numerator in 0..2 {
        for combine in 0..3 {
            let g = crate::enumerate::template(numerator, combine);
            let mut values: Vec<Vec<Val>> = Vec::new();
            for n in &g.nodes {
                let mut out = Vec::new();
                match n.op {
                    Op::Arg(a) => out.push(Val::f64(args[a])),
                    Op::Const { bits } => out.push(Val::F64(bits)),
                    Op::Add(a, b) | Op::Mul(a, b) | Op::Div(a, b) => {
                        for &av in &values[a] {
                            for &bv in &values[b] {
                                for model in MODELS {
                                    let x = binop::<F>(n.op, model, av, bv);
                                    if !out.contains(&x) {
                                        out.push(x);
                                    }
                                }
                            }
                        }
                    }
                }
                values.push(out);
            }
            for &v in &values[g.output] {
                if !roots.contains(&v) {
                    roots.push(v);
                }
            }
        }
    }
    roots
}

pub fn publish_roots<F: Fpu>(roots: &[Val]) -> Vec<u64> {
    let mut published = Vec::new();
    for &v in roots {
        let bits = publication_bits::<F>(v);
        if !published.contains(&bits) {
            published.push(bits);
        }
    }
    published.sort_unstable();
    published
}

pub fn reachable<F: Fpu>(args: &[f64; 5]) -> Vec<u64> {
    publish_roots::<F>(&root_values::<F>(args))
}
