use crate::ported::excel_numeric::*;
use std::io::BufRead;
pub fn hx(s: &str) -> f64 { f64::from_bits(u64::from_str_radix(s, 16).unwrap()) }
pub fn ld(x: f64) -> XExt { xload(x) }
pub fn st(x: &XExt) -> f64 { xstore(x) }
pub fn alnrel(rate: f64) -> f64 {
    if rate.abs() < 0.375 { excel_log1p_atanh(rate) } else { excel_log(excel_x87_add(1.0, rate)) }
}
pub fn rows(file: &str) -> Vec<Vec<f64>> {
    std::io::BufReader::new(std::fs::File::open(file).unwrap()).lines()
        .map(|l| l.unwrap().split_whitespace().map(hx).collect()).collect()
}
pub fn report(file: &str, f: &dyn Fn(&[f64]) -> Vec<(String, f64)>) {
    let mut names: Vec<String> = Vec::new();
    let mut counts: Vec<usize> = Vec::new();
    let mut total = 0usize;
    for p in rows(file) {
        let want = p[p.len() - 1];
        let c = f(&p);
        if names.is_empty() { names = c.iter().map(|x| x.0.clone()).collect(); counts = vec![0; c.len()]; }
        total += 1;
        for (i, (_, v)) in c.iter().enumerate() { if v.to_bits() == want.to_bits() { counts[i] += 1; } }
    }
    println!("rows {total}");
    let mut v: Vec<_> = names.iter().zip(counts.iter()).collect();
    v.sort_by_key(|x| *x.1);
    for (nm, c) in v { println!("{nm:24} {c}"); }
}

fn singles(p: &[f64]) -> Vec<(String, f64)> {
    let (r, n, pv, k) = (p[0], p[1], p[2], p[3]);
    let ls: Vec<(&str, f64)> = vec![
        ("alnrel", alnrel(r)),
        ("lnst", excel_log(excel_x87_add(1.0, r))),
        ("lnpl", excel_log(1.0 + r)),
        ("log1p", excel_log1p(r)),
        ("atanh", excel_log1p_atanh(r)),
    ];
    let mut c: Vec<(String, f64)> = Vec::new();
    let l1 = alnrel(r);
    let emn1 = excel_expm1_spilled(-(n * l1));
    let q1 = excel_x87_div(pv, emn1);
    let pmt1 = q1 * r;
    let m = n - k + 1.0;
    for (ln, l) in ls.iter() {
        let l = *l;
        let emn = excel_expm1_spilled(-(n * l));
        let q = excel_x87_div(pv, emn);
        let pmt = q * r;
        let t = -(m * l);
        let u = excel_exp(t);
        let g = excel_exp(m * l);
        c.push((format!("{ln}_A"), excel_x87_mul(excel_x87_mul(q, u), r)));
        c.push((format!("{ln}_B"), excel_x87_mul(pmt, u)));
        c.push((format!("{ln}_G"), excel_x87_div(pmt, g)));
        c.push((format!("{ln}_pmt1_B"), excel_x87_mul(pmt1, u)));
        c.push((format!("{ln}_pmt1_G"), excel_x87_div(pmt1, g)));
        c.push((format!("{ln}_q1_A"), excel_x87_mul(excel_x87_mul(q1, u), r)));
    }
    for (pn, u, g) in [
        ("crtpow", crt_pow(1.0 + r, -m), crt_pow(1.0 + r, m)),
        ("chain", excel_pow_chain(1.0 + r, -m), excel_pow_chain(1.0 + r, m)),
        ("x87d", excel_pow_x87_direct(1.0 + r, -m), excel_pow_x87_direct(1.0 + r, m)),
    ] {
        c.push((format!("{pn}_pmt1_B"), excel_x87_mul(pmt1, u)));
        c.push((format!("{pn}_pmt1_G"), excel_x87_div(pmt1, g)));
        c.push((format!("{pn}_q1_A"), excel_x87_mul(excel_x87_mul(q1, u), r)));
        c.push((format!("{pn}_pmt1_Bp"), pmt1 * u));
        c.push((format!("{pn}_pmt1_Gp"), pmt1 / g));
    }
    c
}


fn lscan(file: &str) {
    let all = rows(file);
    let mut groups: Vec<Vec<Vec<f64>>> = Vec::new();
    for p in all {
        if let Some(g) = groups.last_mut() { if g[0][0] == p[0] && g[0][1] == p[1] && g[0][2] == p[2] { g.push(p); continue; } }
        groups.push(vec![p]);
    }
    for g in groups {
        let (r, n, pv) = (g[0][0], g[0][1], g[0][2]);
        let pmt = g[0][6];
        let base = alnrel(r);
        println!("loan r={r} n={n} pv={pv} rows={} log1p-alnrel ulps={}", g.len(), (excel_log1p(r).to_bits() as i64) - (base.to_bits() as i64));
        for off in -6i64..=6 {
            let l = f64::from_bits((base.to_bits() as i64 + off) as u64);
            let mut cnt = [0usize; 8];
            for p in &g {
                let want = p[7];
                let m = n - p[3] + 1.0;
                let t = -(m * l);
                let u = excel_exp(t);
                let ue = xexp(&ld(t));
                let te = xchs(&xmul(&ld(m), &ld(l)), XCW);
                let ue2 = xexp(&te);
                let gr = excel_exp(m * l);
                let emn = excel_expm1_spilled(-(n * l));
                let q = excel_x87_div(pv, emn);
                let c = [excel_x87_mul(pmt, u), pmt * u, excel_x87_div(pmt, gr), pmt / gr,
                         st(&xmul(&ld(pmt), &ue)), st(&xmul(&ld(pmt), &ue2)),
                         excel_x87_mul(excel_x87_mul(q, u), r), st(&xdiv(&ld(pmt), &xexp(&ld(m*l))))];
                for i in 0..8 { if c[i].to_bits() == want.to_bits() { cnt[i] += 1; } }
            }
            println!("  off {off:3} B {:4} Bp {:4} G {:4} Gp {:4} E {:4} J {:4} A {:4} H {:4}", cnt[0], cnt[1], cnt[2], cnt[3], cnt[4], cnt[5], cnt[6], cnt[7]);
        }
    }
}


fn resid(file: &str) {
    for p in rows(file) {
        let (r, n, pv, k, pmt, want) = (p[0], p[1], p[2], p[3], p[6], p[7]);
        let l = alnrel(r);
        let m = n - k + 1.0;
        let b = excel_x87_mul(pmt, excel_exp(-(m * l)));
        let emn = excel_expm1_spilled(-(n * l));
        let q = excel_x87_div(pv, emn);
        let a = excel_x87_mul(excel_x87_mul(q, excel_exp(-(m * l))), r);
        println!("{r} {n} {k} {m} {} {}", want.to_bits() as i64 - b.to_bits() as i64, want.to_bits() as i64 - a.to_bits() as i64);
    }
}


pub fn xbits(x: &XExt) -> String { let mut s = String::new(); for b in x.0.iter().rev() { s.push_str(&format!("{b:02x}")); } s }
fn dumpu(file: &str) {
    for p in rows(file) {
        let (r, n, k) = (p[0], p[1], p[3]);
        let l = alnrel(r);
        let m = n - k + 1.0;
        let t = -(m * l);
        let u = excel_exp(t);
        let ue = xexp(&ld(t));
        let g = excel_exp(-t);
        let ge = xexp(&ld(-t));
        let te = xchs(&xmul(&ld(m), &ld(l)), XCW);
        let ue2 = xexp(&te);
        println!("{:016x} {} {:016x} {} {}", u.to_bits(), xbits(&ue), g.to_bits(), xbits(&ge), xbits(&ue2));
    }
}


fn dumpl(file: &str) {
    let mut last = f64::NAN;
    for p in rows(file) { if p[0] != last { last = p[0]; println!("{} {:016x}", p[0], alnrel(p[0]).to_bits()); } }
}


fn rz(file: &str) {
    report(file, &|p: &[f64]| {
        let (r, n, pv, k, pmt) = (p[0], p[1], p[2], p[3], p[6]);
        let l = alnrel(r);
        let m = n - k + 1.0;
        let t = -(m * l);
        let urz = excel_exp_rz(t);
        let u = excel_exp(t);
        let ue = xexp(&ld(t));
        let emn = excel_expm1_spilled(-(n * l));
        let q = excel_x87_div(pv, emn);
        const RZ: u16 = 0x133F | 0x0C00;
        let g = excel_exp(-t);
        let grz = excel_exp_rz(-t);
        vec![
            ("pmt*urz x87".into(), excel_x87_mul(pmt, urz)),
            ("pmt*urz".into(), pmt * urz),
            ("pmt*u rzstore".into(), xstore_rz(&xmul(&ld(pmt), &ld(u)))),
            ("pmt*urz rzstore".into(), xstore_rz(&xmul(&ld(pmt), &ld(urz)))),
            ("pmt*ue rzstore".into(), xstore_rz(&xmul(&ld(pmt), &ue))),
            ("pmt*ue rzall".into(), xstore_rz(&xmulc(&ld(pmt), &ue, RZ))),
            ("q*urz*r".into(), excel_x87_mul(excel_x87_mul(q, urz), r)),
            ("pmt/g".into(), excel_x87_div(pmt, g)),
            ("pmt/grz".into(), excel_x87_div(pmt, grz)),
            ("pmt/g rz".into(), xstore_rz(&xdiv(&ld(pmt), &ld(g)))),
            ("pmt*u".into(), excel_x87_mul(pmt, u)),
        ]
    });
}


fn diffs(file: &str) {
    report(file, &|p: &[f64]| {
        let (r, n, pv, s, e) = (p[0], p[1], p[2], p[3], p[4]);
        let l = alnrel(r);
        let emn = excel_expm1_spilled(-(n * l));
        let q = excel_x87_div(pv, emn);
        let qe = xdiv(&ld(pv), &ld(emn));
        let a = n - e; let b = n - s + 1.0;
        let em = |j: f64| if j == 0.0 { 0.0 } else { excel_expm1_spilled(-(j * l)) };
        let emi = |j: f64| if j == 0.0 { 0.0 } else { excel_expm1_internal(-(j * l)) };
        let u = |j: f64| excel_exp(-(j * l));
        let ue = |j: f64| xexp(&ld(-(j * l)));
        let mut v: Vec<(String, f64)> = vec![
            ("q*(em_a-em_b) x87".into(), excel_x87_mul(q, excel_x87_sub(em(a), em(b)))),
            ("q*(emi_a-emi_b) x87".into(), excel_x87_mul(q, excel_x87_sub(emi(a), emi(b)))),
            ("q*(u_a-u_b) x87".into(), excel_x87_mul(q, excel_x87_sub(u(a), u(b)))),
            ("q*(u_a-u_b) ext".into(), st(&xmul(&ld(q), &xsub(&ld(u(a)), &ld(u(b)))))),
            ("q*(ue_a-ue_b) ext".into(), st(&xmul(&ld(q), &xsub(&ue(a), &ue(b))))),
            ("qe*(ue_a-ue_b) ext".into(), st(&xmul(&qe, &xsub(&ue(a), &ue(b))))),
            ("q*(em_a-em_b) ext".into(), st(&xmul(&ld(q), &xsub(&ld(em(a)), &ld(em(b)))))),
            ("q*ua-q*ub".into(), excel_x87_sub(excel_x87_mul(q, u(a)), excel_x87_mul(q, u(b)))),
            ("q*ema-q*emb".into(), excel_x87_sub(excel_x87_mul(q, em(a)), excel_x87_mul(q, em(b)))),
            ("pv*(em_a-em_b)/emn".into(), excel_x87_div(excel_x87_mul(pv, excel_x87_sub(em(a), em(b))), emn)),
            ("(em_a-em_b)/emn*pv".into(), excel_x87_mul(excel_x87_div(excel_x87_sub(em(a), em(b)), emn), pv)),
        ];
        v
    });
}


fn groups(file: &str) -> Vec<Vec<Vec<f64>>> {
    let mut groups: Vec<Vec<Vec<f64>>> = Vec::new();
    for p in rows(file) {
        if let Some(g) = groups.last_mut() { if g[0][0] == p[0] && g[0][1] == p[1] && g[0][2] == p[2] { g.push(p); continue; } }
        groups.push(vec![p]);
    }
    groups
}
fn cscan(file: &str) {
    for g in groups(file) {
        let (r, n, pv) = (g[0][0], g[0][1], g[0][2]);
        let pmt = g[0][6];
        let l = alnrel(r);
        let emn = excel_expm1_spilled(-(n * l));
        let q = excel_x87_div(pv, emn);
        println!("loan r={r} n={n} pv={pv} rows={}", g.len());
        let mut best: Vec<(usize, String)> = Vec::new();
        for (cn, c0) in [("pmt", pmt), ("q", q)] {
            for off in -40i64..=40 {
                let c = f64::from_bits((c0.to_bits() as i64 + off) as u64);
                let mut cnt = [0usize; 9];
                for p in &g {
                    let want = p[7].to_bits();
                    let m = n - p[3] + 1.0;
                    let t = -(m * l);
                    let u = excel_exp(t);
                    let ue = xexp(&ld(t));
                    let gr = excel_exp(-t);
                    let vals = [excel_x87_mul(c, u), st(&xmul(&ld(c), &ue)), excel_x87_mul(excel_x87_mul(c, u), r),
                        excel_x87_mul(c, excel_x87_mul(u, r)), excel_x87_mul(excel_x87_mul(c, r), u), excel_x87_div(c, gr),
                        st(&xmul(&xmul(&ld(c), &ue), &ld(r))), st(&xmul(&ld(c), &ld(excel_exp_rz(t)))), excel_x87_mul(excel_x87_div(c, gr), r)];
                    for i in 0..9 { if vals[i].to_bits() == want { cnt[i] += 1; } }
                }
                for i in 0..9 { best.push((cnt[i], format!("{cn}{off:+} form{i}"))); }
            }
        }
        best.sort();
        for b in best.iter().rev().take(6) { println!("   {} {}", b.0, b.1); }
    }
}


fn crt(file: &str) {
    report(file, &|p: &[f64]| {
        let (r, n, pv, k, pmt) = (p[0], p[1], p[2], p[3], p[6]);
        let m = n - k + 1.0;
        let mut v: Vec<(String, f64)> = Vec::new();
        let l1 = alnrel(r);
        let emn = excel_expm1_spilled(-(n * l1));
        let q = excel_x87_div(pv, emn);
        for (ln, l) in [("alnrel", l1), ("crtlog", crt_log(1.0 + r)), ("log1p", excel_log1p(r)), ("xlog", excel_log(1.0 + r))] {
            let t = -(m * l);
            for (en, u) in [("crtexp", crt_exp(t)), ("exp", excel_exp(t)), ("portable", exp_portable(t))] {
                v.push((format!("{ln}/{en} pmt*u"), pmt * u));
                v.push((format!("{ln}/{en} pmt*u x87"), excel_x87_mul(pmt, u)));
                v.push((format!("{ln}/{en} q*u*r"), q * u * r));
                v.push((format!("{ln}/{en} q*u*r x87"), excel_x87_mul(excel_x87_mul(q, u), r)));
                v.push((format!("{ln}/{en} pmt/g"), pmt / crt_exp(-t)));
            }
        }
        for (pn, u) in [("crtpow", crt_pow(1.0 + r, -m)), ("crtpowinv", 1.0 / crt_pow(1.0 + r, m))] {
            v.push((format!("{pn} pmt*u"), pmt * u));
            v.push((format!("{pn} pmt*u x87"), excel_x87_mul(pmt, u)));
            v.push((format!("{pn} pmt/g"), pmt / crt_pow(1.0 + r, m)));
            v.push((format!("{pn} q*u*r"), q * u * r));
        }
        v
    });
}


fn show(file: &str) {
    for p in rows(file) {
        let (r, n, pv, k, pmt, want) = (p[0], p[1], p[2], p[3], p[6], p[7]);
        let l = alnrel(r);
        let m = n - k + 1.0;
        let t = -(m * l);
        let u = excel_exp(t);
        let b = excel_x87_mul(pmt, u);
        let ue2 = xexp(&xchs(&xmul(&ld(m), &ld(l)), XCW));
        let j = st(&xmul(&ld(pmt), &ue2));
        // t accumulated: start -(n*l), add l (k-1) times in double
        let mut ta = -(n * l); for _ in 1..(k as i64) { ta += l; }
        let mut tx = xchs(&xmul(&ld(n), &ld(l)), XCW); for _ in 1..(k as i64) { tx = xadd(&tx, &ld(l)); }
        let ua = excel_x87_mul(pmt, excel_exp(ta));
        let ux = st(&xmul(&ld(pmt), &xexp(&tx)));
        let d = |x: f64| x.to_bits() as i64 - want.to_bits() as i64;
        println!("r={r} n={n} k={k} m={m} want={:016x} B{:+} J{:+} acc{:+} accx{:+} t0={:016x} l={:016x}", want.to_bits(), d(b), d(j), d(ua), d(ux), t.to_bits(), l.to_bits());
    }
}


fn lext(file: &str) {
    report(file, &|p: &[f64]| {
        let (r, n, pv, k, pmt) = (p[0], p[1], p[2], p[3], p[6]);
        let m = n - k + 1.0;
        let base = excel_x87_add(1.0, r);
        let lx = xlog(&ld(base));
        let l = alnrel(r);
        let mut v: Vec<(String, f64)> = Vec::new();
        let t1 = -st(&xmul(&ld(m), &lx));
        let tx = xchs(&xmul(&ld(m), &lx), XCW);
        let t2 = st(&tx);
        v.push(("B t=RN(m*Lext)".into(), excel_x87_mul(pmt, excel_exp(t1))));
        v.push(("ext t exp".into(), st(&xmul(&ld(pmt), &xexp(&tx)))));
        v.push(("ext t exp stored u".into(), excel_x87_mul(pmt, st(&xexp(&tx)))));
        v.push(("B".into(), excel_x87_mul(pmt, excel_exp(-(m * l)))));
        // pow-like: 2^(y*log2 x) with y=-m
        let l2 = Native_fyl2x_dummy(m, base);
        v.push(("pow x87 direct".into(), excel_x87_mul(pmt, excel_pow_x87_direct(base, -m))));
        v.push(("pow x87 direct inv".into(), excel_x87_div(pmt, excel_pow_x87_direct(base, m))));
        v.push(("pow chain".into(), excel_x87_mul(pmt, excel_pow_chain(base, -m))));
        v.push(("pow chain inv".into(), excel_x87_div(pmt, excel_pow_chain(base, m))));
        v.push(("crt pow".into(), excel_x87_mul(pmt, crt_pow(base, -m))));
        v.push(("crt pow inv".into(), excel_x87_div(pmt, crt_pow(base, m))));
        let _ = l2;
        v
    });
}
#[allow(non_snake_case)]
fn Native_fyl2x_dummy(_m: f64, _b: f64) -> f64 { 0.0 }


fn funcof(file: &str, r: f64, n: f64) {
    use std::collections::HashMap;
    let l = alnrel(r);
    let emn = excel_expm1_spilled(-(n * l));
    let emi = excel_expm1_internal(-(n * l));
    let lines: Vec<Vec<f64>> = rows(file);
    let keys: Vec<(&str, Box<dyn Fn(f64) -> u64>)> = vec![
        ("q=x87div(pv,em)", Box::new(move |pv| excel_x87_div(pv, emn).to_bits())),
        ("pv/em plain", Box::new(move |pv| (pv / emn).to_bits())),
        ("pv/emi", Box::new(move |pv| excel_x87_div(pv, emi).to_bits())),
        ("pv*r", Box::new(move |pv| excel_x87_mul(pv, r).to_bits())),
        ("pv*r plain", Box::new(move |pv| (pv * r).to_bits())),
        ("-pv", Box::new(move |pv| pv.to_bits())),
    ];
    for (name, key) in keys.iter() {
        let mut d: HashMap<u64, u64> = HashMap::new();
        let mut conflicts = 0;
        for p in &lines {
            let k = key(p[0]);
            let v = p[2].to_bits();
            if let Some(old) = d.get(&k) { if *old != v { conflicts += 1; } } else { d.insert(k, v); }
        }
        println!("{name:20} distinct {} conflicts {}", d.len(), conflicts);
    }
}


fn qdump(file: &str, r: f64, n: f64) {
    let l = alnrel(r);
    let emn = excel_expm1_spilled(-(n * l));
    for p in rows(file) {
        let q = excel_x87_div(p[0], emn);
        let qe = xdiv(&ld(p[0]), &ld(emn));
        println!("{:016x} {:016x} {:016x} {}", q.to_bits(), p[1].to_bits(), p[2].to_bits(), xbits(&qe));
    }
}


fn info(r: f64, n: f64) {
    let l = alnrel(r);
    let emn = excel_expm1_spilled(-(n * l));
    println!("L {:016x} emn {:016x} emi {:016x} u1 {:016x} em1 {:016x}", l.to_bits(), emn.to_bits(), excel_expm1_internal(-(n*l)).to_bits(), excel_exp(-l).to_bits(), excel_expm1_spilled(-l).to_bits());
}


fn qrext(file: &str) {
    report(file, &|p: &[f64]| {
        let (r, n, pv, k) = (p[0], p[1], p[2], p[3]);
        let m = n - k + 1.0;
        let l = alnrel(r);
        let emn = excel_expm1_spilled(-(n * l));
        let q = excel_x87_div(pv, emn);
        let qe = xdiv(&ld(pv), &ld(emn));
        let t = -(m * l);
        let u = excel_exp(t);
        let ue = xexp(&ld(t));
        let qr = xmul(&ld(q), &ld(r));
        let qer = xmul(&qe, &ld(r));
        vec![
            ("(q*r)x*u".into(), st(&xmul(&qr, &ld(u)))),
            ("(q*r)x*ue".into(), st(&xmul(&qr, &ue))),
            ("(qe*r)x*u".into(), st(&xmul(&qer, &ld(u)))),
            ("(qe*r)x*ue".into(), st(&xmul(&qer, &ue))),
            ("q*(r*u)x".into(), st(&xmul(&ld(q), &xmul(&ld(r), &ld(u))))),
            ("q*(r*ue)x".into(), st(&xmul(&ld(q), &xmul(&ld(r), &ue)))),
            ("(pv*r)x/em*u".into(), st(&xmul(&xdiv(&xmul(&ld(pv), &ld(r)), &ld(emn)), &ld(u)))),
            ("(pv*r*u)x/em".into(), st(&xdiv(&xmul(&xmul(&ld(pv), &ld(r)), &ld(u)), &ld(emn)))),
            ("(pv*u*r)x/em".into(), st(&xdiv(&xmul(&xmul(&ld(pv), &ld(u)), &ld(r)), &ld(emn)))),
            ("(pv*u)x/em*r".into(), st(&xmul(&xdiv(&xmul(&ld(pv), &ld(u)), &ld(emn)), &ld(r)))),
            ("(pv*ue)x/em*r".into(), st(&xmul(&xdiv(&xmul(&ld(pv), &ue), &ld(emn)), &ld(r)))),
            ("(pv*r*ue)x/em".into(), st(&xdiv(&xmul(&xmul(&ld(pv), &ld(r)), &ue), &ld(emn)))),
        ]
    });
}


fn type1(file: &str) {
    let mut tot = 0; let mut cnt = [0usize; 6]; let mut ulps = [0i64; 6];
    for p in rows(file) {
        let (r, n, pv, k, pmt1, want) = (p[0], p[1], p[2], p[3], p[6], p[7]);
        if k == 1.0 { continue; }
        let m = n - k + 1.0;
        let l = alnrel(r);
        let u = excel_exp(-(m * l));
        let emn = excel_expm1_spilled(-(n * l));
        let q = excel_x87_div(pv, emn);
        let due = excel_x87_recip(excel_x87_add(1.0, excel_x87_recip(r)));
        let cur = exp_portable(-(m * excel_log1p(r))) * pmt1;
        let c = [excel_x87_mul(pmt1, u), pmt1 * u, cur, excel_x87_mul(excel_x87_mul(q, u), due), excel_x87_mul(p[8], excel_exp(-((m + 1.0) * l))), excel_x87_div(excel_x87_mul(p[8], u), excel_x87_add(1.0, r))];
        tot += 1;
        for i in 0..6 { if c[i].to_bits() == want.to_bits() { cnt[i] += 1; } ulps[i] = ulps[i].max((c[i].to_bits() as i64 - want.to_bits() as i64).abs()); }
    }
    println!("rows {tot} B1 {} B1plain {} current {} Adue {} pmt0*u(m+1) {} pmt0*u/(1+r) {}", cnt[0], cnt[1], cnt[2], cnt[3], cnt[4], cnt[5]);
    println!("max ulps {:?}", ulps);
}


fn ranges(file: &str) {
    let mut names: Vec<String> = Vec::new();
    let mut per: std::collections::BTreeMap<String, Vec<usize>> = Default::default();
    for p in rows(file) {
        let (r, n, pv, s, e, pmt, want) = (p[0], p[1], p[2], p[3], p[4], p[6], p[7]);
        let l = alnrel(r);
        let emn = excel_expm1_spilled(-(n * l));
        let q = excel_x87_div(pv, emn);
        let u = |m: f64| excel_exp(-(m * l));
        let ue = |m: f64| xexp(&ld(-(m * l)));
        let terms: Vec<f64> = (s as i64..=e as i64).map(|k| excel_x87_mul(pmt, u(n - k as f64 + 1.0))).collect();
        let termsa: Vec<f64> = (s as i64..=e as i64).map(|k| excel_x87_mul(excel_x87_mul(q, u(n - k as f64 + 1.0)), r)).collect();
        let mut c: Vec<(String, f64)> = Vec::new();
        let mut a = 0.0; for t in &terms { a += t; } c.push(("B l2r".into(), a));
        let mut a = 0.0; for t in terms.iter().rev() { a += t; } c.push(("B r2l".into(), a));
        let mut a = ld(0.0); for t in &terms { a = xadd(&a, &ld(*t)); } c.push(("B x87acc".into(), st(&a)));
        let mut a = 0.0; for t in &termsa { a += t; } c.push(("A l2r".into(), a));
        let mut a = ld(0.0); for t in &termsa { a = xadd(&a, &ld(*t)); } c.push(("A x87acc".into(), st(&a)));
        // ext summand: pmt * ue, accumulate ext
        let mut a = ld(0.0); for k in s as i64..=e as i64 { a = xadd(&a, &xmul(&ld(pmt), &ld(u(n - k as f64 + 1.0)))); } c.push(("Bext acc".into(), st(&a)));
        let mut a = ld(0.0); for k in s as i64..=e as i64 { a = xadd(&a, &xmul(&ld(pmt), &ue(n - k as f64 + 1.0))); } c.push(("Bext ue acc".into(), st(&a)));
        // pmt * sum(u) : factor out
        let mut a = ld(0.0); for k in s as i64..=e as i64 { a = xadd(&a, &ld(u(n - k as f64 + 1.0))); } c.push(("pmt*sum(u)ext".into(), st(&xmul(&ld(pmt), &a))));
        let mut a = 0.0; for k in s as i64..=e as i64 { a += u(n - k as f64 + 1.0); } c.push(("pmt*sum(u)".into(), excel_x87_mul(pmt, a)));
        let mut a = ld(0.0); for k in s as i64..=e as i64 { a = xadd(&a, &ld(u(n - k as f64 + 1.0))); } c.push(("q*sum(u)*r".into(), excel_x87_mul(excel_x87_mul(q, st(&a)), r)));
        let ua = n - e; let ub = n - s + 1.0;
        let uu = |m: f64| if m == 0.0 { 1.0 } else { u(m) };
        c.push(("q*(ua-ub)".into(), excel_x87_mul(q, excel_x87_sub(uu(ua), uu(ub)))));
        c.push(("pmt*(ua-ub)/r".into(), excel_x87_div(excel_x87_mul(pmt, excel_x87_sub(uu(ua), uu(ub))), r)));
        if names.is_empty() { names = c.iter().map(|x| x.0.clone()).collect(); }
        let key = format!("{r}/{n}");
        let v = per.entry(key).or_insert(vec![0; c.len() + 1]);
        v[c.len()] += 1;
        for (i, (_, x)) in c.iter().enumerate() { if x.to_bits() == want.to_bits() { v[i] += 1; } }
    }
    print!("{:14}", "loan"); for nm in &names { print!(" {:>14}", nm); } println!(" rows");
    for (k, v) in per { print!("{:14}", k); for x in v { print!(" {:>14}", x); } println!(); }
}


fn recur(file: &str) {
    let mut names: Vec<String> = Vec::new();
    let mut per: std::collections::BTreeMap<String, Vec<usize>> = Default::default();
    for p in rows(file) {
        let (r, n, pv, s, e, pmt, want) = (p[0], p[1], p[2], p[3], p[4], p[6], p[7]);
        let l = alnrel(r);
        let g = excel_x87_add(1.0, r);
        let ms = n - s + 1.0;
        let p0 = excel_x87_mul(pmt, excel_exp(-(ms * l)));
        let p0x = xmul(&ld(pmt), &ld(excel_exp(-(ms * l))));
        let cnt = (e - s) as i64 + 1;
        let mut c: Vec<(String, f64)> = Vec::new();
        // R1 double P double sum
        { let mut pp = p0; let mut sum = 0.0; for _ in 0..cnt { sum += pp; pp = pp * g; } c.push(("Pd*g sumd".into(), sum)); }
        { let mut pp = p0; let mut sum = 0.0; for _ in 0..cnt { sum = excel_x87_add(sum, pp); pp = excel_x87_mul(pp, g); } c.push(("Pd*g x87 sumd".into(), sum)); }
        { let mut pp = ld(p0); let mut sum = ld(0.0); for _ in 0..cnt { sum = xadd(&sum, &pp); pp = xmul(&pp, &ld(g)); } c.push(("Px*g sumx".into(), st(&sum))); }
        { let mut pp = p0x; let mut sum = ld(0.0); for _ in 0..cnt { sum = xadd(&sum, &pp); pp = xmul(&pp, &ld(g)); } c.push(("P0x Px*g sumx".into(), st(&sum))); }
        { let mut pp = p0; let mut sum = ld(0.0); for _ in 0..cnt { sum = xadd(&sum, &ld(pp)); pp = excel_x87_mul(pp, g); } c.push(("Pd*g sumx".into(), st(&sum))); }
        { let mut pp = ld(p0); let mut sum = 0.0; for _ in 0..cnt { sum = excel_x87_add(sum, st(&pp)); pp = xmul(&pp, &ld(g)); } c.push(("Px*g sumd".into(), sum)); }
        // P += P*r form
        { let mut pp = p0; let mut sum = 0.0; for _ in 0..cnt { sum = excel_x87_add(sum, pp); pp = excel_x87_add(pp, excel_x87_mul(pp, r)); } c.push(("Pd+=Pd*r sumd".into(), sum)); }
        { let mut pp = ld(p0); let mut sum = ld(0.0); for _ in 0..cnt { sum = xadd(&sum, &pp); pp = xadd(&pp, &xmul(&pp, &ld(r))); } c.push(("Px+=Px*r sumx".into(), st(&sum))); }
        // t increment
        { let mut t = -(ms * l); let mut sum = 0.0; for _ in 0..cnt { sum = excel_x87_add(sum, excel_x87_mul(pmt, excel_exp(t))); t += l; } c.push(("t+=L sumd".into(), sum)); }
        { let mut t = -(ms * l); let mut sum = ld(0.0); for _ in 0..cnt { sum = xadd(&sum, &ld(excel_x87_mul(pmt, excel_exp(t)))); t += l; } c.push(("t+=L sumx".into(), st(&sum))); }
        // divide form P/v
        { let v = excel_x87_div(1.0, g); let mut pp = p0; let mut sum = 0.0; for _ in 0..cnt { sum = excel_x87_add(sum, pp); pp = excel_x87_div(pp, v); } c.push(("Pd/v sumd".into(), sum)); }
        // backward from end
        { let me = n - e + 1.0; let mut pp = excel_x87_mul(pmt, excel_exp(-(me * l))); let mut sum = 0.0; for _ in 0..cnt { sum = excel_x87_add(sum, pp); pp = excel_x87_div(pp, g); } c.push(("back Pd/g sumd".into(), sum)); }
        { let me = n - e + 1.0; let mut pp = ld(excel_x87_mul(pmt, excel_exp(-(me * l)))); let mut sum = ld(0.0); for _ in 0..cnt { sum = xadd(&sum, &pp); pp = xdiv(&pp, &ld(g)); } c.push(("back Px/g sumx".into(), st(&sum))); }
        let me = n - e + 1.0;
        let ue_ = excel_exp(-(me * l));
        { let mut pp = xmul(&ld(pmt), &ld(ue_)); let mut sum = ld(0.0); for _ in 0..cnt { sum = xadd(&sum, &pp); pp = xdiv(&pp, &ld(g)); } c.push(("back Pex/g sumx".into(), st(&sum))); }
        { let mut pp = ld(excel_x87_mul(pmt, ue_)); let mut sum = 0.0; for _ in 0..cnt { sum = excel_x87_add(sum, st(&pp)); pp = xdiv(&pp, &ld(g)); } c.push(("back Px/g sumd".into(), sum)); }
        { let v = xdiv(&xone(), &ld(g)); let mut pp = ld(excel_x87_mul(pmt, ue_)); let mut sum = ld(0.0); for _ in 0..cnt { sum = xadd(&sum, &pp); pp = xmul(&pp, &v); } c.push(("back Px*vx sumx".into(), st(&sum))); }
        { let v = ld(excel_x87_div(1.0, g)); let mut pp = ld(excel_x87_mul(pmt, ue_)); let mut sum = ld(0.0); for _ in 0..cnt { sum = xadd(&sum, &pp); pp = xmul(&pp, &v); } c.push(("back Px*vd sumx".into(), st(&sum))); }
        { let mut pp = ld(ue_); let mut sum = ld(0.0); for _ in 0..cnt { sum = xadd(&sum, &pp); pp = xdiv(&pp, &ld(g)); } c.push(("back u/g sumx *pmt".into(), st(&xmul(&sum, &ld(pmt))))); }
        { let mut pp = ld(ue_); let mut sum = ld(0.0); for _ in 0..cnt { sum = xadd(&sum, &pp); pp = xdiv(&pp, &ld(g)); } c.push(("back u/g sumd *pmt".into(), excel_x87_mul(st(&sum), pmt))); }
        { let mut pp = xexp(&ld(-(me * l))); let mut sum = ld(0.0); for _ in 0..cnt { sum = xadd(&sum, &pp); pp = xdiv(&pp, &ld(g)); } c.push(("back uex/g sumx *pmt".into(), st(&xmul(&sum, &ld(pmt))))); }
        if names.is_empty() { names = c.iter().map(|x| x.0.clone()).collect(); }
        let key = format!("{r}/{n}");
        let v = per.entry(key).or_insert(vec![0; c.len() + 1]);
        v[c.len()] += 1;
        for (i, (_, x)) in c.iter().enumerate() { if x.to_bits() == want.to_bits() { v[i] += 1; } }
    }
    for (i, nm) in names.iter().enumerate() { print!("{:16}", nm); for (_, v) in per.iter() { print!(" {:>4}", v[i]); } println!(); }
    print!("{:16}", "rows"); for (k, v) in per.iter() { print!(" {:>4}", v[names.len()]); } println!(); for k in per.keys() { print!("{k} "); } println!();
}


fn recur2(file: &str) {
    let mut names: Vec<String> = Vec::new();
    let mut per: std::collections::BTreeMap<String, Vec<usize>> = Default::default();
    for p in rows(file) {
        let (r, n, pv, s, e, pmt, want) = (p[0], p[1], p[2], p[3], p[4], p[6], p[7]);
        let l = alnrel(r);
        let g = excel_x87_add(1.0, r);
        let emn = excel_expm1_spilled(-(n * l));
        let q = excel_x87_div(pv, emn);
        let cnt = (e - s) as i64 + 1;
        let me = n - e + 1.0;
        let ue_ = excel_exp(-(me * l));
        let mut pp = ld(ue_); let mut sum = ld(0.0); for _ in 0..cnt { sum = xadd(&sum, &pp); pp = xdiv(&pp, &ld(g)); }
        let sd = st(&sum);
        let mut c: Vec<(String, f64)> = Vec::new();
        c.push(("S*pmt".into(), excel_x87_mul(sd, pmt)));
        c.push(("(q*S)*r".into(), excel_x87_mul(excel_x87_mul(q, sd), r)));
        c.push(("q*(S*r)".into(), excel_x87_mul(q, excel_x87_mul(sd, r))));
        c.push(("(S*r)*q".into(), excel_x87_mul(excel_x87_mul(sd, r), q)));
        c.push(("(q*r)x*S".into(), st(&xmul(&xmul(&ld(q), &ld(r)), &ld(sd)))));
        c.push(("(q*S)x*r".into(), st(&xmul(&xmul(&ld(q), &ld(sd)), &ld(r)))));
        c.push(("Sx*pmt".into(), st(&xmul(&sum, &ld(pmt)))));
        c.push(("(Sx*q)*r".into(), st(&xmul(&xmul(&sum, &ld(q)), &ld(r)))));
        // multiply-by-v variants with v double
        let v = excel_x87_div(1.0, g);
        let mut pp2 = ld(ue_); let mut sum2 = ld(0.0); for _ in 0..cnt { sum2 = xadd(&sum2, &pp2); pp2 = xmul(&pp2, &ld(v)); }
        c.push(("Sv*pmt".into(), excel_x87_mul(st(&sum2), pmt)));
        c.push(("(q*Sv)*r".into(), excel_x87_mul(excel_x87_mul(q, st(&sum2)), r)));
        // pp stored each step
        let mut pp3 = ue_; let mut sum3 = ld(0.0); for _ in 0..cnt { sum3 = xadd(&sum3, &ld(pp3)); pp3 = excel_x87_div(pp3, g); }
        c.push(("Sstored*pmt".into(), excel_x87_mul(st(&sum3), pmt)));
        c.push(("(q*Sst)*r".into(), excel_x87_mul(excel_x87_mul(q, st(&sum3)), r)));
        let qe = xdiv(&ld(pv), &ld(emn));
        c.push(("q*(Sd*r)x".into(), st(&xmul(&ld(q), &xmul(&ld(sd), &ld(r))))));
        c.push(("qe*(Sd*r)st".into(), st(&xmul(&qe, &ld(excel_x87_mul(sd, r))))));
        c.push(("qe*(Sd*r)x".into(), st(&xmul(&qe, &xmul(&ld(sd), &ld(r))))));
        c.push(("(qe*Sd)x*r".into(), st(&xmul(&xmul(&qe, &ld(sd)), &ld(r)))));
        c.push(("pv*(Sd*r)/em".into(), excel_x87_div(excel_x87_mul(pv, excel_x87_mul(sd, r)), emn)));
        c.push(("(Sd*r)*pv/em x".into(), st(&xdiv(&xmul(&ld(excel_x87_mul(sd, r)), &ld(pv)), &ld(emn)))));
        c.push(("(Sd*r)/em*pv".into(), excel_x87_mul(excel_x87_div(excel_x87_mul(sd, r), emn), pv)));
        c.push(("(Sd*r)x/em x*pv".into(), st(&xmul(&xdiv(&xmul(&ld(sd), &ld(r)), &ld(emn)), &ld(pv)))));
        let loopx = |f0: XExt| { let mut f = f0; let mut sum = ld(0.0); for _ in 0..cnt { sum = xadd(&sum, &f); f = xdiv(&f, &ld(g)); } sum };
        c.push(("V1 q*st(S(u*r))".into(), excel_x87_mul(q, st(&loopx(xmul(&ld(ue_), &ld(r)))))));
        c.push(("V1b q*st(S(st(u*r)))".into(), excel_x87_mul(q, st(&loopx(ld(excel_x87_mul(ue_, r)))))));
        let ppmt_e = excel_x87_mul(excel_x87_mul(q, ue_), r);
        c.push(("V2 S(PPMT_e)".into(), st(&loopx(ld(ppmt_e)))));
        c.push(("V4 S(qx*u*r)".into(), st(&loopx(xmul(&xmul(&ld(q), &ld(ue_)), &ld(r))))));
        c.push(("V5 S(q*u st)*r".into(), excel_x87_mul(st(&loopx(ld(excel_x87_mul(q, ue_)))), r)));
        c.push(("V6 S((q*u)x)*r".into(), excel_x87_mul(st(&loopx(xmul(&ld(q), &ld(ue_)))), r)));
        c.push(("V7 S(q*u st) x r".into(), st(&xmul(&loopx(ld(excel_x87_mul(q, ue_))), &ld(r)))));
        if names.is_empty() { names = c.iter().map(|x| x.0.clone()).collect(); }
        let key = format!("{r}/{n}");
        let v = per.entry(key).or_insert(vec![0; c.len() + 1]);
        v[c.len()] += 1;
        for (i, (_, x)) in c.iter().enumerate() { if x.to_bits() == want.to_bits() { v[i] += 1; } }
    }
    for (i, nm) in names.iter().enumerate() { print!("{:16}", nm); for (_, v) in per.iter() { print!(" {:>4}", v[i]); } println!(); }
    print!("{:16}", "rows"); for (_, v) in per.iter() { print!(" {:>4}", v[names.len()]); } println!(); for k in per.keys() { print!("{k} "); } println!();
}


fn miss3(file: &str) {
    for p in rows(file) {
        let (r, n, pv, s, e, pmt, want) = (p[0], p[1], p[2], p[3], p[4], p[6], p[7]);
        if r != 3.0 { continue; }
        let l = alnrel(r);
        let g = excel_x87_add(1.0, r);
        let emn = excel_expm1_spilled(-(n * l));
        let q = excel_x87_div(pv, emn);
        let cnt = (e - s) as i64 + 1;
        let me = n - e + 1.0;
        let ue_ = excel_exp(-(me * l));
        let mut pp = ld(ue_); let mut sum = ld(0.0); for _ in 0..cnt { sum = xadd(&sum, &pp); pp = xdiv(&pp, &ld(g)); }
        let sd = st(&sum);
        let mut sdd = 0.0; let mut ppd = ue_; for _ in 0..cnt { sdd = excel_x87_add(sdd, ppd); ppd = ppd / g; }
        let a = excel_x87_mul(q, excel_x87_mul(sd, r));
        let b = excel_x87_mul(q, excel_x87_mul(sdd, r));
        let d = |x: f64| x.to_bits() as i64 - want.to_bits() as i64;
        if a != want { println!("s={s} e={e} cnt={cnt} H{:+} Hdblsum{:+} S={:016x} Sdd={:016x}", d(a), d(b), sd.to_bits(), sdd.to_bits()); }
    }
}


fn lscan2(file: &str) {
    for g in groups(file) {
        let (r, n, pv) = (g[0][0], g[0][1], g[0][2]);
        let base = alnrel(r);
        let emn0 = excel_expm1_spilled(-(n * base));
        let q0 = excel_x87_div(pv, emn0);
        println!("loan r={r} n={n} pv={pv} rows={}", g.len());
        for off in -8i64..=8 {
            let l = f64::from_bits((base.to_bits() as i64 + off) as u64);
            let emn = excel_expm1_spilled(-(n * l));
            let q = excel_x87_div(pv, emn);
            let mut cnt = [0usize; 4];
            for p in &g {
                let want = p[7];
                let m = n - p[3] + 1.0;
                let u = excel_exp(-(m * l));
                let c = [excel_x87_mul(q0, excel_x87_mul(u, r)), excel_x87_mul(q, excel_x87_mul(u, r)),
                         excel_x87_mul(q0, st(&xmul(&xexp(&xchs(&xmul(&ld(m), &ld(l)), XCW)), &ld(r)))),
                         excel_x87_mul(q0, st(&xmul(&xexp(&ld(-(m * l))), &ld(r))))];
                for i in 0..4 { if c[i].to_bits() == want.to_bits() { cnt[i] += 1; } }
            }
            println!("  off {off:3} C(q0) {:4} C(q') {:4} Cext-t {:4} Cue {:4}", cnt[0], cnt[1], cnt[2], cnt[3]);
        }
    }
}


fn biases() {
    for r in [0.25f64, 0.5, 0.625, 0.75, 0.875, 1.0, 1.5, 3.0, 0.1, 0.3, 2.0] {
        let l = alnrel(r);
        let u1 = excel_exp(-l);
        let ur = excel_x87_mul(u1, r);
        let em1 = excel_expm1_spilled(-l);
        println!("{r} u1 {:016x} ur {:016x} em1 {:016x} lnx {:016x} l {:016x}", u1.to_bits(), ur.to_bits(), em1.to_bits(), excel_log(1.0+r).to_bits(), l.to_bits());
    }
}


pub fn v1(r: f64, n: f64, pv: f64, s: f64, e: f64) -> f64 {
    let l = alnrel(r);
    let g = excel_x87_add(1.0, r);
    let emn = excel_expm1_spilled(-(n * l));
    let q = excel_x87_div(pv, emn);
    let cnt = (e - s) as i64 + 1;
    let me = n - e + 1.0;
    let mut f = xmul(&ld(excel_exp(-(me * l))), &ld(r));
    let mut sum = ld(0.0);
    let gx = xadd(&xone(), &ld(r));
    let _ = g;
    for _ in 0..cnt { sum = xadd(&sum, &f); f = xdiv(&f, &gx); }
    excel_x87_mul(q, st(&sum))
}
fn cur(r: f64, n: f64, pv: f64, s: f64, e: f64) -> f64 {
    let l = alnrel(r);
    let emn = excel_expm1_spilled(-(n * l));
    let q = excel_x87_div(pv, emn);
    let mut sum = 0.0;
    for k in s as i64..=e as i64 { sum += excel_x87_mul(excel_x87_mul(q, excel_exp(-((n - k as f64 + 1.0) * l))), r); }
    sum
}
fn cmp(file: &str) {
    let mut st_: std::collections::BTreeMap<String, [i64; 6]> = Default::default();
    for p in rows(file) {
        let (r, n, pv, s, e, want) = (p[0], p[1], p[2], p[3], p[4], p[p.len()-1]);
        let a = v1(r, n, pv, s, e); let b = cur(r, n, pv, s, e);
        let da = (a.to_bits() as i64 - want.to_bits() as i64).abs(); let db = (b.to_bits() as i64 - want.to_bits() as i64).abs();
        let x = st_.entry(format!("{r}/{n}")).or_insert([0; 6]);
        x[0] += 1; x[1] += (da == 0) as i64; x[2] += (db == 0) as i64; x[3] = x[3].max(da); x[4] = x[4].max(db); x[5] += (da <= db) as i64;
    }
    println!("loan rows v1exact curexact v1max curmax v1<=cur");
    for (k, v) in st_ { println!("{k} {:?}", v); }
}


fn expvar(file: &str) {
    report(file, &|p: &[f64]| {
        let (r, n, pv, s, e) = (p[0], p[1], p[2], p[3], p[4]);
        let l = alnrel(r);
        let emn = excel_expm1_spilled(-(n * l));
        let q = excel_x87_div(pv, emn);
        let cnt = (e - s) as i64 + 1;
        let me = n - e + 1.0;
        let t = -(me * l);
        let gx = xadd(&xone(), &ld(r));
        let run = |u0: XExt| { let mut f = xmul(&u0, &ld(r)); let mut sum = ld(0.0); for _ in 0..cnt { sum = xadd(&sum, &f); f = xdiv(&f, &gx); } excel_x87_mul(q, st(&sum)) };
        let run2 = |u0: XExt| { let mut f = u0; let mut sum = ld(0.0); for _ in 0..cnt { sum = xadd(&sum, &f); f = xdiv(&f, &gx); } excel_x87_mul(q, excel_x87_mul(st(&sum), r)) };
        vec![
            ("excel_exp".into(), run(ld(excel_exp(t)))),
            ("exp ext".into(), run(xexp(&ld(t)))),
            ("crt_exp".into(), run(ld(crt_exp(t)))),
            ("exp_portable".into(), run(ld(exp_portable(t)))),
            ("crt_pow(g,-m)".into(), run(ld(crt_pow(1.0 + r, -me)))),
            ("x87 pow direct".into(), run(ld(excel_pow_x87_direct(1.0 + r, -me)))),
            ("1/crt_pow".into(), run(ld(1.0 / crt_pow(1.0 + r, me)))),
            ("run2 excel_exp".into(), run2(ld(excel_exp(t)))),
            ("run2 crt_exp".into(), run2(ld(crt_exp(t)))),
            ("crt_exp(crtlog)".into(), run(ld(crt_exp(-(me * crt_log(1.0 + r)))))),
        ]
    });
}


fn pc53(file: &str) {
    report(file, &|p: &[f64]| {
        let (r, n, pv, s, e) = (p[0], p[1], p[2], p[3], p[4]);
        let mut v: Vec<(String, f64)> = Vec::new();
        for (cwn, cw) in [("pc64", 0x133Fu16), ("pc53", 0x123Fu16)] {
            let l = alnrel(r);
            let emn = excel_expm1_spilled(-(n * l));
            let q = excel_x87_div(pv, emn);
            let cnt = (e - s) as i64 + 1;
            let me = n - e + 1.0;
            let t = -(me * l);
            let gx = xaddc(&xone(), &ld(r), cw);
            let u0 = xexpc(&ld(t), cw);
            let mut f = xmulc(&u0, &ld(r), cw); let mut sum = ld(0.0);
            for _ in 0..cnt { sum = xaddc(&sum, &f, cw); f = xdivc(&f, &gx, cw); }
            v.push((format!("{cwn} q*S"), xstorec(&xmulc(&ld(q), &ld(xstorec(&sum, cw)), cw), cw)));
            v.push((format!("{cwn} q*Sx"), xstorec(&xmulc(&ld(q), &sum, cw), cw)));
            // u stored to double first
            let u0d = ld(xstorec(&u0, cw));
            let mut f = xmulc(&u0d, &ld(r), cw); let mut sum = ld(0.0);
            for _ in 0..cnt { sum = xaddc(&sum, &f, cw); f = xdivc(&f, &gx, cw); }
            v.push((format!("{cwn} ud q*S"), xstorec(&xmulc(&ld(q), &ld(xstorec(&sum, cw)), cw), cw)));
            // L in this cw too
            let lx = xlogc(&ld(1.0 + r), cw);
            let l2 = xstorec(&lx, cw);
            let t2 = -(me * l2);
            let u2 = ld(xstorec(&xexpc(&ld(t2), cw), cw));
            let mut f = xmulc(&u2, &ld(r), cw); let mut sum = ld(0.0);
            for _ in 0..cnt { sum = xaddc(&sum, &f, cw); f = xdivc(&f, &gx, cw); }
            v.push((format!("{cwn} lnL ud q*S"), xstorec(&xmulc(&ld(q), &ld(xstorec(&sum, cw)), cw), cw)));
        }
        v
    });
}


fn bank0(file: &str) {
    let (mut tot, mut a, mut b, mut lose, mut gain) = (0, 0, 0, 0, 0);
    let mut maxa = 0i64; let mut maxb = 0i64;
    for p in rows(file) {
        if p[5] != 0.0 { continue; }
        let (r, n, pv, s, e, want) = (p[0], p[1].trunc(), p[2], p[3].trunc(), p[4].trunc(), p[7]);
        if !(r > 0.0) { continue; }
        let x = v1(r, n, pv, s, e); let y = cur(r, n, pv, s, e);
        tot += 1;
        let xa = x.to_bits() == want.to_bits(); let yb = y.to_bits() == want.to_bits();
        a += xa as i32; b += yb as i32; if yb && !xa { lose += 1; } if xa && !yb { gain += 1; }
        maxa = maxa.max((x.to_bits() as i64 - want.to_bits() as i64).abs()); maxb = maxb.max((y.to_bits() as i64 - want.to_bits() as i64).abs());
    }
    println!("type0 rows {tot} v1 {a} cur {b} v1-gains {gain} v1-losses {lose} maxulp v1 {maxa} cur {maxb}");
}


fn pmt1(r: f64, n: f64, pv: f64) -> f64 {
    let l = alnrel(r);
    let emn = excel_expm1_spilled(-(n * l));
    let q = excel_x87_div(pv, emn);
    excel_x87_mul(q, due(r))
}
fn due(r: f64) -> f64 { let rc = excel_x87_recip(r); if rc.is_finite() { excel_x87_recip(excel_x87_add(1.0, rc)) } else { r } }
fn cur1(r: f64, n: f64, pv: f64, s: f64, e: f64) -> f64 {
    let pay = pmt1(r, n, pv);
    let mut sum = 0.0;
    for k in s as i64..=e as i64 {
        let p = if k == 1 { pay } else {
            let m = n - k as f64 + 1.0;
            exp_portable(-(m * excel_log1p(r))) * pay
        };
        sum += p;
    }
    sum
}
fn t1forms(r: f64, n: f64, pv: f64, s: f64, e: f64) -> Vec<f64> {
    let l = alnrel(r);
    let emn = excel_expm1_spilled(-(n * l));
    let q = excel_x87_div(pv, emn);
    let pay = pmt1(r, n, pv);
    let d = due(r);
    let mut out = Vec::new();
    // F0: sum of B1 terms (pmt1*u), double, k=1 term = pay
    let mut s0 = 0.0; let mut s1 = 0.0;
    for k in s as i64..=e as i64 {
        let m = n - k as f64 + 1.0;
        let u = excel_exp(-(m * l));
        s0 += if k == 1 { pay } else { excel_x87_mul(pay, u) };
        s1 += if k == 1 { pay } else { excel_x87_mul(excel_x87_mul(q, u), d) };
    }
    out.push(s0); out.push(s1);
    // V1-like for type 1: loop over k=e..max(s,2), f = u_e*due ext, /gx; plus k=1 term pay if s==1
    let gx = xadd(&xone(), &ld(r));
    let s2 = if s < 2.0 { 2.0 } else { s };
    let mut tot = 0.0;
    if e >= s2 {
        let cnt = (e - s2) as i64 + 1;
        let me = n - e + 1.0;
        let mut f = xmul(&ld(excel_exp(-(me * l))), &ld(d)); let mut sum = ld(0.0);
        for _ in 0..cnt { sum = xadd(&sum, &f); f = xdiv(&f, &gx); }
        tot = excel_x87_mul(q, st(&sum));
    }
    out.push(if s < 2.0 { excel_x87_add(tot, pay) } else { tot });
    out.push(if s < 2.0 { excel_x87_add(pay, tot) } else { tot });
    // V1 with r factor then due: q*st(S(u*r))/(1+r)?
    let mut tot2 = 0.0;
    if e >= s2 {
        let cnt = (e - s2) as i64 + 1;
        let me = n - e + 1.0;
        let mut f = xmul(&ld(excel_exp(-(me * l))), &ld(r)); let mut sum = ld(0.0);
        for _ in 0..cnt { sum = xadd(&sum, &f); f = xdiv(&f, &gx); }
        tot2 = excel_x87_div(excel_x87_mul(q, st(&sum)), excel_x87_add(1.0, r));
    }
    out.push(if s < 2.0 { excel_x87_add(tot2, pay) } else { tot2 });
    out
}
fn bank1(file: &str) {
    let names = ["sumB1", "sumAdue", "V1due tot+pay", "V1due pay+tot", "V1r/(1+r)"];
    let mut tot = 0; let mut c = 0; let mut k = [0usize; 5]; let mut lose = [0usize; 5]; let mut mx = [0i64; 6];
    for p in rows(file) {
        if p[5] != 1.0 { continue; }
        let (r, n, pv, s, e, want) = (p[0], p[1].trunc(), p[2], p[3].trunc(), p[4].trunc(), p[p.len()-1]);
        if !(r > 0.0) { continue; }
        tot += 1;
        let y = cur1(r, n, pv, s, e); let yb = y.to_bits() == want.to_bits(); c += yb as usize;
        mx[5] = mx[5].max((y.to_bits() as i64 - want.to_bits() as i64).abs());
        for (i, x) in t1forms(r, n, pv, s, e).iter().enumerate() {
            let xa = x.to_bits() == want.to_bits(); k[i] += xa as usize; if yb && !xa { lose[i] += 1; }
            mx[i] = mx[i].max((x.to_bits() as i64 - want.to_bits() as i64).abs());
        }
    }
    println!("type1 rows {tot} cur {c} maxulp cur {}", mx[5]);
    for i in 0..5 { println!("  {:16} {} losses {} max {}", names[i], k[i], lose[i], mx[i]); }
}


fn uimp(file: &str, rsel: f64) {
    for p in rows(file) {
        let (r, n, pv, k, want) = (p[0], p[1], p[2], p[3], p[7]);
        if r != rsel { continue; }
        let l = alnrel(r);
        let emn = excel_expm1_spilled(-(n * l));
        let q = excel_x87_div(pv, emn);
        let m = n - k + 1.0;
        let t = -(m * l);
        let u = excel_exp(t);
        println!("{} {:016x} {:016x} {:016x} {} {:016x}", m, q.to_bits(), u.to_bits(), want.to_bits(), xbits(&xexp(&ld(t))), t.to_bits());
    }
}


fn v1q(q: f64, r: f64, n: f64, s: f64, e: f64) -> f64 {
    let l = alnrel(r);
    let cnt = (e - s) as i64 + 1;
    let me = n - e + 1.0;
    let gx = xadd(&xone(), &ld(r));
    let mut f = xmul(&ld(excel_exp(-(me * l))), &ld(r));
    let mut sum = ld(0.0);
    for _ in 0..cnt { sum = xadd(&sum, &f); f = xdiv(&f, &gx); }
    excel_x87_mul(q, st(&sum))
}
fn fwd(file: &str) {
    report(file, &|p: &[f64]| {
        let (r, n, pv, s, e) = (p[0], p[1], p[2], p[3], p[4]);
        let l = alnrel(r);
        let gn = excel_exp(n * l);
        let emn = excel_expm1_spilled(-(n * l));
        let ep = excel_expm1_spilled(n * l);
        let q0 = excel_x87_div(pv, emn);
        let qs: Vec<(&str, f64)> = vec![
            ("q=pv/em", q0),
            ("-(pv*G)/(G-1)", -excel_x87_div(excel_x87_mul(pv, gn), excel_x87_sub(gn, 1.0))),
            ("-(pv*G)/ep", -excel_x87_div(excel_x87_mul(pv, gn), ep)),
            ("-pv/(1-1/G)", -excel_x87_div(pv, excel_x87_sub(1.0, excel_x87_div(1.0, gn)))),
            ("-(pv*G)x/(G-1)", -st(&xdiv(&xmul(&ld(pv), &ld(gn)), &xsub(&ld(gn), &xone())))),
            ("-pv*G/ep x", -st(&xdiv(&xmul(&ld(pv), &ld(gn)), &ld(ep)))),
            ("-pv/(ep/G)", -excel_x87_div(pv, excel_x87_div(ep, gn))),
            ("pv/(-ep/G) x", -st(&xdiv(&ld(pv), &xdiv(&ld(ep), &ld(gn))))),
        ];
        qs.iter().map(|(nm, q)| (nm.to_string(), v1q(*q, r, n, s, e))).collect()
    });
}


fn toff(file: &str) {
    for p in rows(file) {
        let (r, n, pv, k, want) = (p[0], p[1], p[2], p[3], p[7]);
        let l = alnrel(r);
        let emn = excel_expm1_spilled(-(n * l));
        let q = excel_x87_div(pv, emn);
        let m = n - k + 1.0;
        let t0 = -(m * l);
        let mut hits = String::new();
        for j in -4i64..=4 {
            let t = f64::from_bits((t0.to_bits() as i64 + j) as u64);
            let u = excel_exp(t);
            let pr = excel_x87_mul(q, excel_x87_mul(u, r));
            if pr == want { hits.push_str(&format!("{j:+} ")); }
        }
        println!("r={r} m={m} hits {hits}");
    }
}


fn binv(file: &str) {
    let rows = rows(file);
    let m0s: Vec<(&str, Box<dyn Fn(f64, f64) -> f64>)> = vec![
        ("crtpow", Box::new(|q: f64, n: f64| crt_pow(q, n))),
        ("exp(n*log1p(-p))", Box::new(|q: f64, n: f64| excel_exp(n * alnrel(q - 1.0)))),
        ("exp(n*ln q)", Box::new(|q: f64, n: f64| excel_exp(n * excel_log(q)))),
        ("powchain", Box::new(|q: f64, n: f64| excel_pow_chain(q, n))),
        ("x87direct", Box::new(|q: f64, n: f64| excel_pow_x87_direct(q, n))),
        ("repeated", Box::new(|q: f64, n: f64| { let mut r = 1.0; for _ in 0..(n as i64) { r *= q; } r })),
    ];
    for (mn, m0) in m0s.iter() {
        for fac in 0..4 {
            for strict in [false, true] {
                let mut ok = 0; let mut tot = 0;
                for p in &rows {
                    let (n, pp, alpha, want) = (p[0].trunc(), p[1], p[2], p[3]);
                    if n > 1000.0 { continue; }
                    let q = 1.0 - pp;
                    tot += 1;
                    let mut mass = m0(q, n);
                    let mut sum = 0.0; let mut k = 0.0; let mut res = n;
                    let ratio = pp / q;
                    while k < n {
                        sum += mass;
                        if (strict && sum > alpha) || (!strict && sum >= alpha) { res = k; break; }
                        k += 1.0;
                        mass = match fac { 0 => mass * ((n - k + 1.0) / k * pp / q), 1 => mass * (n - k + 1.0) / k * pp / q, 2 => mass * ((n - k + 1.0) / k) * ratio, _ => mass * (n - k + 1.0) * pp / (k * q) };
                    }
                    if res == want { ok += 1; }
                }
                println!("{mn:18} fac{fac} strict={strict} {ok}/{tot}");
            }
        }
    }
}

pub fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args[1].as_str() {
        "singles" => report(&args[2], &singles),
        "lscan" => lscan(&args[2]),
        "resid" => resid(&args[2]),
        "dumpu" => dumpu(&args[2]),
        "dumpl" => dumpl(&args[2]),
        "rz" => rz(&args[2]),
        "diffs" => diffs(&args[2]),
        "cscan" => cscan(&args[2]),
        "crt" => crt(&args[2]),
        "show" => show(&args[2]),
        "lext" => lext(&args[2]),
        "binv" => binv(&args[2]),
        "toff" => toff(&args[2]),
        "fwd" => fwd(&args[2]),
        "uimp" => uimp(&args[2], args[3].parse().unwrap()),
        "bank1" => bank1(&args[2]),
        "bank0" => bank0(&args[2]),
        "pc53" => pc53(&args[2]),
        "expvar" => expvar(&args[2]),
        "cmp" => cmp(&args[2]),
        "biases" => biases(),
        "lscan2" => lscan2(&args[2]),
        "miss3" => miss3(&args[2]),
        "recur2" => recur2(&args[2]),
        "recur" => recur(&args[2]),
        "ranges" => ranges(&args[2]),
        "type1" => type1(&args[2]),
        "qrext" => qrext(&args[2]),
        "info" => info(args[2].parse().unwrap(), args[3].parse().unwrap()),
        "qdump" => qdump(&args[2], args[3].parse().unwrap(), args[4].parse().unwrap()),
        "funcof" => funcof(&args[2], args[3].parse().unwrap(), args[4].parse().unwrap()),
        _ => panic!("mode"),
    }
}
