//! First contradiction elimination, followed by whole-family reachability.
use crate::ported::excel_numeric::Fpu;
use crate::{Row, dsl::*, enumerate, eval::*, json_write};
use serde_json::{Value, json};
use std::{
    fs,
    io::{BufWriter, Write},
    path::Path,
    sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    },
    time::Instant,
};

pub fn run<F: Fpu + 'static>(rows: &[Row], output: &str) {
    fs::create_dir_all(output).unwrap();
    let out = Path::new(output);
    let candidates = enumerate::family();
    let mut programs = BufWriter::new(fs::File::create(out.join("programs.jsonl")).unwrap());
    for c in &candidates {
        serde_json::to_writer(&mut programs, c).unwrap();
        writeln!(programs).unwrap();
    }
    programs.flush().unwrap();
    let grammar = hash_bytes(&serde_json::to_vec(&candidates).unwrap());
    let mut alive: Vec<_> = (0..candidates.len()).collect();
    let mut kills = BufWriter::new(fs::File::create(out.join("first-kills.jsonl")).unwrap());
    let mut statistics = BufWriter::new(fs::File::create(out.join("witness-kills.jsonl")).unwrap());
    let mut kill_summary = Vec::new();
    let mut total_evaluations = 0usize;
    for (order, r) in rows.iter().enumerate() {
        let mut next = Vec::new();
        let before = alive.len();
        for id in alive {
            let c = &candidates[id];
            let got = outcome::<F>(eval_graph::<F>(&c.graph, &r.args()));
            total_evaluations += 1;
            if got == r.expected() {
                next.push(id);
            } else {
                serde_json::to_writer(
                    &mut kills,
                    &json!({"candidate_id":id,"candidate_hash":c.hash(),"order":order,
                    "trace_index":r.trace_index,"witness":r.raw,"expected":r.expected(),"got":got}),
                )
                .unwrap();
                writeln!(kills).unwrap();
            }
        }
        let killed = before - next.len();
        alive = next;
        let stats = json!({"order":order,"trace_index":r.trace_index,"id":r.raw["id"],"source":r.raw["source"],
            "reading":r.raw["reading"],"killed":killed,"remaining":alive.len(),"prior_gap":r.prior_gap});
        serde_json::to_writer(&mut statistics, &stats).unwrap();
        writeln!(statistics).unwrap();
        if killed > 0 {
            println!("kill {stats}");
            kill_summary.push(stats);
        }
        if order < 32 || killed > 0 {
            kills.flush().unwrap();
            statistics.flush().unwrap();
            json_write(
                out.join("elimination-checkpoint.json"),
                &json!({"next_order":order+1,"alive_ids":alive,"grammar_hash":grammar,"complete":order+1==rows.len()}),
            );
        }
    }
    let survivors: Vec<_> = alive.iter().map(|id| &candidates[*id]).collect();
    json_write(out.join("survivors.json"), &survivors);
    json_write(
        out.join("elimination-checkpoint.json"),
        &json!({"next_order":rows.len(),"alive_ids":alive,"grammar_hash":grammar,"complete":true}),
    );
    println!(
        "elimination complete: {} programs, {} survivors, {total_evaluations} exact evaluations",
        candidates.len(),
        survivors.len()
    );
    let start = Instant::now();
    let next = AtomicUsize::new(0);
    let progress = AtomicUsize::new(0);
    let workers = 8usize.min(rows.len().max(1));
    let counts: Vec<_> = std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for worker in 0..workers {
            let next = &next;
            let progress = &progress;
            handles.push(scope.spawn(move || {
                let mut f=BufWriter::new(fs::File::create(out.join(format!("reachability-{worker}.jsonl"))).unwrap());
                let mut certificates=BufWriter::new(fs::File::create(out.join(format!("unreachable-roots-{worker}.jsonl"))).unwrap());
                let mut unreachable=0;let mut prior_gap_unreachable=0;let mut processed=0;let mut outcomes=0;
                loop {
                    let order=next.fetch_add(1,Ordering::Relaxed);if order>=rows.len(){break}
                    let r=&rows[order];let roots=root_values::<F>(&r.args());
                    let outputs=publish_roots::<F>(&roots);outcomes+=outputs.len();
                    let reachable=outputs.iter().any(|b|bits_outcome(*b)==r.expected());
                    unreachable+=(!reachable) as usize;prior_gap_unreachable+=(!reachable && r.prior_gap) as usize;
                    let record=json!({"order":order,"trace_index":r.trace_index,"id":r.raw["id"],"source":r.raw["source"],
                        "reachable":reachable,"prior_gap":r.prior_gap,"expected":r.expected(),
                        "outputs":outputs.iter().map(|b|bits_outcome(*b)).collect::<Vec<_>>()});
                    serde_json::to_writer(&mut f,&record).unwrap();writeln!(f).unwrap();processed+=1;
                    if !reachable {
                        serde_json::to_writer(&mut certificates,&json!({"order":order,"trace_index":r.trace_index,
                            "expected":r.expected(),"raw_roots":roots.iter().map(|v|v.exact_encoding()).collect::<Vec<_>>() })).unwrap();
                        writeln!(certificates).unwrap();
                    }
                    let done=progress.fetch_add(1,Ordering::Relaxed)+1;
                    if done%10_000==0 {f.flush().unwrap();println!("reachability {done}/{} {:.1}s",rows.len(),start.elapsed().as_secs_f64());}
                }
                (processed,unreachable,prior_gap_unreachable,outcomes)
            }));
        }
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let totals = counts.iter().fold((0, 0, 0, 0), |a, b| {
        (a.0 + b.0, a.1 + b.1, a.2 + b.2, a.3 + b.3)
    });
    assert_eq!(totals.0, rows.len());
    let summary = json!({"grammar_hash":grammar,"programs":candidates.len(),"scoped_recorded_rows":rows.len(),
        "survivors":survivors.len(),"exact_elimination_evaluations":total_evaluations,"nonzero_kill_witnesses":kill_summary,
        "whole_family_reachable":rows.len()-totals.1,"whole_family_unreachable":totals.1,
        "prior_gap_unreachable":totals.2,"published_outcomes_enumerated":totals.3,
        "reachability_seconds":start.elapsed().as_secs_f64(),"workers":workers,"complete":true,
        "selection":"first exact contradiction; no tolerance, rankings, coefficients or fitted switches",
        "backend":std::any::type_name::<F>(),"fresh_validation_rows":0});
    json_write(out.join("summary.json"), &summary);
    println!("{summary}");
}

/// One declared structural extension. This streams first contradictions and
/// does not build a prediction matrix or rank contradicted candidates.
pub fn associations<F: Fpu>(rows: &[Row], output: &str) {
    association_programs::<F>(
        rows,
        output,
        enumerate::association_family(),
        "4 numerators x 8 combines x 3^8 plus one 7-site rearrangement x 3^7. All annotations; root duplicates retained. No claim of all possible PMT programs.",
    );
}

pub fn timing_corrections<F: Fpu>(rows: &[Row], output: &str) {
    association_programs::<F>(
        rows,
        output,
        enumerate::timing_family(),
        "4 numerators x 3 subtractive timing corrections x 3^10. All arithmetic annotations; exact negations remain resident PC64. Frozen binary64 helper, not a claim about other helpers or all PMT programs.",
    );
}

fn association_programs<F: Fpu>(
    rows: &[Row],
    output: &str,
    candidates: Vec<enumerate::AssociationCandidate>,
    coverage: &str,
) {
    let out = Path::new(output);
    fs::create_dir(out).unwrap(); // never overwrite a previous experiment
    let mut programs = BufWriter::new(
        fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(out.join("programs.jsonl"))
            .unwrap(),
    );
    let mut buckets = std::collections::HashMap::<String, Vec<&Graph>>::new();
    for c in &candidates {
        let hash = hash_bytes(&serde_json::to_vec(&c.graph).unwrap());
        let bucket = buckets.entry(hash).or_default();
        assert!(!bucket.contains(&&c.graph), "duplicate serialized graph");
        bucket.push(&c.graph);
        serde_json::to_writer(&mut programs, c).unwrap();
        writeln!(programs).unwrap();
    }
    programs.flush().unwrap();
    let grammar = hash_bytes(&serde_json::to_vec(&candidates).unwrap());
    let mut alive: Vec<_> = (0..candidates.len()).collect();
    let mut kills = BufWriter::new(
        fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(out.join("first-kills.jsonl"))
            .unwrap(),
    );
    let mut kill_summary = Vec::new();
    let mut total_evaluations = 0usize;
    let mut visited = 0usize;
    for (order, row) in rows.iter().enumerate() {
        if alive.is_empty() {
            break;
        }
        let before = alive.len();
        let mut next = Vec::new();
        for id in alive {
            let c = &candidates[id];
            let root = eval_graph::<F>(&c.graph, &row.args());
            let got = outcome::<F>(root);
            total_evaluations += 1;
            if got == row.expected() {
                next.push(id);
            } else {
                // Each first kill retains every operation's full-width value.
                let trace = crate::eval::trace_graph::<F>(&c.graph, &row.args());
                serde_json::to_writer(
                    &mut kills,
                    &json!({"candidate_id":id,"order":order,
                    "trace_index":row.trace_index,"witness":row.raw,"expected":row.expected(),
                    "got":got,"root":root.exact_encoding(),"trace":trace}),
                )
                .unwrap();
                writeln!(kills).unwrap();
            }
        }
        alive = next;
        visited = order + 1;
        if before != alive.len() {
            let stats = json!({"order":order,"trace_index":row.trace_index,"id":row.raw["id"],
                "killed":before-alive.len(),"remaining":alive.len()});
            println!("association kill {stats}");
            kill_summary.push(stats);
        }
        kills.flush().unwrap();
        json_write(
            out.join("checkpoint.json"),
            &json!({"next_order":visited,
            "alive_ids":alive,"grammar_hash":grammar,"complete":alive.is_empty() || visited==rows.len()}),
        );
    }
    kills.flush().unwrap();
    let survivors: Vec<_> = alive.iter().map(|id| &candidates[*id]).collect();
    json_write(out.join("survivors.json"), &survivors);
    let summary = json!({"complete":true,"programs":candidates.len(),"survivors":survivors.len(),
        "grammar_hash":grammar,"scoped_rows_available":rows.len(),"rows_visited":visited,
        "exact_evaluations":total_evaluations,"first_kill_counts":kill_summary,
        "backend":std::any::type_name::<F>(),"fresh_validation_rows":0,
        "scope":"Original finite discount scope; frozen binary64 log1p/expm1 are conditional hypotheses",
        "coverage":coverage,
        "selection":"First exact contradiction; no match-count ranking, closeness beam, fitted constants or guards",
        "questions":"Only nonempty survivor sets justify distinguishing questions; requests would be frozen separately from predictions"});
    json_write(out.join("summary.json"), &summary);
    println!("{summary}");
}
