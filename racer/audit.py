"""Reconcile both backend runs, scoped identities, first kills and reachability.

All Excel input/answer/provenance fields come from the existing witness trace.
The union of model predictions is a reachability constraint, never a pass rate.
"""
from pathlib import Path
import collections
import hashlib
import json
import sys
import gzip


def sha(path):
    h=hashlib.sha256()
    with path.open('rb')as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def stream(path):
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt')as f:
        for line in f:yield json.loads(line)


def kill_audit(path,rows,programs,statistics):
    seen=bytearray(programs);counts=collections.Counter();n=0
    for k in stream(path):
        i=k['candidate_id'];assert 0<=i<programs and not seen[i],i;seen[i]=1;n+=1
        r=rows[k['order']]
        assert k['trace_index']==r['trace_index']
        assert k['expected']==r['raw']['expected']!=k['got']
        if 'witness'in k:assert k['witness']==r['raw']
        counts[k['order']]+=1
    assert n==programs and all(seen),(n,programs)
    remaining=programs;stats_count=0;zeros=0
    for stat in stream(statistics):
        order=stat['order'];r=rows[order]
        assert stat['trace_index']==r['trace_index']and stat['id']==r['raw']['id']
        assert stat['killed']==counts[order]
        remaining-=stat['killed'];assert remaining==stat['remaining']
        zeros+=stat['killed']==0;stats_count+=1
    assert stats_count==len(rows)and remaining==0
    return dict(programs_killed_once=n,zero_kill_witnesses=zeros,
        kills_by_order={str(k):v for k,v in counts.items()if v})


def audit_r7(directory,native_name,software_name,kind):
    root=Path(directory);native=root/native_name;software=root/software_name
    ns=json.loads((native/'summary.json').read_text());ss=json.loads((software/'summary.json').read_text())
    assert ns['programs']==ss['programs']and ns['grammar_hash']==ss['grammar_hash']
    if kind=='boundary':
        rows=[r for r in stream(root/'frozen/scoped.jsonl')if r['leaves'][2]in(0,1<<63)]
        kill='all-fv0-first-kills.jsonl';stats='all-fv0-witness-kills.jsonl'
    else:
        rows=list(stream(root/'poisson-frozen/scoped.jsonl'))
        if kind=='count-one':rows=[r for r in rows if r['count']==1]
        kill='first-kills.jsonl';stats='witness-kills.jsonl'
    def hashes_match(name):return sha(native/name)==sha(software/name)
    elimination=kill_audit(native/kill,rows,ns['programs'],native/stats)
    same_kills=hashes_match(kill)
    if not same_kills:kill_audit(software/kill,rows,ss['programs'],software/stats)
    n=collected(native,'reachability-*.jsonl');s=collected(software,'reachability-*.jsonl')
    assert set(n)==set(s)==set(range(len(rows)))
    nr=collected(native,'unreachable-roots-*.jsonl');sr=collected(software,'unreachable-roots-*.jsonl')
    assert len(nr)==ns['whole_family_unreachable'];assert len(sr)==ss['whole_family_unreachable']
    output_diffs=[];reach_diffs=[];root_diffs=[];by_source=collections.defaultdict(collections.Counter)
    destination=root/f'audit-{native_name}';destination.mkdir(exist_ok=True)
    with (destination/'unreachable-witnesses.jsonl').open('x')as gaps:
        for i,r in enumerate(rows):
            a,b=n[i],s[i]
            for result in [a,b]:
                assert result['trace_index']==r['trace_index']and result['id']==r['raw']['id']
                assert result['expected']==r['raw']['expected']
                assert result['reachable']==(result['expected']in result['outputs'])
            if a['outputs']!=b['outputs']:output_diffs.append(dict(order=i,witness=r['raw'],native=a['outputs'],software=b['outputs']))
            if a['reachable']!=b['reachable']:reach_diffs.append(i)
            if i in nr and i in sr and set(nr[i]['raw_roots'])!=set(sr[i]['raw_roots']):root_diffs.append(i)
            by_source[r['raw']['source']]['rows']+=1
            by_source[r['raw']['source']]['native_unreachable']+=not a['reachable']
            by_source[r['raw']['source']]['software_unreachable']+=not b['reachable']
            if not a['reachable']:
                assert i in nr
                gaps.write(json.dumps(dict(order=i,trace_index=r['trace_index'],witness=r['raw'],outputs=a['outputs']),separators=(',',':'))+'\n')
    with(destination/'backend-output-differences.jsonl').open('x')as f:
        for d in output_diffs:f.write(json.dumps(d,separators=(',',':'))+'\n')
    summary=dict(rows=len(rows),programs=ns['programs'],native_unreachable=len(nr),software_unreachable=len(sr),
        first_kill_files_equal=same_kills,witness_kill_files_equal=hashes_match(stats),elimination=elimination,
        published_output_set_differences=len(output_diffs),reachability_class_differences=reach_diffs,
        unreachable_identity_symmetric_difference=sorted(set(nr)^set(sr)),
        unreachable_root_set_differences=root_diffs,unreachable_by_source=by_source,
        claim='Finite family sets/first kills only. Individually reachable targets are constraints, never calculated Excel passes.',
        fresh_validation_rows=0)
    (destination/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))


def read(path):
    return [json.loads(line) for line in path.open()]


def collected(directory, pattern):
    records = {}
    for path in sorted(directory.glob(pattern)):
        for record in read(path):
            assert record["order"] not in records
            records[record["order"]] = record
    return records


def audit(directory):
    root = Path(directory)
    rows = read(root / "frozen/scoped.jsonl")
    native = root / "native"
    software = root / "software"
    ns = json.loads((native / "summary.json").read_text())
    ss = json.loads((software / "summary.json").read_text())
    assert {k: v for k, v in ns.items() if k not in ["backend", "reachability_seconds"]} == {
        k: v for k, v in ss.items() if k not in ["backend", "reachability_seconds"]}
    n = collected(native, "reachability-*.jsonl")
    s = collected(software, "reachability-*.jsonl")
    assert n == s
    assert set(n) == set(range(len(rows)))
    nr = collected(native, "unreachable-roots-*.jsonl")
    sr = collected(software, "unreachable-roots-*.jsonl")
    assert nr == sr
    assert len(nr) == ns["whole_family_unreachable"]
    for file in ["first-kills.jsonl", "programs.jsonl", "survivors.json", "witness-kills.jsonl"]:
        assert (native / file).read_bytes() == (software / file).read_bytes(), file
    kills = read(native / "first-kills.jsonl")
    assert len(kills) == len({k["candidate_id"] for k in kills}) == 39366
    assert {k["candidate_id"] for k in kills} == set(range(39366))
    counts = collections.Counter()
    by_source = collections.defaultdict(collections.Counter)
    for k in kills:
        r = rows[k["order"]]
        assert k["witness"] == r["raw"]
        assert k["trace_index"] == r["trace_index"]
        assert k["expected"] == r["raw"]["expected"] != k["got"]
        counts[(k["order"], r["raw"]["id"])] += 1
    statistics = read(native / "witness-kills.jsonl")
    assert len(statistics) == len(rows)
    remaining = 39366
    with (root / "unreachable-witnesses.jsonl").open("w") as gaps, (root / "reachability-ordered.jsonl").open("w") as reachable:
        for order, r in enumerate(rows):
            record = n[order]
            assert record["trace_index"] == r["trace_index"]
            assert record["id"] == r["raw"]["id"] and record["source"] == r["raw"]["source"]
            assert record["expected"] == r["raw"]["expected"]
            assert record["reachable"] == (record["expected"] in record["outputs"])
            stats = statistics[order]
            assert stats["killed"] == counts[(order, r["raw"]["id"])]
            remaining -= stats["killed"]
            assert stats["remaining"] == remaining
            assert stats["id"] == r["raw"]["id"] and stats["trace_index"] == r["trace_index"]
            assert stats["reading"] == r["raw"]["reading"]
            by_source[r["raw"]["source"]]["scoped"] += 1
            by_source[r["raw"]["source"]]["unreachable"] += not record["reachable"]
            reachable.write(json.dumps(record, separators=(",", ":")) + "\n")
            if not record["reachable"]:
                assert order in nr
                gaps.write(json.dumps({"order": order, "witness": r["raw"], "frozen_leaves": r["leaves"],
                    "prior_gap": r["prior_gap"], "family_outputs": record["outputs"]}, separators=(",", ":")) + "\n")
    assert remaining == 0
    summary = {"native_software_prediction_differences": 0, "native_software_raw_root_differences": 0,
               "scoped_rows_reconciled": len(rows), "programs_killed_once": len(kills),
               "zero_kill_witnesses": len(rows) - len(counts),
               "nonzero_kill_witnesses": [{"order": o, "id": ident, "killed": count} for (o, ident), count in counts.items() if count],
               "unreachable_by_source": dict(by_source),
               "unreachable_witnesses_sha256": hashlib.sha256((root / "unreachable-witnesses.jsonl").read_bytes()).hexdigest(),
               "no_predictions_scored_as_excel_answers": True}
    # Accesses to Counter during reconciliation add zero entries; recount directly.
    summary["zero_kill_witnesses"] = sum(s["killed"] == 0 for s in statistics)
    (root / "audit-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


def audit_questions(directory, packet_name):
    root=Path(directory);packet=root/packet_name
    requests=json.loads((packet/'requests.json').read_text())
    assert set(requests)=={'function','witnesses'} and requests['function']=='POISSON.DIST'
    questions=requests['witnesses'];assert 0<len(questions)<=64
    ids=[q['id']for q in questions];assert len(set(ids))==len(ids)
    keys=[]
    for q in questions:
        assert set(q)=={'id','args'} and len(q['args'])==3
        assert q['args'][0]=='0x3ff0000000000000'
        assert q['args'][2] is False or q['args'][2]=='0x0000000000000000'
        keys.append(json.dumps(q['args'],separators=(',',':')))
    assert len(set(keys))==len(keys)
    old={json.dumps(r['args'],separators=(',',':'))for r in stream(root/'poisson-application-before.jsonl')}
    candidates=set(json.loads((packet/'hypotheses.json').read_text())['candidate_ids'])
    predictions=[]
    for backend in ['native','software']:
        records=list(stream(packet/f'predictions-{backend}.jsonl'))
        assert [r['id']for r in records]==ids
        encoded=[]
        for record in records:
            assert record['source_build_cohort']=='20527'
            outcomes={}
            for part in record['partitions']:
                assert part['candidate_ids']
                for c in part['candidate_ids']:
                    assert c not in outcomes;outcomes[c]=part['outcome']
            assert set(outcomes)==candidates
            encoded.append(outcomes)
        predictions.append(encoded)
    differences=sum(a[c]!=b[c]for a,b in zip(*predictions)for c in candidates)
    controls=json.loads((packet/'selection.json').read_text())['control_ids']
    assert len(controls)==10 and set(controls)<=set(ids)
    assert any(q['args'][1]=='0x4086800000000000'for q in questions) #720 exactly
    report=dict(requests=len(questions),means=len({q['args'][1]for q in questions}),answers_in_requests=0,
        new_typed_inputs_vs_decoded_trace=sum(k not in old for k in keys),rechecks_vs_decoded_trace=sum(k in old for k in keys),
        newness_scope='Against49674 decoded original inputs; four unread source records remain unresolved failures.',
        candidate_programs=len(candidates),prediction_partition_records=2*len(questions),
        candidate_prediction_backend_differences=differences,source_build_cohort='20527',
        global_family_survived=False,version_drift_confirmed=False,calculation_kernel_changed=False,
        request_sha256=sha(packet/'requests.json'))
    (packet/'validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__ == "__main__":
    if sys.argv[1]=="questions":audit_questions(sys.argv[2],sys.argv[3])
    elif len(sys.argv)>2:audit_r7(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
    else:audit(sys.argv[1])
