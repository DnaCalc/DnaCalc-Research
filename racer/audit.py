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


def pmt_observation_questions(trace, certificates, destination):
    """Freeze stratified rechecks; no historical value enters the request file."""
    import datetime
    import struct
    from collections import defaultdict, Counter
    rows = list(stream(Path(trace)))
    gaps = {r['id']: r for r in stream(Path(certificates))}
    assert len(gaps) == 734
    def number(s): return struct.unpack('>d', bytes.fromhex(s[2:]))[0]
    def stratum(r):
        a = list(map(number, r['args']))
        return (r['source'], a[1] == 1, a[3] != 0,
                'tiny' if a[0] < 1e-6 else 'small' if a[0] < .01 else 'ordinary',
                a[4] != 0, a[1].is_integer(), a[2] < 0)
    def choose(pool, count):
        cells = defaultdict(list)
        for r in pool:
            cells[stratum(r)].append(r)
        for values in cells.values():
            values.sort(key=lambda r: hashlib.sha256(json.dumps(r, sort_keys=True).encode()).hexdigest())
        chosen = []
        for index in range(count):
            available = [k for k, v in cells.items() if v]
            key = min(available, key=lambda k: (sum(stratum(r) == k for r in chosen), str(k)))
            chosen.append(cells[key].pop(0))
        return chosen
    selected = []
    for build, quota in [('20131', 16), ('20228', 8)]:
        pool = [r for r in rows if r['id'] in gaps and r['reading']['excel_build'] == build]
        for r in pool:
            g = gaps[r['id']]
            assert r['args'][0][2:].lower() == g['rate_bits'].lower()
            assert r['expected'][2:].lower() == g['target_bits'].lower()
        selected += [(r, 'certified_gap') for r in choose(pool, quota)]
    for build in ['20131', '20228']:
        pool = [r for r in rows if r['id'] not in gaps and r['actual'] == r['expected']
                and r['reading']['excel_build'] == build and len(r['args']) == 5
                and all(isinstance(a, str) and a.startswith('0x') for a in r['args'])
                and number(r['args'][0]) > 0 and number(r['args'][1]) > 0
                and number(r['args'][4]) in (0, 1)]
        selected += [(r, 'passing_control') for r in choose(pool, 4)]
    out = Path(destination)
    assert not (out/'requests-PMT.json').exists()
    out.mkdir(parents=True, exist_ok=True)
    questions, selection, predictions = [], [], []
    def add(index, kind):
        old, role = selected[index]
        args = old['args'].copy()
        transport = 'cell_value2_bulk'
        if kind == 'literal': transport = 'formula_literals'
        elif kind == 'logical_timing': args[4] = number(args[4]) != 0
        elif kind == 'text_rate': args[0] = format(number(args[0]), '.17g')
        q = dict(id=f'r10-pmt-{len(questions):02d}', args=args, input_transport=transport)
        questions.append(q)
        selection.append(dict(id=q['id'], original_id=old['id'], original_source=old['source'],
            original_reading=old['reading'], sample_role=role, variant=kind,
            paired_numeric_id=f'r10-pmt-{index:02d}', stratum=list(stratum(old)),
            argument_types=['logical' if isinstance(a, bool) else 'binary64' if a.startswith('0x') else 'text' for a in args],
            requested=args, submitted='typed Value2 cells' if transport == 'cell_value2_bulk' else
            '17-decimal-place scientific literals in PMT and echo cells; adapter emits exact strings',
            read_back=None, effective_callee_inputs='not directly observable'))
        predictions.append(dict(id=q['id'], hypothesis='build-invariant PMT and transparent transport/coercion',
            predicted=old['expected'], historical_observation=old,
            status='conditional prediction, not a new Excel answer; coercion variants have no historical observation of those types'))
    for index in range(32): add(index, 'numeric')
    # Equal build quotas; include both historic sources and passing controls.
    for index in list(range(0, 8)) + list(range(16, 20)) + [24, 25, 28, 29]: add(index, 'literal')
    for index in [0, 2, 4, 6, 16, 18, 20, 22]: add(index, 'logical_timing')
    for index in [1, 3, 5, 7, 17, 19, 21, 23]: add(index, 'text_rate')
    assert len(questions) == len({q['id'] for q in questions}) == 64
    request = dict(function='PMT', selection='R10 observation audit; intentional rechecks, do not use --new-only', witnesses=questions)
    (out/'requests-PMT.json').write_text(json.dumps(request, indent=2)+'\n')
    (out/'PMT-selection.json').write_text(json.dumps(dict(method='deterministic balanced round-robin strata, hash order within each stratum; no numeric fitting',
        unique_historical_rows=32, gaps=24, passing_controls=8, numeric_build_counts={'20131':20,'20228':12},
        legacy_transport='not recorded; not presumed Value2', target_build='20527; actual build must be recorded and compared', questions=selection), indent=2)+'\n')
    (out/'PMT-predictions-historical.jsonl').write_text(''.join(json.dumps(p, sort_keys=True)+'\n' for p in predictions))
    (out/'README.txt').write_text("""Round 10: PMT observation audit (64 calls, zero new answers).
From the updated DnaCalc master on Windows:
 cargo xtask witness read --func PMT --inputs <packet>/requests-PMT.json
Do not pass --new-only. These are intentional rechecks and transport pairs.
Run once, retain the frozen predictions, readings file and entire worker directory.
Compare all 32 numeric rows with PMT-predictions-historical.jsonl by ID.
Then compare paired literals/logical timing/text-rate results with their numeric row.
Only unchanged read-back inputs support same-input build comparisons. Literal echo
cells expose the decimal parser, not effective callee internals. Changed inputs and
unread results remain failed audit obligations, never semantic contradictions.
The adapter retains input_audits even for changed/refused rows. Preserve its raw
request.json/result.json. Record actual build, CPU/bitness and settings; do not
silently attribute build/CPU/surface differences to any one cause. No old observation
is overwritten. Historical transport and Compatibility Version are unknown.
32 numeric calls: 24 certified gaps (16 build20131,8 build20228),8 exact passing
controls (4 per build). Plus16 literal pairs,8 logical-timing,8 text-rate pairs.
Only requests-PMT.json is submitted. Predictions and old answers are separate.
This audit cannot claim build20527 reproduction until Windows returns readings.
""")
    files = [out/name for name in ['requests-PMT.json', 'PMT-selection.json', 'PMT-predictions-historical.jsonl', 'README.txt']]
    manifest = dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), function='PMT',
        calls=64, new_excel_answers=0, input_trace_sha256=sha(Path(trace)), certificate_sha256=sha(Path(certificates)),
        files_sha256={p.name:sha(p) for p in files}, assumptions='conditional legacy transport/build invariance; not confirmed')
    (out/'PMT-manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps(manifest))


def pmt_group_questions(scoped, certificates, destination):
    """Recheck sufficient contradiction rows and neighboring amount/timing inputs.

    These are observation audits. No candidate survives the source constraints,
    and no old answer is supplied to the Excel request.
    """
    import datetime
    import math
    import struct
    root=Path(certificates)
    selected={}
    proof_files=[]
    for family in ('timing-complement','annuity-factor'):
        for path in sorted((root/family).glob('*-first-kills.jsonl')):
            proof_files.append(path)
            for certificate in stream(path):
                for witness in certificate['proof']['witnesses']:
                    index=witness['trace_index'];raw=witness['witness']
                    assert index not in selected or selected[index]==raw
                    selected[index]=raw
    assert selected and len(selected)<=16
    originals={}
    for row in stream(Path(scoped)):
        if row['trace_index']in selected:
            assert row['raw']==selected[row['trace_index']]
            originals[row['trace_index']]=row
    assert originals.keys()==selected.keys()
    out=Path(destination);out.mkdir(parents=True,exist_ok=False)
    def number(bits):return struct.unpack('>d',bytes.fromhex(bits[2:]))[0]
    def encoded(x):return '0x'+struct.pack('>d',x).hex()
    questions=[];bindings=[];historical=[];seen=set()
    # The priority is fixed before predictions: every sufficient original row,
    # each opposite timing, then adjacent values, then binary scalings.
    for variant in ('original','opposite_timing','previous_amount','next_amount','half_amount','double_amount'):
        for index,raw in sorted(selected.items()):
            args=raw['args'].copy();pv=number(args[2])
            if variant=='opposite_timing':args[4]=encoded(float(number(args[4])==0))
            elif variant=='previous_amount':args[2]=encoded(math.nextafter(pv,-math.inf))
            elif variant=='next_amount':args[2]=encoded(math.nextafter(pv,math.inf))
            elif variant=='half_amount':args[2]=encoded(pv/2)
            elif variant=='double_amount':args[2]=encoded(pv*2)
            key=tuple(args)
            if key in seen or len(questions)==64:continue
            seen.add(key);ident=f'r12-pmt-group-{len(questions):02d}'
            questions.append(dict(id=ident,args=args,input_transport='cell_value2_bulk'))
            bindings.append(dict(id=ident,variant=variant,source_trace_index=index,
                original_id=raw['id'],original_source=raw['source'],original_reading=raw['reading'],
                source_arguments=raw['args'],requested_arguments=args,
                source_input_readback=None,effective_callee_inputs='not directly observable'))
            if variant=='original':historical.append(dict(id=ident,historical_observation=raw,
                status='Original recorded observation only; verify current read-back inputs before semantic comparison.'))
    assert len(questions)==64 and len(historical)==len(selected)
    request=dict(function='PMT',selection='R12 shared-helper observation audit; intentional rechecks, do not use --new-only',witnesses=questions)
    with(out/'requests-PMT.json').open('x')as f:json.dump(request,f,indent=2);f.write('\n')
    with(out/'selection.json').open('x')as f:json.dump(dict(method='Sufficient first-contradiction rows, fixed amount/timing variant priority; no ranking by prediction or fit.',
        original_rows=len(selected),calls=len(questions),bindings=bindings),f,indent=2);f.write('\n')
    with(out/'historical-observations.jsonl').open('x')as f:
        for row in historical:f.write(json.dumps(row,sort_keys=True)+'\n')
    with(out/'README.txt').open('x')as f:f.write('''Round12 PMT shared-helper input audit:64 calls, zero new answers.
On current DnaCalc master on Windows:
 cargo xtask witness read --func PMT --inputs <packet>/requests-PMT.json
Do not pass --new-only. Preserve the readings and entire raw worker directory,
including input_audits, request.json and result.json. Submit only requests-PMT.json.
The13 original rows establish exact contradictions for the earlier shared-helper
projections, conditional on their recorded inputs. Older callee inputs were not
read back. Compare original rows with historical-observations.jsonl only after
verifying the new Value2 read-back. Opposite timings, adjacent amounts and binary
scalings audit amount dependence and transport; their outputs have no Excel label
until this capture. Predictions remain separate and never count as observations.
Record actual build, CPU/bitness, calculation settings and Compatibility Version.
Neither a changed input nor a failed read may be normalized into a semantic match.
No version drift, transport defect, helper or complete PMT program is established
by preparing this packet. The older R11 structural and cumulative packets remain
unanswered and unchanged.
''')
    files=sorted(p for p in out.iterdir()if p.is_file())
    manifest=dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),calls=64,
        new_excel_answers=0,function='PMT',source_trace_sha256=sha(Path(scoped)),
        proof_sha256={str(p.relative_to(root)):sha(p)for p in proof_files},
        files_sha256={p.name:sha(p)for p in files},
        limitations='Observation/transport audit only; no surviving candidate is claimed. Historical effective inputs remain unknown.')
    with(out/'manifest.json').open('x')as f:json.dump(manifest,f,indent=2);f.write('\n')
    print(json.dumps(manifest))


def structural_questions(out):
    """R11 domain-boundary and conversion-route observations, not fitted bodies."""
    import math
    import struct
    import datetime
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    bits=lambda x:'0x'+struct.pack('>d',x).hex()
    batches={f:[]for f in ('BESSELI','BESSELJ','CONCATENATE','CONCAT')}
    hypotheses=[];selection=[]
    def add(function,args,role,hypothesis=None):
        qid=f'r11-struct-{function.lower()}-{len(batches[function]):02}'
        batches[function].append(dict(id=qid,args=args,input_transport='cell_value2_bulk'))
        selection.append(dict(id=qid,role=role,exposure='all previously inspected banks and corpus are discovery; this is a route/domain audit, not validation of a complete numeric candidate'))
        if hypothesis is not None:hypotheses.append(dict(id=qid,kind_only_boundary_hypothesis=hypothesis,not_a_complete_value_prediction=True))
    limits={'BESSELI':709.7,'BESSELJ':float(2**27)}
    for function,limit in limits.items():
        extra=709.7000001 if function=='BESSELI' else 134217727.5
        xs=[math.nextafter(limit,0.0),limit,math.nextafter(limit,math.inf),extra]
        for order in (0,1,2):
            for x in xs:
                refused=abs(x)>limit if function=='BESSELI'else abs(x)>=limit
                add(function,[bits(x),bits(float(order))],
                    'adjacent binary64 boundary/control; orders distinguish a common argument bound from a body-specific failure',
                    'error:Num'if refused else'finite number; exact bits unspecified')
        for x in xs:
            refused=abs(x)>limit if function=='BESSELI'else abs(x)>=limit
            add(function,[bits(-x),bits(0.0)],'negative-argument boundary/control',
                'error:Num'if refused else'finite number; exact bits unspecified')
    values=[5758361573115.375,521310271410497.5,-6278515587045.875,
        543212677140.0625,5.248519552626765e-6,7200547821627.375,
        999999999999999.5,1234567890123455.0,123456789012345.5,4037429.935952025]
    for function in ('CONCATENATE','CONCAT'):
        for i,x in enumerate(values):
            add(function,[bits(x)],('known ampersand discrepancy'if i<5 else'known conversion control')+
                '; compare the actual function route with the existing ampersand capture, without changing rounding rules')
    for function,rows in batches.items():
        (out/f'requests-{function}.json').write_text(json.dumps(dict(function=function,witnesses=rows),indent=2)+'\n')
    (out/'selection.json').write_text(json.dumps(selection,indent=2)+'\n')
    (out/'boundary-hypotheses.json').write_text(json.dumps(hypotheses,indent=2)+'\n')
    (out/'RUN-STRUCTURAL-R11.txt').write_text('''Round 11 structural/number-text audit: 52 calls, no new answers.
From the updated DnaCalc master on Windows, substitute this packet's path:
 cargo xtask witness read --func BESSELI --inputs <packet>/requests-BESSELI.json
 cargo xtask witness read --func BESSELJ --inputs <packet>/requests-BESSELJ.json
 cargo xtask witness read --func CONCATENATE --inputs <packet>/requests-CONCATENATE.json
 cargo xtask witness read --func CONCAT --inputs <packet>/requests-CONCAT.json
Do not use --new-only. Preserve every raw result and input audit. Score immutable
registered-route predictions first. Boundary hypotheses predict kind only and
must never be reported as exact function matches. They ask whether an argument
admission rule has an independently distinguishable adjacent-input transition;
no numerical body or coefficient is selected from their results.
''')
    print(json.dumps({f:len(rows)for f,rows in batches.items()}))


def bank_delta(before, after, output, mode='file'):
    """Compare sorted canonical witness identities while retaining new coverage."""
    before,after,output=map(Path,(before,after,output))
    if mode not in ('file','stream'):
        raise ValueError('bank delta mode must be file or stream')
    frozen={str(p):sha(p)for p in ((before,after,Path(__file__))if mode=='file'else(Path(__file__),))}
    def key(row):
        selection=row.get('result_index')
        suffix=''if selection is None else f'@{selection[0]},{selection[1]}'
        return row['function'],json.dumps(row['args'],ensure_ascii=False,sort_keys=True,separators=(',',':'))+suffix
    def ordered(path):
        previous=None
        for row in stream(path):
            current=key(row)
            if previous is not None and previous>=current:
                raise ValueError(f'canonical witness order/uniqueness differs: {previous[:300]!r} >= {current[:300]!r}')
            previous=current;yield current,row
    a,b=iter(ordered(before)),iter(ordered(after));old,new=next(a,None),next(b,None)
    counts=collections.Counter();examples=[]
    while old is not None or new is not None:
        if old is None or new is not None and new[0]<old[0]:
            counts['added_inputs']+=1;counts['added_matched']+=new[1]['expected']==new[1]['actual']
            counts['added_failed']+=new[1]['expected']!=new[1]['actual'];new=next(b,None);continue
        if new is None or old[0]<new[0]:
            counts['removed_inputs']+=1;old=next(a,None);continue
        x,y=old[1],new[1];counts['shared_inputs']+=1
        counts['unchanged_outputs']+=x['actual']==y['actual']
        counts['changed_outputs']+=x['actual']!=y['actual']
        counts['changed_expected']+=x['expected']!=y['expected']
        counts['representative_metadata_changes']+={k:v for k,v in x.items()if k not in ('actual','expected')}!={k:v for k,v in y.items()if k not in ('actual','expected')}
        mx,my=x['actual']==x['expected'],y['actual']==y['expected']
        counts['shared_before_matches']+=mx;counts['shared_after_matches']+=my
        counts['shared_gains']+=my and not mx;counts['shared_losses']+=mx and not my
        if (x['actual']!=y['actual']or x['expected']!=y['expected'])and len(examples)<20:
            examples.append(dict(before=x,after=y))
        old,new=next(a,None),next(b,None)
    if frozen!={p:sha(Path(p))for p in frozen}:raise ValueError('bank/source changed during comparison')
    result=dict(complete=True,counts=dict(counts),examples=examples,input_sha256=frozen,
        limits='Counts distinguish input coverage from gains/losses on shared inputs. Metadata '
               'changes retain the selected capture source. Unread rows have no trace and '
               'still fail in each original judge summary; no calculated-only pass rate.')
    if mode=='stream':
        result['trace_mode']='streamed once; input bindings and judge completion supplied by bank-intake receipt'
    with output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result));return result


def bank_intake(executable, checkout, before_receipt, output):
    """Compare two evidence views through the existing judge, without trace dumps."""
    import os
    import subprocess
    import tempfile
    executable,checkout,before_receipt,output=(Path(p).resolve()for p in(executable,checkout,before_receipt,output))
    output.mkdir(parents=True,exist_ok=False)
    baseline=json.loads(before_receipt.read_text())
    if baseline.get('complete')is not True:
        raise ValueError('incomplete pre-intake evidence binding')
    old_scope=baseline['bank_scope'];evidence=Path(old_scope['root'])
    if old_scope.get('status') or old_scope.get('dirty_or_untracked_sha256'):
        raise ValueError('intake comparison requires a clean pre-intake evidence revision')
    frozen={str(p):sha(p)for p in(executable,before_receipt,Path(__file__))}
    def git(*args):
        return subprocess.check_output(['git','-C',str(evidence),*args],text=True,timeout=60).strip()
    old_paths=set(git('ls-tree','-r','--name-only',old_scope['revision']).splitlines())
    changed=set(git('diff','--name-only',old_scope['revision']).splitlines())
    if (changed&old_paths)-{'.gitattributes','README.md'}:
        raise ValueError('pre-intake bank files changed')
    def state():
        readings={str(p.relative_to(evidence)):sha(p)for p in sorted((evidence/'readings').rglob('*.json'))}
        return dict(revision=git('rev-parse','HEAD'),status=git('status','--porcelain=v1','--untracked-files=all'),
                    readings_sha256=readings)
    new_scope=state()
    if not all(new_scope['readings_sha256'].get(name)==h for name,h in old_scope['readings_sha256'].items()):
        raise ValueError('original reading missing or changed during intake')
    added=sorted(new_scope['readings_sha256'].keys()-old_scope['readings_sha256'].keys())
    if not added:
        raise ValueError('empty intake comparison')
    commands=[]
    with tempfile.TemporaryDirectory(prefix='dnacalc-intake-')as scratch:
        root=Path(scratch);view=root/'before'/'OxFunc-Evidence';view.mkdir(parents=True)
        for source in evidence.iterdir():
            if source.name not in ('.git','readings'):
                (view/source.name).symlink_to(source,target_is_directory=source.is_dir())
        for name in old_scope['readings_sha256']:
            path=view/name;path.parent.mkdir(parents=True,exist_ok=True);path.symlink_to(evidence/name)
        before,after=root/'before.fifo',root/'after.fifo'
        for path in(before,after):os.mkfifo(path,0o600)
        children=[]
        with (output/'before-judge.log').open('x')as old_log,(output/'judge.log').open('x')as new_log:
            try:
                for trace,legacy,log in((before,view.parent,old_log),(after,evidence.parent,new_log)):
                    command=[str(executable),'witness','judge','--bank','full','--trace',str(trace)]
                    env=dict(os.environ,DNACALC_LEGACY_ROOT=str(legacy))
                    commands.append(dict(argv=command,legacy_root=str(legacy)))
                    # ubs:ignore — selected local binary, argv only; these two owned children are bounded and always reaped.
                    children.append(subprocess.Popen(command,cwd=checkout,env=env,stdout=log,stderr=subprocess.STDOUT))
                delta=bank_delta(before,after,output/'delta.json','stream')
                if any(child.wait(timeout=3600)for child in children):
                    raise ValueError('bank intake judge failed')
            finally:
                for child in children:
                    if child.poll()is None:
                        child.terminate()
                        try:child.wait(timeout=5)
                        except subprocess.TimeoutExpired:child.kill();child.wait()
    if frozen!={p:sha(Path(p))for p in frozen}or state()!=new_scope:
        raise ValueError('intake source, executable or evidence changed during comparison')
    import re
    def totals(name):
        match=re.search(r'; (\d+) of (\d+) rows matched, (\d+) failed .*; (\d+) unread\)',(output/name).read_text())
        if match is None:raise ValueError('missing intake judge totals')
        matched,total,failed,unread=map(int,match.groups())
        if failed!=total-matched:raise ValueError('incorrect intake judge failure count')
        return dict(matched=matched,total=total,failed=failed,unread=unread)
    old,new=totals('before-judge.log'),totals('judge.log');c=delta['counts']
    if old!={key:baseline[key]for key in('matched','total','failed','unread')}:
        raise ValueError('before-view judge differs from the bound pre-intake totals')
    if c.get('removed_inputs',0)or c['shared_inputs']+old['unread']!=old['total']or c['shared_inputs']+c.get('added_inputs',0)+new['unread']!=new['total']:
        raise ValueError('intake coverage does not reconcile')
    if c['shared_before_matches']!=old['matched']or c['shared_after_matches']+c.get('added_matched',0)!=new['matched']:
        raise ValueError('intake matches do not reconcile')
    result=dict(complete=True,before=old,after=new,counts=c,commands=commands,
        before_bank_scope=old_scope,after_bank_scope=new_scope,added_reading_files=added,
        input_sha256=frozen,logs_sha256={p.name:sha(p)for p in output.iterdir()if p.is_file()},
        limits='Same application executable on two evidence views. Existing cache/data files are shared unchanged; the before view exposes exactly the originally bound reading files. '
               'All identities/answers and typed actuals use the existing canonical judge and bank-delta reconciliation. Added coverage and selected-capture metadata are separate from product gains. '
               'Unread observations remain failures. Two FIFO traces are consumed once and never retained; no target or bulk bank copy.')
    with(output/'receipt.json').open('x')as file:json.dump(result,file,indent=2);file.write('\n')
    print(json.dumps(dict(before=old,after=new,counts=c)));return result


def bank_replay(before, executable, checkout, output, mode='file', function=None):
    """Replay through the existing judge/comparator without retaining a new trace."""
    import os
    import subprocess
    import tempfile
    before, executable, checkout, output = (Path(p).resolve() for p in
        (before, executable, checkout, output))
    if mode not in ('file', 'live'):
        raise ValueError('bank replay mode must be file or live')
    output.mkdir(parents=True, exist_ok=False)
    sources = {str(p): sha(p) for p in (before, executable, Path(__file__))}
    legacy_root = Path(os.environ['DNACALC_LEGACY_ROOT']) if 'DNACALC_LEGACY_ROOT' in os.environ else checkout.parent
    evidence = (legacy_root if legacy_root.is_absolute() else checkout/legacy_root).resolve()/'OxFunc-Evidence'
    def evidence_state():
        def git(*args):
            return subprocess.check_output(['git', '-C', str(evidence), *args], timeout=60)
        changed = set(git('diff', '--name-only', '-z', 'HEAD').decode().split('\0'))
        changed.update(git('ls-files', '--others', '--exclude-standard', '-z').decode().split('\0'))
        dirty = {name: sha(evidence/name) if (evidence/name).is_file() else None
                 for name in sorted(changed) if name}
        return dict(root=str(evidence), revision=git('rev-parse', 'HEAD').decode().strip(),
            status=git('status', '--porcelain=v1', '--untracked-files=all').decode(),
            dirty_or_untracked_sha256=dirty,
            readings_sha256={str(p.relative_to(evidence)): sha(p) for p in sorted((evidence/'readings').rglob('*.json'))})
    bank_scope = evidence_state()
    commands = []
    with tempfile.TemporaryDirectory(prefix='dnacalc-bank-') as scratch:
        fifo = Path(scratch)/'trace.fifo'
        os.mkfifo(fifo, 0o600)
        old_trace = before
        old_judge = None
        if mode == 'live':
            old_trace = Path(scratch)/'before.fifo'
            os.mkfifo(old_trace, 0o600)
            old_judge = [str(before), 'witness', 'judge', '--bank', 'full', '--trace', str(old_trace)]
            if function:
                old_judge.extend(('--func', function))
        compare = [str(executable), 'witness', 'compare', '--before', str(old_trace), '--after', str(fifo)]
        judge = [str(executable), 'witness', 'judge', '--bank', 'full', '--trace', str(fifo)]
        if function:
            judge.extend(('--func', function))
        commands.extend((compare, judge))
        if old_judge is not None:
            commands.append(old_judge)
        with (output/'compare.log').open('x') as compare_log, (output/'judge.log').open('x') as judge_log:
            # ubs:ignore — operator-selected local executable and argv lists, never a shell.
            reader = subprocess.Popen(compare, cwd=checkout, stdout=compare_log, stderr=subprocess.STDOUT)
            writer = None
            old_writer = None
            old_log = None
            try:
                if old_judge is not None:
                    old_log = (output/'before-judge.log').open('x')
                    # ubs:ignore — operator-selected local executable and argv, never a shell.
                    old_writer = subprocess.Popen(old_judge, cwd=checkout, stdout=old_log, stderr=subprocess.STDOUT)
                # ubs:ignore — operator-selected local executable and argv lists, never a shell.
                writer = subprocess.Popen(judge, cwd=checkout, stdout=judge_log, stderr=subprocess.STDOUT)
                writer_code = writer.wait(timeout=3600)
                old_code = old_writer.wait(timeout=60) if old_writer is not None else 0
                reader_code = reader.wait(timeout=60)
                if writer_code or reader_code or old_code:
                    raise ValueError(f'judge/comparison failed ({old_code}, {writer_code}, {reader_code}); retain all logs')
            finally:
                for child in (writer, old_writer, reader):
                    if child is not None and child.poll() is None:
                        child.terminate()
                        try:
                            child.wait(timeout=5)
                        except subprocess.TimeoutExpired:
                            child.kill()
                            child.wait()
                if old_log is not None:
                    old_log.close()
    if sources != {p: sha(Path(p)) for p in sources}:
        raise ValueError('bank, executable or orchestration source changed during replay')
    if bank_scope != evidence_state():
        raise ValueError('evidence revision, dirty files or readings changed during replay')
    import re
    compared = (output/'compare.log').read_text()
    summary = re.search(r'xtask witness compare: (\d+) -> (\d+) matched, \+(\d+) / -(\d+); (\d+) actuals changed', compared)
    rows = re.search(r'xtask witness compare: (\d+) decoded rows;', compared)
    judged = re.search(r'; (\d+) of (\d+) rows matched, (\d+) failed .*; (\d+) unread\)', (output/'judge.log').read_text())
    if summary is None or rows is None or judged is None or int(rows[1]) == 0:
        raise ValueError('missing nonempty judge/comparison totals')
    old, new, gains, losses, changed = map(int, summary.groups())
    matched, total, failed, unread = map(int, judged.groups())
    if new != matched or total != int(rows[1])+unread or failed != total-matched:
        raise ValueError('judge and comparator totals do not reconcile')
    if mode == 'live':
        previous = re.search(r'; (\d+) of (\d+) rows matched, (\d+) failed .*; (\d+) unread\)', (output/'before-judge.log').read_text())
        if previous is None or tuple(map(int, previous.groups())) != (old, total, total-old, unread):
            raise ValueError('baseline judge and comparison totals do not reconcile')
    result = dict(complete=True, decoded_rows=int(rows[1]), before_matched=old, matched=matched,
        total=total, failed=failed, unread=unread, gains=gains, losses=losses, changed_outputs=changed,
        commands=commands, mode=mode, function=function, input_sha256=sources, bank_scope=bank_scope,
        logs_sha256={p.name: sha(p) for p in output.glob('*.log')},
        limits='Existing application judge and verdict comparator; all identity/answer fields compared. '
               'Unread observations remain failures. The FIFO trace is consumed once and never retained.')
    with (output/'receipt.json').open('x') as file:
        json.dump(result, file, indent=2)
        file.write('\n')
    print(json.dumps(result))
    return result


def corpus_delta(before, after, original_catalogue_after, output):
    """Compare fixed observations; '-' uses every matching catalogue session."""
    override_catalogue = str(original_catalogue_after) != '-'
    before, after, original_catalogue_after, output = map(Path,
        (before, after, original_catalogue_after, output))
    document = lambda path: json.loads(path.read_text())
    summary = document(before/'summary.json')
    counts = collections.Counter()
    changes, pairs, workbooks = [], [], set()
    hashes = {str(before/'summary.json'): sha(before/'summary.json'), str(Path(__file__)): sha(Path(__file__))}
    for old_report in sorted(before.glob('workbooks/**/report.json')):
        relative = old_report.relative_to(before)
        new_report = original_catalogue_after/'report.json' if override_catalogue and 'calc-func_catalogue_oracle.xlsx' in relative.parts else after/relative
        if not override_catalogue:
            captures = [p.parent.parent/(p.parent.name+'.capture.json') for p in (old_report, new_report)]
            if any(p.exists() for p in captures):
                if not all(p.is_file() for p in captures) or sha(captures[0]) != sha(captures[1]):
                    raise ValueError('capture observations changed: '+str(relative))
                for p in captures: hashes[str(p)] = sha(p)
        paths = [old_report, new_report, old_report.with_name('evaluations.json'), new_report.with_name('evaluations.json')]
        old_triage, new_triage = old_report.with_name('triage.json'), new_report.with_name('triage.json')
        for p in (old_triage, new_triage): hashes[str(p)] = sha(p)
        for p in paths: hashes[str(p)] = sha(p)
        old, new, old_evaluations, new_evaluations = map(document, paths)
        key = lambda r: (r['sheet'], r['row'], r['column'])
        policy = {key(r): r for r in document(old_triage)['cells']}
        new_policy = {key(r): r for r in document(new_triage)['cells']}
        a = {key(r): r for r in old['cells']}
        b = {key(r): r for r in new['cells']}
        if a.keys() != b.keys(): raise ValueError('original formula inventory changed: '+str(relative))
        x = {key(r): r for r in old_evaluations['evaluations']}
        y = {key(r): r for r in new_evaluations['evaluations']}
        if x.keys() != y.keys(): raise ValueError('original evaluation coverage changed: '+str(relative))
        workbooks.add(relative.parts[1])
        counts['sessions'] += 1
        counts['formula_cells'] += len(a)
        counts['evaluations'] += len(x)
        for k in a:
            if (a[k]['verdict'] == 'awaiting_capture') != (b[k]['verdict'] == 'awaiting_capture'):
                raise ValueError('capture availability changed; score coverage separately: '+str(k))
            excluded = bool(policy[k].get('external_links') or policy[k].get('circular_rule'))
            if excluded != bool(new_policy[k].get('external_links') or new_policy[k].get('circular_rule')):
                raise ValueError('original comparison context changed: '+str(k))
            if excluded:
                counts['external_or_circular'] += 1
                continue
            was, now = a[k]['verdict'] == 'matched', b[k]['verdict'] == 'matched'
            counts['before_matched'] += was
            counts['after_matched'] += now
            counts['gains'] += now and not was
            counts['losses'] += was and not now
            if a[k]['verdict'] != b[k]['verdict']:
                changes.append(dict(workbook=str(relative), cell=k, before=a[k], after=b[k],
                                    original_evaluation=x.get(k), current_evaluation=y.get(k)))
        for k in x:
            counts['changed_values'] += x[k].get('returned') != y[k].get('returned')
        pairs.append(dict(before=str(old_report), after=str(new_report), formula_cells=len(a)))
    counts['workbooks'] = len(workbooks)
    expected = summary['totals']['counts']
    if counts['workbooks'] != summary['totals']['workbooks'] or counts['formula_cells'] != expected['formula_cells'] or counts['before_matched'] != expected['pass']:
        raise ValueError('original corpus totals do not reconcile: '+str(dict(counts)))
    if any(sha(Path(p)) != h for p,h in hashes.items()): raise ValueError('corpus/source changed during comparison')
    result = dict(complete=True, counts=dict(counts), changes=changes, pairs=pairs, inputs_sha256=hashes,
        external=expected['external_links'], circular=expected['circular_rule'],
        before_fail=expected['fail'], after_fail=expected['fail']-counts['gains']+counts['losses'],
        limits='Fixed formula inventory, answers and workbook contexts, including every comparison session. '
               'New capture coverage is separate from product gains; no answer or scope-dependent pass is synthesized.')
    with output.open('x')as file:json.dump(result,file,indent=2);file.write('\n')
    print(json.dumps(dict(counts=dict(counts), before_fail=result['before_fail'], after_fail=result['after_fail'])))
    return result


if __name__ == "__main__":
    if sys.argv[1]=="bank-intake":bank_intake(*sys.argv[2:]);raise SystemExit(0)
    if sys.argv[1]=="corpus-delta":corpus_delta(*sys.argv[2:]);raise SystemExit(0)
    if sys.argv[1]=="bank-delta":bank_delta(*sys.argv[2:]);raise SystemExit(0)
    if sys.argv[1]=="bank-replay":bank_replay(*sys.argv[2:]);raise SystemExit(0)
    if sys.argv[1]=="structural-questions":structural_questions(sys.argv[2])
    elif sys.argv[1]=="pmt-group-questions":pmt_group_questions(*sys.argv[2:])
    elif sys.argv[1]=="pmt-observation-questions":pmt_observation_questions(*sys.argv[2:])
    elif sys.argv[1]=="questions":audit_questions(sys.argv[2],sys.argv[3])
    elif len(sys.argv)>2:audit_r7(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
    else:audit(sys.argv[1])
