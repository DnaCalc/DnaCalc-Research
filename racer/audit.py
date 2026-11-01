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


if __name__ == "__main__":
    if sys.argv[1]=="structural-questions":structural_questions(sys.argv[2])
    elif sys.argv[1]=="pmt-group-questions":pmt_group_questions(*sys.argv[2:])
    elif sys.argv[1]=="pmt-observation-questions":pmt_observation_questions(*sys.argv[2:])
    elif sys.argv[1]=="questions":audit_questions(sys.argv[2],sys.argv[3])
    elif len(sys.argv)>2:audit_r7(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
    else:audit(sys.argv[1])
