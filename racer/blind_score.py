"""Score frozen R7 predictions before a separate exact-elimination step.

Only consumes original xtask judge traces, never a second oracle importer.
All requested identities must have a build-20527 reading from the fresh packet.
Predictions remain research data and never enter an answer bank.
"""
import collections as C
import csv
import datetime
import hashlib
import json
from pathlib import Path
import sys


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_packet(packet):
    manifest=json.loads((packet/'manifest.json').read_text())
    for name,digest in manifest['files_sha256'].items():
        assert sha(packet/name)==digest, ('frozen file changed',name)
    requests=json.loads((packet/'requests.json').read_text())['witnesses']
    hypotheses=json.loads((packet/'hypotheses.json').read_text())
    assert len(requests)==38 and len(hypotheses['candidate_ids'])==52908
    return requests,hypotheses,manifest


def identity(row):
    return json.dumps([row['args'],row.get('result_index')],sort_keys=True,separators=(',',':'))


def score(packet,trace,out):
    requests,hypotheses,manifest=load_packet(packet)
    keys={identity(r):r['id']for r in requests}
    assert len(keys)==len(requests)
    answers={}
    for i,line in enumerate(open(trace)):  # ubs:ignore — finite one-shot iterator; fully consumed and released before process exit
        row=json.loads(line)
        if row['function']!='POISSON.DIST' or identity(row) not in keys:
            continue
        # Require the newly captured packet, not a same-input legacy reading.
        if str(row['reading'].get('excel_build'))!='20527' or not row['source'].endswith('POISSON.DIST-calcvba7.json'):
            continue
        qid=keys[identity(row)]
        if qid in answers:
            assert answers[qid]['expected']==row['expected'], ('conflicting fresh answers',qid)
        answers[qid]=dict(trace_index=i,**row)
    missing=[r['id']for r in requests if r['id']not in answers]
    assert not missing, ('missing/unread fresh requests; do not score legacy answers',missing)
    out.mkdir(parents=True,exist_ok=False)
    scores={cid:C.Counter() for cid in hypotheses['candidate_ids']}
    per_question=[]
    signatures={cid:[]for cid in scores}
    # No survivor elimination occurs in this scoring pass.
    for line in open(packet/'predictions-native.jsonl'):  # ubs:ignore — finite one-shot iterator; fully consumed and released before process exit
        pred=json.loads(line);qid=pred['id'];expected=answers[qid]['expected']
        seen=set();matched=refused=0
        for partition in pred['partitions']:
            actual=partition['outcome']
            status='matched' if actual==expected else 'refused' if actual.startswith('fault:') else 'mismatched'
            for cid in partition['candidate_ids']:
                assert cid in scores and cid not in seen
                seen.add(cid);scores[cid][status]+=1;signatures[cid].append(actual)
            matched+=len(partition['candidate_ids'])*(status=='matched')
            refused+=len(partition['candidate_ids'])*(status=='refused')
        assert seen==scores.keys()
        per_question.append(dict(id=qid,expected=expected,matched_hypotheses=matched,
            failed_hypotheses=len(scores)-matched,refused_hypotheses=refused,
            answer=answers[qid],prediction_partitions=len(pred['partitions'])))
    assert len(per_question)==38
    with (out/'per-hypothesis-blind-scores.csv').open('x',newline='')as f:
        w=csv.writer(f);w.writerow(['candidate_id','pass','fail','mismatched','refused','missing','total','pass_rate','prediction_signature_sha256'])
        for cid,c in sorted(scores.items()):
            assert sum(c.values())==38
            w.writerow([cid,c['matched'],38-c['matched'],c['mismatched'],c['refused'],0,38,c['matched']/38,
                hashlib.sha256(json.dumps(signatures[cid],separators=(',',':')).encode()).hexdigest()])
    histogram=C.Counter(c['matched']for c in scores.values())
    summary=dict(stage='blind scoring only; no elimination performed',
        scored_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        frozen_packet_manifest_sha256=sha(packet/'manifest.json'),trace_sha256=sha(trace),
        prediction_sha256=sha(packet/'predictions-native.jsonl'),
        hypotheses=len(scores),question_count=38,score_histogram={str(k):v for k,v in sorted(histogram.items())},
        per_hypothesis_file='per-hypothesis-blind-scores.csv',
        evidence_policy='Every hypothesis scored on all requests; refuses/missing count as failures. Both backend frozen predictions have identical hashes.',
        questions=per_question)
    (out/'blind-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items()if k!='questions'}))


def eliminate(packet,scoring,out):
    requests,hypotheses,_=load_packet(packet)
    summary=json.loads((scoring/'blind-summary.json').read_text())
    assert summary['stage'].startswith('blind scoring only')
    answers={q['id']:q for q in summary['questions']}
    alive=set(hypotheses['candidate_ids']);first={};per_question=[]
    out.mkdir(parents=True,exist_ok=False)
    with (out/'first-kills.jsonl').open('x')as f:
        for line in open(packet/'predictions-native.jsonl'):  # ubs:ignore — finite one-shot iterator; fully consumed and released before process exit
            row=json.loads(line);qid=row['id'];expected=answers[qid]['expected'];before=len(alive)
            for partition in row['partitions']:
                if partition['outcome']==expected:continue
                for cid in sorted(alive.intersection(partition['candidate_ids'])):
                    record=dict(candidate_id=cid,question=qid,expected=expected,got=partition['outcome'],
                        answer=answers[qid]['answer'])
                    f.write(json.dumps(record,separators=(',',':'))+'\n');first[cid]=record
                    alive.remove(cid)
            per_question.append(dict(id=qid,killed=before-len(alive),remaining=len(alive)))
    # Reconcile with the already-written blind per-hypothesis scores.
    with (scoring/'per-hypothesis-blind-scores.csv').open()as f:
        for r in csv.DictReader(f):
            assert (int(r['candidate_id'])in alive)==(int(r['pass'])==38)
    assert len(alive)+len(first)==52908
    (out/'survivors.json').write_text(json.dumps(sorted(alive))+'\n')
    result=dict(stage='exact elimination after blind scoring',
        eliminated_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        blind_summary_sha256=sha(scoring/'blind-summary.json'),initial=52908,
        survivors=len(alive),killed=len(first),questions=per_question,
        claim='Scoped finite-family elimination only; no algorithm identification or product promotion.')
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


def score_round10(packet, evidence, out):
    """Score the immutable R10 packet against existing witness-reader exports.

    No candidate elimination or amended prediction is performed here. Literal
    observations remain transport audits, even when an echo input was unchanged.
    """
    from fractions import Fraction as Q
    import struct
    scorer_hash = sha(Path(__file__))
    frozen = {}
    for name in ('PMT-manifest.json', 'NORM-manifest.json', 'NORM-lattice-manifest.json'):
        manifest = json.loads((packet/name).read_text())
        for file, digest in manifest['files_sha256'].items():
            assert sha(packet/file) == digest, ('frozen packet changed', file)
            frozen[file] = digest
    provenance = json.loads((packet/'reader-provenance.json').read_text())
    for file, digest in provenance['application_prediction_hashes'].items():
        assert sha(packet/file) == digest, ('application prediction changed', file)
        frozen[file] = digest
    # The Windows copy was retained without opening its prediction files.
    archive = evidence/'probes/2026-11-01-calcvba-run8'
    for file, digest in frozen.items():
        assert sha(archive/'packet'/file) == digest, ('Windows packet differs', file)
    # Frozen requests/IDs are checked before any answer file is read.
    requests = {fn: json.loads((packet/file).read_text())['witnesses'] for fn,file in (
        ('PMT','requests-PMT.json'), ('NORM.S.DIST','requests-NORM.S.DIST-all.json'))}
    assert len(requests['PMT']) == 64 and len(requests['NORM.S.DIST']) == 24
    assert all(len({q['id'] for q in rows}) == len(rows) for rows in requests.values())
    manifest = json.loads((archive/'blind-manifest.json').read_text())
    for name, digest in manifest['sha256'].items():
        assert sha(archive/name) == digest, ('answer provenance hash differs', name)
    answers = {}
    source_hashes = {}
    for fn in requests:
        path = evidence/'readings/excel-20527'/(fn+'.json')
        source_hashes[str(path)] = sha(path)
        doc = json.loads(path.read_text())
        assert doc['unread'] == []
        assert doc['capture_provenance']['environment']['excel_build'] == '20527'
        semantic = {r['id']: r for r in doc['witnesses']}
        audits = {r['id']: r for r in doc['transport_audit_readings']}
        assert not (semantic.keys() & audits.keys())
        assert semantic.keys() | audits.keys() == {q['id'] for q in requests[fn]}
        for q in requests[fn]:
            literal = q['input_transport'] == 'formula_literals'
            row = (audits if literal else semantic)[q['id']]
            assert row['args'] == q['args'] and row.get('result_index') == q.get('result_index')
            if literal:
                # This is the reader's raw outcome; never a semantic-bank import.
                assert row['semantic_witness'] is False
                raw = row['observed_raw']
                assert set(raw) == {'number'}, ('unsupported audit outcome', raw)
                observed = '0x' + raw['number']
            else:
                observed = row['expected_bits']
            answers[q['id']] = dict(observed=observed, semantic=not literal,
                input_admission_failed=bool(row.get('reading_failure')), reading=row)
    out.mkdir(parents=True, exist_ok=False)
    details = []
    selection = {r['id']: r for r in json.loads((packet/'PMT-selection.json').read_text())['questions']}
    historical = {r['id']: r for r in map(json.loads,(packet/'PMT-predictions-historical.jsonl').read_text().splitlines())}
    application = {r['id']: r for r in map(json.loads,(packet/'PMT-predictions-application.jsonl').read_text().splitlines())}
    pmt_counts = C.defaultdict(C.Counter)
    for q in requests['PMT']:
        qid=q['id']; observed=answers[qid]['observed']; meta=selection[qid]
        h=historical[qid]['predicted']; a=application[qid]['calc_func_prediction']
        row=dict(id=qid, variant=meta['variant'], sample_role=meta['sample_role'],
            paired_numeric_id=meta['paired_numeric_id'], requested_args=q['args'],
            **answers[qid], historical_prediction=h, application_prediction=a,
            historical_matches=observed==h, application_matches=observed==a)
        details.append(row)
        for label, predicted in [('historical',h),('application',a)]:
            c=pmt_counts[(meta['variant'],label)]
            c['total']+=1
            c['raw_output_coincidences']+=observed==predicted
            passed=not row['input_admission_failed'] and observed==predicted
            c['pass']+=passed;c['fail']+=not passed
            c['input_admission_failed']+=row['input_admission_failed']
    normal_application={r['id']:r for r in map(json.loads,(packet/'NORM-predictions-application.jsonl').read_text().splitlines())}
    normal_counts=C.Counter(); profiles={i:C.Counter()for i in range(5)}
    for q in requests['NORM.S.DIST']:
        actual=answers[q['id']]['observed']; predicted=normal_application[q['id']]['calc_func_prediction']
        details.append(dict(id=q['id'],**answers[q['id']],application_prediction=predicted,application_matches=actual==predicted))
        normal_counts['total']+=1;normal_counts['pass']+=actual==predicted;normal_counts['fail']+=actual!=predicted
    base=list(map(json.loads,(packet/'NORM-predictions.jsonl').read_text().splitlines()))
    anchors=[]; relation_rows=[]
    for anchor,companion in zip(base[::2],base[1::2]):
        negative=answers[anchor['id']]['observed']; positive=answers[companion['id']]['observed']
        prerequisite=negative==anchor['predicted']
        anchors.append(dict(id=anchor['id'],observed=negative,predicted=anchor['predicted'],matches=prerequisite))
        for alternative in companion['alternatives']:
            pid=alternative['profile_id'];member=positive in alternative['possible_outcomes']
            c=profiles[pid];c['questions']+=1;c['failed_prerequisite']+=not prerequisite
            c['pass']+=prerequisite and member;c['fail']+=not prerequisite or not member
            relation_rows.append(dict(profile_id=pid,negative=anchor['id'],positive=companion['id'],
                prerequisite_met=prerequisite,possible_outcomes=alternative['possible_outcomes'],
                observed=positive,membership=member,passed=prerequisite and member))
    def rational_outcome(s):
        return Q(struct.unpack('>d',bytes.fromhex(s[2:]))[0])
    def contains(interval,x):
        lo,hi=(Q(*map(int,interval[k]))for k in ('lower','upper'))
        return (x>lo or x==lo and interval['lower_closed']) and (x<hi or x==hi and interval['upper_closed'])
    lattice=json.loads((packet/'NORM-lattice-predictions.json').read_text());lattice_rows=[]
    for pair in lattice['pairs']:
        negative=answers[pair['negative_request']]['observed'];positive=answers[pair['positive_request']]['observed']
        applicable=[(i,c)for i,c in enumerate(lattice['cases'])if contains(c['negative_outcome_interval'],rational_outcome(negative))]
        assert len(applicable)<=1
        for pid in profiles:
            possible=applicable[0][1]['alternatives'][pid]['possible_outcomes']if applicable else []
            passed=bool(applicable)and positive in possible;c=profiles[pid]
            c['questions']+=1;c['failed_prerequisite']+=not applicable;c['pass']+=passed;c['fail']+=not passed
            lattice_rows.append(dict(profile_id=pid,negative=pair['negative_request'],positive=pair['positive_request'],
                negative_observed=negative,positive_observed=positive,case=applicable[0][0]if applicable else None,
                prerequisite_met=bool(applicable),possible_outcomes=possible,passed=passed))
    assert sha(Path(__file__))==scorer_hash
    summary=dict(stage='blind scoring only; no elimination performed',
        scored_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scorer_sha256=scorer_hash,
        frozen_prediction_hashes=frozen,answer_source_hashes=source_hashes,
        calls=88,semantic_calls=sum(r['semantic']for r in answers.values()),literal_audits=sum(not r['semantic']for r in answers.values()),
        literal_failed_admission=sum(r['input_admission_failed']for r in answers.values()),
        PMT=[dict(variant=k[0],prediction=k[1],**v)for k,v in sorted(pmt_counts.items())],
        NORM_application=normal_counts,NORM_anchor_rechecks=anchors,
        NORM_conditional_profile_scores=profiles,
        limitations='Five normal relations are necessary shared-latent constraints, not full programs. Failed prerequisites remain failures of this frozen experiment, not semantic contradictions of an inapplicable relation. Literal pass means echo admission succeeded and raw output coincided, an audit-only success; all changed-input obligations fail even if their outputs coincide. Literal outcomes are never semantic matches or same-input semantic contradictions; effective callee inputs unavailable.')
    (out/'blind-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (out/'per-question.json').write_text(json.dumps(details,indent=2)+'\n')
    (out/'normal-relations.json').write_text(json.dumps(dict(base=relation_rows,lattice=lattice_rows),indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items()if k not in ('frozen_prediction_hashes','answer_source_hashes')}))


def score_routes(packet,trace,out):
    """Score all frozen route predictions from the existing canonical judge trace.

    Missing and unread requests fail. Kind-only boundary obligations are scored
    separately and never become exact numerical matches. No search/elimination.
    """
    import math
    import struct
    scorer_hash=sha(Path(__file__))
    manifest=json.loads((packet/'manifest.json').read_text())
    for name,digest in manifest['files_sha256'].items():
        assert sha(packet/name)==digest,('frozen packet changed',name)
    requests={};predictions={}
    for name in manifest['files_sha256']:
        if name.startswith('requests-'):
            doc=json.loads((packet/name).read_text());fn=doc['function']
            for q in doc['witnesses']:
                assert q['id']not in requests and q['input_transport']=='cell_value2_bulk'
                requests[q['id']]=dict(function=fn,**q)
        if name.startswith(('application-predictions-','registered-predictions-')):
            for line in (packet/name).read_text().splitlines():
                p=json.loads(line);assert p['id']not in predictions;predictions[p['id']]=p
    assert len(requests)==manifest['calls']and requests.keys()==predictions.keys()
    for qid,q in requests.items():assert predictions[qid]['args']==q['args']
    answers={};source_hash=sha(trace)
    with trace.open()as f:
        for index,line in enumerate(f):
            row=json.loads(line);qid=row.get('id')
            if qid not in requests:continue
            q=requests[qid]
            assert row['function']==q['function']and row['args']==q['args']and row.get('result_index')==q.get('result_index'),('input binding differs',qid)
            assert row['reading'].get('excel_build')and row['source']
            assert not row.get('not_excel_evidence'),('synthetic trace is not an answer',qid)
            if qid in answers:assert answers[qid]['expected']==row['expected'],('conflicting readings',qid)
            answers[qid]=dict(trace_index=index,**row)
    totals=C.defaultdict(C.Counter);details=[]
    for qid,q in requests.items():
        p=predictions[qid]['calc_func_prediction'];answer=answers.get(qid)
        if answer is None:status='missing'
        elif p.startswith('fault:'):status='refused'
        elif p==answer['expected']:status='matched'
        else:status='mismatched'
        c=totals[q['function']];c[status]+=1;c['total']+=1;c['pass']+=status=='matched';c['fail']+=status!='matched'
        details.append(dict(request=q,prediction=p,status=status,answer=answer))
    boundary=[];boundary_counts=C.Counter()
    path=packet/'boundary-hypotheses.json'
    if path.name in manifest['files_sha256']:
        for h in json.loads(path.read_text()):
            qid=h['id'];answer=answers.get(qid);wanted=h['kind_only_boundary_hypothesis']
            observed=answer['expected']if answer else None
            finite=(isinstance(observed,str)and observed.startswith('0x')and len(observed)==18 and
                math.isfinite(struct.unpack('>d',bytes.fromhex(observed[2:]))[0]))
            passed=observed is not None and (finite if wanted.startswith('finite number')else observed==wanted)
            boundary.append(dict(hypothesis=h,observed=observed,passed=passed));boundary_counts['total']+=1;boundary_counts['pass']+=passed;boundary_counts['fail']+=not passed
    assert sha(trace)==source_hash  # ubs:ignore — public evidence identity, no secret comparison
    assert sha(Path(__file__))==scorer_hash  # ubs:ignore — public source identity, no secret comparison
    out.mkdir(parents=True,exist_ok=False)
    summary=dict(stage='blind route scoring only; no search or elimination',
        scored_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scorer_sha256=scorer_hash,
        packet_manifest_sha256=sha(packet/'manifest.json'),trace_sha256=source_hash,
        calls=len(requests),available_readings=len(answers),functions=dict(totals),
        boundary_kind_obligations=boundary_counts,
        limits='Canonical existing judge trace only. Missing/unread obligations fail. Kind-only outcomes remain separate from exact value scores. No effective-callee-input observation is inferred. Evaluate frozen coupled relations only after retaining this blind score.')
    (out/'blind-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (out/'per-question.json').write_text(json.dumps(details,indent=2)+'\n')
    (out/'boundary-kind-obligations.json').write_text(json.dumps(boundary,indent=2)+'\n')
    print(json.dumps(summary))


if __name__=='__main__':
    action,packet,source,out=sys.argv[1:]
    {'score':score,'eliminate':eliminate,'score-round10':score_round10,'score-routes':score_routes}[action](Path(packet),Path(source),Path(out))
