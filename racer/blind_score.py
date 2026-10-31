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
    for i,line in enumerate(open(trace)):
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
    for line in open(packet/'predictions-native.jsonl'):
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
        for line in open(packet/'predictions-native.jsonl'):
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


if __name__=='__main__':
    action,packet,source,out=sys.argv[1:]
    {'score':score,'eliminate':eliminate}[action](Path(packet),Path(source),Path(out))
