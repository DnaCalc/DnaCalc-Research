"""Deliberately corrupted R13 proofs/packets must never receive a receipt."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

from x87 import check
from questions_r13 import check as check_questions, check_readings, check_annuity_questions


def shadow(source,destination,replaced):
    destination.mkdir(parents=True,exist_ok=False)
    for src in source.rglob('*'):
        relative=src.relative_to(source);dst=destination/relative
        if src.is_dir():dst.mkdir(exist_ok=True)
        elif relative.as_posix()not in replaced:dst.symlink_to(src.resolve())


def run(proof,packet,trace,gaps,library,output):
    proof,packet,output=map(Path,(proof,packet,output))
    output.mkdir(parents=True,exist_ok=False)
    results=[]
    def rejected(label,function,message,receipt):
        try:function()
        except ValueError as error:
            if message not in str(error):raise
            if receipt.exists():raise ValueError('failed control wrote a successful receipt')
            results.append(dict(control=label,rejected=True,error=str(error),success_receipt_exists=False))
        else:raise ValueError('corrupted control was accepted: '+label)
    bad=output/'bad-helper-output';shadow(proof,bad,{'rows/000.json'})
    with (proof/'rows/000.json').open('rb')as src,(bad/'rows/000.json').open('xb')as dst:
        first=src.read(1<<20);marker=b'"values": "';at=first.index(marker)+len(marker)
        first=first[:at]+(b'A'if first[at:at+1]!=b'A'else b'B')+first[at+1:]
        dst.write(first);shutil.copyfileobj(src,dst)
    receipt=output/'bad-helper-receipt.json'
    rejected('changed one native helper output byte',lambda:check(bad,trace,library,receipt),
             'helper instruction/store output differs',receipt)
    bad=output/'missing-first-kill';shadow(proof,bad,{'first-kills.jsonl'})
    with (proof/'first-kills.jsonl').open()as src,(bad/'first-kills.jsonl').open('x')as dst:
        next(src);shutil.copyfileobj(src,dst)
    receipt=output/'missing-kill-receipt.json'
    rejected('omitted one first-kill record',lambda:check(bad,trace,library,receipt),
             'missing representative first-kill trace',receipt)
    bad=output/'wrong-function-graph';shutil.copytree(packet,bad)
    path=bad/'predictions/native-elementary.jsonl'
    rows=[json.loads(line)for line in path.read_text().splitlines()]
    first=rows[0];first['graph']=dict(nodes=[dict(op='input',a=0)],outputs=[0,0],sites=0)
    request=next(r for r in json.loads((bad/'requests-LN.json').read_text())['witnesses']if r['id']==first['id'])
    first['predictions_by_store_mask']=[request['args'][0]]
    path.write_text(''.join(json.dumps(row)+'\n'for row in rows))
    manifest=json.loads((bad/'manifest.json').read_text())
    manifest['files_sha256']['predictions/native-elementary.jsonl']=hashlib.sha256(path.read_bytes()).hexdigest()
    (bad/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    receipt=output/'wrong-function-receipt.json'
    rejected('self-consistent wrong function graph with updated checksum',
             lambda:check_questions(bad,trace,gaps,library,receipt),
             'standalone function graph differs',receipt)
    result=dict(complete=True,controls=results,all_rejected=True,
                checker_source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in (
                    Path(__file__),Path(__file__).with_name('x87.py'),Path(__file__).with_name('questions_r13.py'),Path(__file__).with_name('exact.py'))},
                input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in (
                    proof/'summary.json',packet/'manifest.json',Path(trace),Path(gaps),Path(library))},
                limits='Constructed corruptions, never Excel observations. Original proof and packet bytes are unchanged; shadow copies/symlinks retain failures.')
    (output/'controls.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result));return result


def intake_controls(packet, results, output):
    """A changed value, changed input or omitted reading batch must be rejected."""
    packet,results,output=map(Path,(packet,results,output));output.mkdir(parents=True,exist_ok=False)
    controls=[]
    function = 'VAR.S' if (packet/'requests-VAR.S.json').is_file() else 'CUMIPMT'
    def run(label,change,message):
        shadow=output/label;shutil.copytree(results,shadow);change(shadow)
        receipt=output/(label+'-must-not-exist.json')
        try:check_readings(packet,shadow,receipt)
        except ValueError as error:
            if message not in str(error):raise
            if receipt.exists():raise ValueError('corrupted intake received a success receipt')
            controls.append(dict(control=label,rejected=True,error=str(error),success_receipt_exists=False))
        else:raise ValueError('corrupted intake was accepted')
    def changed_answer(root):
        path=root/f'OxFunc-Evidence/readings/excel-20527/{function}.json'
        doc=json.loads(path.read_text());doc['witnesses'][0]['expected_bits']='0x3ff0000000000000'
        path.write_text(json.dumps(doc))
    def changed_input(root):
        path=next((root/'worker/read').glob(function+'-*/result.json'))
        doc=json.loads(path.read_text());doc['rows'][0]['input_audit']['read_back'][0]['number']='3ff0000000000000'
        path.write_text(json.dumps(doc))
    def missing_batch(root):
        (root/f'OxFunc-Evidence/readings/excel-20527/{function}.json').unlink()
    def wrong_formula(root):
        request=next((root/'worker/read').glob(function+'-*/request.json'))
        result=request.with_name('result.json')
        asked=json.loads(request.read_text());readings=json.loads(result.read_text())
        asked['rows'][0]['formula']='=PMT(RC[1],RC[2],RC[3])'
        readings['rows'][0]['input_audit']['formula']=asked['rows'][0]['formula']
        request.write_text(json.dumps(asked));result.write_text(json.dumps(readings))
    run('changed-canonical-answer',changed_answer,'reader/raw result differs')
    run('changed-readback-input',changed_input,'transport changed requested operands')
    run('omitted-reader-batch',missing_batch,'missing reader batch')
    run('self-consistent-wrong-function',wrong_formula,'requested function/formula differs')
    result=dict(complete=True,all_rejected=True,controls=controls,
        checker_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in (
            Path(__file__),Path(__file__).with_name('questions_r13.py'))},
        limits='Corrupted private copies only; original packet and supplied result files unchanged.')
    (output/'controls.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


def next_packet_controls(packet, trace, output):
    packet,trace,output=map(Path,(packet,trace,output));output.mkdir(parents=True,exist_ok=False)
    controls=[]
    for label,kind,message in [('wrong-count','coupled payment/principal/interest preimages','relation count differs'),
            ('wrong-loan','published payment partition audit','partition loan/period differs')]:
        shadow=output/label;shutil.copytree(packet,shadow)
        path=shadow/'conditional-relations.json';relations=json.loads(path.read_text())
        relation=next(r for r in relations if r['kind']==kind)
        if label=='wrong-count':relation['count']+=1
        else:relation['PMT']='r15-pmt-001'
        path.write_text(json.dumps(relations,indent=2)+'\n')
        manifest=shadow/'manifest.json';doc=json.loads(manifest.read_text())
        doc['files_sha256'][path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        manifest.write_text(json.dumps(doc,indent=2)+'\n')
        receipt=output/(label+'-must-not-exist.json')
        try:check_annuity_questions(shadow,trace,receipt)
        except ValueError as error:
            if message not in str(error):raise
            if receipt.exists():raise ValueError('forged relation received a success receipt')
            controls.append(dict(control=label,rejected=True,error=str(error),success_receipt_exists=False))
        else:raise ValueError('forged relation accepted after checksum update')
    result=dict(complete=True,all_rejected=True,controls=controls,
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in (
            Path(__file__),Path(__file__).with_name('questions_r13.py'))},
        limits='Private forged packet copies with updated checksums; the original packet is unchanged.')
    (output/'controls.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


def round16_controls(artifacts, packet, readings, workbook, output):
    """Self-consistent forgeries must fail semantic checks, not just hashes."""
    from geomean import check_exact_roots, verify_static_reference
    from regression import check_rank
    from bessel import check_domains
    from questions_r13 import check_round16_questions
    from variance import check as check_variance
    import zipfile
    artifacts,packet,output=map(Path,(artifacts,packet,output));output.mkdir()
    results=[]
    def rejected(label,function,message):
        receipt=output/(label+'-must-not-exist.json')
        try:function(receipt)
        except ValueError as error:
            if message not in str(error):raise
            if receipt.exists():raise ValueError('forged control received a receipt')
            results.append(dict(control=label,rejected=True,error=str(error),success_receipt_exists=False))
        else:raise ValueError('forged control was accepted: '+label)
    def read_lines(path):return [json.loads(s)for s in path.read_text().splitlines()]
    def lines(path,rows):path.write_text(''.join(json.dumps(r)+'\n'for r in rows))
    source=artifacts/'geomean-root-references'
    bad=output/'wrong-root';shadow(source,bad,{'first-kills.jsonl'})
    rows=read_lines(source/'first-kills.jsonl');rows[0]['got']=f"0x{int(rows[0]['got'],16)^1:016x}";lines(bad/'first-kills.jsonl',rows)
    rejected('wrong-root',lambda receipt:check_exact_roots(artifacts/'geomean-with-references.jsonl',bad,receipt),'incorrect product/root')
    bad2=output/'missing-root-kill';shadow(source,bad2,{'first-kills.jsonl'});lines(bad2/'first-kills.jsonl',rows[1:])
    rejected('missing-root-kill',lambda receipt:check_exact_roots(artifacts/'geomean-with-references.jsonl',bad2,receipt),'unsearched root program')
    # Pair a different package with the unchanged original capture, updating all
    # proposed operand and package hashes. The capture lineage must still reject.
    row=json.loads((artifacts/'geomean-references/cases.jsonl').read_text().splitlines()[0])
    binding=row['reference_binding'];original=Path(binding['package']);forged=output/'unrelated-package.xlsx'
    with zipfile.ZipFile(original)as src,zipfile.ZipFile(forged,'w')as dst:
        for info in src.infolist():
            data=src.read(info.filename)
            if info.filename=='xl/worksheets/sheet1.xml':
                data=data.replace(b'<c r="A1"><v>1</v></c>',b'<c r="A1"><v>2</v></c>')
            dst.writestr(info,data)
    binding['package']=str(forged);binding['package_sha256']=hashlib.sha256(forged.read_bytes()).hexdigest()
    binding['values'][0]='0x4000000000000000';row['args'][0][0][0]='0x4000000000000000'
    rejected('unrelated-reference-package',lambda receipt:verify_static_reference(row),'capture/package lineage differs')
    source=artifacts/'rank-geometry';bad3=output/'false-rank-survivor'
    shadow(source,bad3,{'first-kills.jsonl','survivors.json'})
    rows=read_lines(source/'first-kills.jsonl');identity=rows[0]['candidate_id']
    lines(bad3/'first-kills.jsonl',rows[1:]);alive=json.loads((source/'survivors.json').read_text());alive.append(identity)
    (bad3/'survivors.json').write_text(json.dumps(alive)+'\n')
    rejected('false-rank-survivor',lambda receipt:check_rank(artifacts/'all-other/constraints.jsonl',bad3,receipt),'invalid rank survivor')
    source=artifacts/'bessel-domains';bad4=output/'forged-domain-count';shadow(source,bad4,{'results.json'})
    rows=json.loads((source/'results.json').read_text());target=next(r for r in rows if r['rejects_recorded_numeric'])
    target['rejects_recorded_numeric']=0;target['first_early_contradiction']=None
    (bad4/'results.json').write_text(json.dumps(rows)+'\n')
    rejected('forged-domain-count',lambda receipt:check_domains(artifacts/'all-other-v2/constraints.jsonl',bad4,receipt),'incorrect domain contradiction')
    source=artifacts/'all-AVERAGE-native';bad5=output/'empty-variance';shadow(source,bad5,{'summary.json'})
    summary=json.loads((source/'summary.json').read_text());summary.update(observations=0,executions=0,survivors=75)
    (bad5/'summary.json').write_text(json.dumps(summary)+'\n');empty=output/'empty.jsonl';empty.write_text('')
    rejected('empty-variance',lambda receipt:check_variance(empty,bad5,receipt),'empty variance observation scope')
    bad6=output/'unknown-variance-function';shadow(source,bad6,{'summary.json'})
    summary['function']='SYNTHETIC_UNKNOWN';(bad6/'summary.json').write_text(json.dumps(summary)+'\n')
    rejected('unknown-variance-function',lambda receipt:check_variance(empty,bad6,receipt),'unknown variance function')
    def rehash(root):
        m=json.loads((root/'manifest.json').read_text())
        m['files_sha256']={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()for p in root.rglob('*')if p.is_file()and p.name!='manifest.json'}
        m['catalogue_section_sha256']=hashlib.sha256((root/'sheet-dcj9k_r16.txt').read_bytes()).hexdigest()
        (root/'manifest.json').write_text(json.dumps(m,indent=2)+'\n')
    bad7=output/'wrong-packet-input';shutil.copytree(packet,bad7)
    path=bad7/'requests-VAR.S.json';doc=json.loads(path.read_text());qid=doc['witnesses'][0]['id'];doc['witnesses'][0]['args'][0]='0x0000000000000000';path.write_text(json.dumps(doc)+'\n')
    path=bad7/'application-predictions-VAR.S.jsonl';rows=read_lines(path);rows[0]['args'][0]='0x0000000000000000';lines(path,rows)
    path=bad7/'selection.json';rows=json.loads(path.read_text());next(r for r in rows if r['id']==qid)['requested_args'][0]='0x0000000000000000';path.write_text(json.dumps(rows)+'\n');rehash(bad7)
    rejected('wrong-packet-input',lambda receipt:check_round16_questions(bad7,readings,workbook,receipt),'scalar derivation/count differs')
    bad8=output/'wrong-rank-question';shutil.copytree(packet,bad8)
    path=bad8/'catalogue-questions.json';doc=json.loads(path.read_text());q=next(q for q in doc['questions']if q['reason']=='rank and removed-column audit');q['formula']=q['formula'].replace('LINEST(','LOGEST(');at=q['row'];replacement=q['formula'];path.write_text(json.dumps(doc)+'\n')
    path=bad8/'sheet-dcj9k_r16.txt';text=path.read_text().splitlines();text[at]=replacement;path.write_text('\n'.join(text)+'\n')
    path=bad8/'application-predictions-catalogue.jsonl';rows=read_lines(path);rows[at-1]['formula']=replacement;lines(path,rows);rehash(bad8)
    rejected('wrong-rank-question',lambda receipt:check_round16_questions(bad8,readings,workbook,receipt),'catalogue derivation changed')
    result=dict(complete=True,all_rejected=True,controls=results,
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in (
            Path(__file__),Path(__file__).with_name('geomean.py'),Path(__file__).with_name('regression.py'),
            Path(__file__).with_name('bessel.py'),Path(__file__).with_name('variance.py'),Path(__file__).with_name('questions_r13.py'))},
        limits='Constructed corruptions in private copies only; original captures, packages, answers, search outputs and question packet unchanged. Recomputed hashes do not authorize changed semantics.')
    (output/'controls.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='intake':intake_controls(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='next-packet':next_packet_controls(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='round16':round16_controls(*sys.argv[2:]);raise SystemExit(0)
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('proof','packet','trace','gaps','library','output'):p.add_argument(key)
    a=p.parse_args();run(a.proof,a.packet,a.trace,a.gaps,a.library,a.output)
