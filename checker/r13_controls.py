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
        path=root/'OxFunc-Evidence/readings/excel-20527/CUMIPMT.json'
        doc=json.loads(path.read_text());doc['witnesses'][0]['expected_bits']='0x3ff0000000000000'
        path.write_text(json.dumps(doc))
    def changed_input(root):
        path=next((root/'worker/read').glob('CUMIPMT-*/result.json'))
        doc=json.loads(path.read_text());doc['rows'][0]['input_audit']['read_back'][0]['number']='3ff0000000000000'
        path.write_text(json.dumps(doc))
    def missing_batch(root):
        (root/'OxFunc-Evidence/readings/excel-20527/CUMIPMT.json').unlink()
    def wrong_formula(root):
        request=next((root/'worker/read').glob('CUMIPMT-*/request.json'))
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


if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='intake':intake_controls(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='next-packet':next_packet_controls(*sys.argv[2:]);raise SystemExit(0)
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('proof','packet','trace','gaps','library','output'):p.add_argument(key)
    a=p.parse_args();run(a.proof,a.packet,a.trace,a.gaps,a.library,a.output)
