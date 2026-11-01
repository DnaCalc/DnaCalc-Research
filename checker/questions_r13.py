"""Independent preservation, input-derivation and prediction audit for R13."""
import argparse
from collections import Counter
import csv
from pathlib import Path
import hashlib
import json

from exact import BINARY64 as B, EXTENDED64 as E, binary64
from x87 import Exact, Instructions, bits64, decode, encode, require


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def records(path):return [json.loads(line)for line in Path(path).read_text().splitlines()]


def score_readings(packet, results, executable, output):
    """Recalculate verified WitnessSets through the existing registered route."""
    import subprocess
    packet, results, executable, output = map(Path, (packet, results, executable, output))
    executable = executable.resolve()
    output.mkdir(parents=True, exist_ok=False)
    binary_hash = sha(executable)
    intake = check_readings(packet, results, output/'intake.json')
    rows, counts = [], Counter()
    for i, batch in enumerate(intake['batches']):
        path = results/batch['reading']
        require(sha(path) == batch['reading_sha256'], 'reading changed before scoring')
        doc = read(path)
        prediction = output/f'prediction-{i:02d}.jsonl'
        subprocess.run([str(executable), 'invoke-inputs', str(path), str(prediction)], check=True,
                       timeout=300)  # ubs:ignore — explicit compiled executable and argv, no shell; binary hash bound below
        actuals = records(prediction)
        require(len(actuals) == len(doc['witnesses']), 'missing calculated result')
        for original, actual in zip(doc['witnesses'], actuals):
            require(original['id'] == actual['id'] and original['args'] == actual['args'],
                    'calculated identity differs')
            env = doc['capture_provenance']['environment']
            rows.append(dict(id=original['id'], function=doc['function'], args=original['args'],
                result_index=original.get('result_index'), expected=original['expected_bits'],
                actual=actual['calc_func_prediction'],
                reading=dict(bitness=env['excel_bitness'], cpu=env['cpu_id'],
                             excel_build=env['excel_build'], platform=env['platform']),
                source='readings/excel-20527/r13-'+path.name,
                capture_provenance=doc['capture_provenance'], source_sha256=batch['reading_sha256']))
            status = 'matched' if rows[-1]['expected'] == rows[-1]['actual'] else 'failed'
            counts[status] += 1
        require(sha(path) == batch['reading_sha256'], 'reading changed during scoring')
    require(sha(executable) == binary_hash and len(rows) == 234, 'scoring execution changed')
    (output/'observations.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows))
    functions = {}
    for function in sorted({r['function'] for r in rows}):
        group = [r for r in rows if r['function'] == function]
        passed = sum(r['actual'] == r['expected'] for r in group)
        functions[function] = dict(total=len(group), matched=passed, failed=len(group)-passed)
    result = dict(complete=True, total=len(rows), **counts, functions=functions,
        executable_sha256=binary_hash, observations_sha256=sha(output/'observations.jsonl'),
        intake_sha256=sha(output/'intake.json'),
        limits='Every supplied observation is scored, including repeated inputs. Existing '
               'registered invoke route; typed failures count as failures. Value2 read-back '
               'does not expose effective callee inputs. Predictions never become answers.')
    (output/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


def check_readings(packet, results, output):
    """Verify the existing reader exports and raw transport; never create answers."""
    packet, results = Path(packet), Path(results)
    received = results/'packet'
    manifest = read(packet/'manifest.json')
    files = {p.relative_to(packet).as_posix(): sha(p)
             for p in packet.rglob('*') if p.is_file()}
    require(files == read(results/'packet-sha256.json'), 'received packet manifest differs')
    require(files == {p.relative_to(received).as_posix(): sha(p)
                      for p in received.rglob('*') if p.is_file()}, 'received packet bytes differ')
    for name, expected in manifest['files_sha256'].items():
        require(files[name] == expected, 'original packet hash differs: '+name)
    # Preserve/hash the whole supplied result tree, including predictions. Only
    # canonical witness-reader exports, not those predictions, enter the bank.
    frozen = {p.relative_to(results).as_posix(): sha(p)
              for p in results.rglob('*') if p.is_file()}
    exports = {}
    for path in (results/'OxFunc-Evidence/readings/excel-20527').glob('*.json'):
        doc = read(path)
        source = doc['capture_provenance']['frozen_model_predictions']['inputs_sha256']
        require(source not in exports, 'duplicate batch export')
        exports[source] = path, doc
    batches, seen, total, predictions = [], set(), 0, Counter()
    for batch in read(packet/'batches.json'):
        request_path = packet/batch['path']
        request = read(request_path)
        request_hash = sha(request_path)
        require(request_hash in exports, 'missing reader batch: '+batch['path'])
        path, doc = exports.pop(request_hash)
        function = batch['function']
        require(request['function'] == doc['function'] == function, 'function differs')
        require(not doc['unread'] and not doc['transport_audit_readings'], 'unread/audit obligation')
        provenance = doc['capture_provenance']
        model = provenance['frozen_model_predictions']
        model_name = model['frozen_path'].replace('\\', '/').split('/')[-1]
        model_path = results/'worker/frozen'/model_name
        require(sha(model_path) == model['frozen_sha256'], 'frozen worker prediction changed')
        worker = results/'worker/read'/(model_name[:-5]+'-1')
        raw, submitted = read(worker/'result.json'), read(worker/'request.json')
        require(raw['fatal'] is None and not raw['failures'], 'worker failed')
        require(submitted['function'] == function, 'submitted function differs')
        require(raw['automation_security_before'] == raw['automation_security_after'] == 1,
                'worker security settings changed')
        require(raw['excel_file_version'] == '16.0.20527.20088' and
                raw['excel_machine'] == 'x86-64' and 'Ryzen 7 4800H' in raw['cpu'],
                'unexpected capture environment')
        settings = raw['calculation_settings']
        require(settings['calculation'] == -4135 and settings['PrecisionAsDisplayed'] is False
                and settings['Date1904'] is False and settings['CompatibilityVersion'] == 2,
                'unexpected calculation settings')
        require(provenance['calculation_settings'] == settings, 'settings provenance differs')
        require(provenance['environment']['excel_file_version'] == raw['excel_file_version']
                and provenance['environment']['cpu_id'] == raw['cpu'], 'environment provenance differs')
        count = batch['calls']
        require(all(len(rows) == count for rows in (request['witnesses'], doc['witnesses'],
                doc['input_audits'], submitted['rows'], raw['rows'])), 'batch coverage differs')
        for asked, canonical, audit, call, answer in zip(request['witnesses'], doc['witnesses'],
                doc['input_audits'], submitted['rows'], raw['rows']):
            qid = asked['id']
            require(qid not in seen, 'duplicate packet ID'); seen.add(qid)
            require(qid == canonical['id'] == audit['id'] == call['id'] == answer['id'],
                    'request/answer order differs')
            require(asked['args'] == canonical['args'] and
                    asked.get('result_index') == canonical.get('result_index'), 'typed identity differs')
            # This finite packet contains only explicit binary64 numbers.
            require(all(isinstance(x, str) and len(x) == 18 and x.startswith('0x')
                        for x in asked['args']), 'unexpected packet argument type')
            arguments = [dict(number=x[2:]) for x in asked['args']]
            observed = answer['input_audit']
            formula = '='+function+'('+','.join(f'RC[{i}]'for i in range(1,len(arguments)+1))+')'
            require(call['formula']==formula and call['argument_literals'] is None and
                    observed['argument_literals'] is None,'requested function/formula differs')
            require(call['args'] == arguments and observed['requested'] == arguments and
                    observed['submitted_value2'] == arguments and observed['read_back'] == arguments,
                    'transport changed requested operands')
            require(answer['input_kept'] is True and audit['reading_failure'] is None and
                    audit['audit'] == observed and audit['observed_raw'] == answer['value'],
                    'canonical/raw audit differs')
            require(observed['transport'] == 'cell_value2_bulk' and
                    observed['formula'] == call['formula'] and
                    observed['effective_callee_inputs'] == 'not directly observable',
                    'transport surface differs')
            # Audit the existing reader's encoding without importing a second bank.
            value = answer['value']
            if set(value) == {'number'}:
                expected = '0x'+value['number']
            elif set(value) == {'error_code'} and value['error_code'] == 2036:
                expected = 'error:Num'
            elif set(value) == {'text'}:
                expected = 'text'+':'+bytes.fromhex(value['text']).decode('utf-16-be')
            else:
                raise ValueError('unexpected raw outcome in finite packet')
            require(canonical['expected_bits'] == expected, 'reader/raw result differs')
            predictions['total'] += 1
            predictions['matched' if canonical['calc_func_prediction'] == expected else 'failed'] += 1
        total += count
        batches.append(dict(**batch, reading=path.relative_to(results).as_posix(),
                            reading_sha256=sha(path), request_sha256=sha(request_path)))
    require(not exports and total == 234 and len(batches) == 13, 'missing/excess batch')
    require(frozen == {p.relative_to(results).as_posix(): sha(p)
                      for p in results.rglob('*') if p.is_file()}, 'result tree changed')
    result = dict(complete=True, calls=total, batches=batches, result_files_sha256=frozen,
                  packet_files_sha256=files, frozen_windows_application=dict(predictions),
                  limits='Existing reader exports only; Value2 transport audited. Effective callee '
                         'operands remain unobservable. Native CPU measurements and predictions '
                         'are never Excel answers. Repeated captures are observations, not distinct inputs.')
    with Path(output).open('x') as f:
        json.dump(result, f, indent=2); f.write('\n')
    return result


def check_annuity_questions(packet, trace, output):
    """Independently bind the next nine-function packet to original observations."""
    import struct
    packet, trace = Path(packet), Path(trace)
    manifest=read(packet/'manifest.json')
    require(sha(trace)==manifest['source_trace_sha256'],'source observations changed')
    for name, expected in manifest['files_sha256'].items():
        require(sha(packet/name)==expected,'question file changed: '+name)
    originals={r['id']:r for r in records(trace)}
    anchors=['r11-pmt-00','r11-pmt-02','r12-pmt-group-00','r12-pmt-group-05']
    require(manifest['anchors']==anchors,'anchor selection changed')
    bits=lambda x:'0x'+struct.pack('>d',x).hex()
    numeric=lambda x:struct.unpack('>d',bytes.fromhex(x[2:]))[0]
    z,o,neg=bits(0.0),bits(1.0),bits(-1.0)
    requests={};functions=Counter()
    for batch in read(packet/'batches.json'):
        doc=read(packet/batch['path']);f=doc['function']
        require(batch['function']==f and batch['calls']==len(doc['witnesses']),'batch coverage differs')
        predictions=records(packet/f'application-predictions-{f}.jsonl')
        require(len(predictions)==batch['calls'],'prediction count differs')
        for row,prediction in zip(doc['witnesses'],predictions):
            require(set(row)=={'id','args','input_transport'} and row['input_transport']=='cell_value2_bulk',
                'a request contains an answer or different transport')
            require(row['id'] not in requests,'duplicate request ID')
            require(prediction['id']==row['id'] and prediction['args']==row['args'] and
                prediction['not_excel_evidence'] is True,'prediction identity/evidence status differs')
            requests[row['id']]=(f,row['args']);functions[f]+=1
    require(dict(functions)==manifest['calls_by_function'] and len(requests)==manifest['calls']==236,
        'missing/excess question')
    require(set(functions)=={'PMT','IPMT','PPMT','CUMPRINC','CUMIPMT','FV','PV','NPER','RATE'},
        'annuity surface missing')
    seen=set()
    for selection in read(packet/'selection.json'):
        qid=selection['id'];f,args=requests[qid];seen.add(qid)
        anchor=selection['anchor_id'];require(anchor in anchors,'unplanned anchor')
        source=originals[anchor];r,n,pv,fv,_=source['args'];t=selection['timing']
        require(source['function']=='PMT'and fv==z and t in (0,1),'anchor shape differs')
        tb=bits(float(t));role=selection['role']
        require(selection['function']==f,'selection function differs')
        if role=='cashflow scaling control':
            require(f=='PMT'and selection['scale']in (1,2,-1),'undeclared cashflow scale')
            expected=[r,n,bits(numeric(pv)*selection['scale']),z,tb]
        elif role in ('standalone interest','standalone principal'):
            require(f==('IPMT'if role.endswith('interest')else'PPMT'),'standalone function differs')
            k=selection['period'];require(1<=k<=numeric(n),'period outside loan')
            expected=[r,bits(float(k)),n,pv,z,tb]
        elif role in ('cumulative principal','cumulative interest'):
            require(f==('CUMPRINC'if role.endswith('principal')else'CUMIPMT'),'cumulative function differs')
            a,b=selection['start'],selection['end'];require(1<=a<=b<=numeric(n),'interval outside loan')
            expected=[r,n,pv,bits(float(a)),bits(float(b)),tb]
        elif role in ('unit growth or discount control','unit annuity control','adjacent period growth control'):
            require(f in ('FV','PV'),'wrong unit control function')
            period=bits(numeric(n)-1)if role.startswith('adjacent')else n
            require(numeric(period)>=1,'invalid adjacent period')
            expected=[r,period,neg if role=='unit annuity control'else z,
                      z if role=='unit annuity control'else neg,tb]
        else:
            payment=originals[selection['payment_source_id']]
            require(payment['function']=='PMT'and payment['args']==[r,n,pv,z,tb],
                'published-payment operand has a different source')
            p=payment['expected'];binary64(p)
            if role=='published-payment cross function control':
                require(f in ('FV','PV'),'wrong cross function control')
                expected=[r,n,p,pv if f=='FV'else z,tb]
            elif role=='inverse with original published payment':
                require(f=='NPER','wrong period inverse control');expected=[r,p,pv,z,tb]
            else:
                require(role=='explicit guess sensitivity control'and f=='RATE','unknown transformation')
                require(selection['guess']in ('original rate','mathematical 0.1'),'unknown guess')
                expected=[n,p,pv,z,tb,r if selection['guess']=='original rate'else bits(0.1)]
        require(args==expected,'question derivation differs: '+qid)
    require(seen==set(requests),'question lacks original-input derivation')
    for relation in read(packet/'conditional-relations.json'):
        for key in ('PMT','IPMT','PPMT','CUMPRINC','CUMIPMT'):
            if key in relation:require(requests[relation[key]][0]==key,'relation has wrong function')
        kind=relation['kind']
        if kind=='published payment partition audit':
            p=requests[relation['PMT']][1];i=requests[relation['IPMT']][1];q=requests[relation['PPMT']][1]
            require(i==q and p==[i[0],i[2],i[3],i[4],i[5]],'payment partition loan/period differs')
        elif kind=='coupled payment/principal/interest preimages':
            p=requests[relation['PMT']][1];c=requests[relation['CUMPRINC']][1];i=requests[relation['CUMIPMT']][1]
            require(c==i and p==c[:3]+[z,c[5]],'cumulative relation loan/interval differs')
            require(relation['count']==numeric(c[4])-numeric(c[3])+1,'cumulative relation count differs')
        else:
            require(kind=='singleton surface audit','unknown conditional relation')
            c=requests[relation['CUMPRINC']][1];i=requests[relation['CUMIPMT']][1]
            p=requests[relation['PPMT']][1];q=requests[relation['IPMT']][1]
            require(c==i and p==q and c[3]==c[4]==p[1]and
                    c[:3]==[p[0],p[2],p[3]]and p[4]==z and c[5]==p[5],
                    'singleton surface loan/period differs')
    result=dict(complete=True,calls=236,functions=dict(functions),anchors=anchors,new_excel_answers=0,
        packet_manifest_sha256=sha(packet/'manifest.json'),trace_sha256=sha(trace),
        checker_sha256=sha(Path(__file__)),
        limits='Original numeric arguments and fixed algebraic transformations only. Published '
               'PMT used as an authored inverse-call operand, never a fitted private helper. '
               'Predictions remain separate. PowerShell reader execution awaits Windows.')
    with Path(output).open('x')as file:json.dump(result,file,indent=2);file.write('\n')
    return result


def check_rate_newton(trace_directory, directory, output):
    """Independent rational replay of the existing bounded one-step RATE search."""
    from fractions import Fraction
    from exact import operation
    directory, trace_directory = Path(directory), Path(trace_directory)
    frozen={str(p):sha(p)for p in (Path(__file__),trace_directory/'RATE.jsonl',
        trace_directory/'FV.jsonl',directory/'grammar.json',directory/'observations.json',
        directory/'first-kills.jsonl',directory/'survivors.json',directory/'summary.json')}
    models=['strict','stored_pc64','resident_pc64','stored_pc53','resident_pc53']
    grammar=read(directory/'grammar.json')
    require(grammar['programs']==625 and grammar['models']==models and grammar['sites']==[
        'derivative PV+timing*payment','residual target_FV-observed_FV(guess)',
        'residual/derivative','guess-step'],'RATE grammar differs')
    futures={json.dumps([r['args'],r['reading']],sort_keys=True):r for r in records(trace_directory/'FV.jsonl')}
    expected_cases=[]
    for r in records(trace_directory/'RATE.jsonl'):
        if len(r['args'])!=6:continue
        try:a=list(map(binary64,r['args']));target=binary64(r['expected'])
        except (ValueError,TypeError):continue
        if a[0]!=1 or a[2]+(a[1]if int(a[4])!=0 else 0)==0:continue
        key=json.dumps([[r['args'][5],r['args'][0],r['args'][1],r['args'][2],r['args'][4]],r['reading']],sort_keys=True)
        f=futures.get(key)
        if f is None:continue
        try:binary64(f['expected'])
        except (ValueError,TypeError):continue
        expected_cases.append(dict(RATE=r,FV=f))
    require(expected_cases and read(directory/'observations.json')==expected_cases,'RATE/FV observation coverage differs')
    actual=records(directory/'first-kills.jsonl');position=0;alive=set(range(625));nodes=0
    for index,case in enumerate(expected_cases):
        r,f=case['RATE'],case['FV'];a=list(map(binary64,r['args']))
        for code in sorted(alive):
            digits=[(code//5**i)%5 for i in range(4)]
            d=operation(a[2],a[1]if int(a[4])!=0 else Fraction(0),'add',models[digits[0]])
            residual=operation(a[3],-binary64(f['expected']),'add',models[digits[1]])
            step=operation(residual,d,'div',models[digits[2]])
            y=operation(a[5],-step,'add',models[digits[3]]);published=B.round(y);nodes+=4
            require(published is not None,'RATE control produced nonfinite arithmetic')
            if published==binary64(r['expected']):continue
            require(position<len(actual),'missing RATE contradiction');record=actual[position];position+=1
            require(record['program']==code and record['models']==digits and record['case']==index
                and record['original']==case and record['expected']==r['expected']
                and binary64(record['got'])==published,'RATE first-kill identity differs')
            require(record['exact_intermediates']=={k:[str(v.numerator),str(v.denominator)]
                for k,v in [('derivative',d),('residual',residual),('step',step),('next',y)]},
                'RATE arithmetic differs')
            alive.remove(code)
    require(position==len(actual)and read(directory/'survivors.json')==sorted(alive),'RATE program accounting differs')
    summary=read(directory/'summary.json')
    require(summary['programs']==625 and summary['cases']==len(expected_cases)and summary['survivors']==len(alive),
        'RATE summary differs')
    require(frozen=={p:sha(p)for p in frozen},'RATE source/evidence changed')
    result=dict(complete=True,programs=625,cases=len(expected_cases),survivors=len(alive),
        exact_nodes_checked=nodes,input_sha256=frozen,
        limits='One analytic-derivative Newton step with the recorded FV(guess) publication '
               'as a constrained call leaf. Internal RATE evaluations, iterations and stopping '
               'rules remain unknown. This does not eliminate standalone PV/FV programs.')
    with Path(output).open('x')as file:json.dump(result,file,indent=2);file.write('\n')
    return result


def standalone_graph(function,args,reciprocal):
    """Independently construct the advertised worksheet-boundary hypothesis."""
    require(type(reciprocal)is bool,'reciprocal flag is not logical')
    require(len(args)==(2 if function=='POWER'else 1),'wrong standalone arity')
    nodes=[dict(op='input',a=i)for i in range(len(args))]
    one=len(nodes);nodes.append(dict(op='const',bits='3ff0000000000000'))
    def emit(kind,left=0,right=0,spill=True):
        n=dict(op=kind,a=left,b=right)
        if spill:n['site']=sum('site'in n for n in nodes)
        nodes.append(n);return len(nodes)-1
    if function=='LN':
        ln2=emit('fldln2',spill=False);result=emit('fyl2x',ln2,0)
    else:
        if function=='EXP':
            log2e=emit('fldl2e',spill=False);argument=emit('mul',0,log2e)
        else:
            require(function=='POWER','unknown standalone function')
            argument=emit('fyl2x',1,0)
        negative=bool(int(args[0 if function=='EXP'else 1],16)>>63)
        if reciprocal and negative:argument=emit('neg',argument,spill=False)
        integer=emit('frndint',argument,spill=False)
        fraction=emit('sub',argument,integer)
        reduced=emit('f2xm1',fraction)
        mantissa=emit('add',one,reduced)
        result=emit('fscale',mantissa,integer)
        if reciprocal and negative:result=emit('div',one,result)
    return dict(nodes=nodes,outputs=[result,result],sites=sum('site'in n for n in nodes))


def check(packet,trace,gaps,library,output):
    packet=Path(packet);root=Path(__file__).resolve().parents[1]
    manifest=read(packet/'manifest.json');batches=read(packet/'batches.json')
    require(manifest['new_excel_answers']==0,'a prepared packet cannot contain new Excel answers')
    files={str(p.relative_to(packet)):sha(p)for p in packet.rglob('*')if p.is_file()and p!=packet/'manifest.json'}
    require(files==manifest['files_sha256'],'packet file list or checksum differs')
    frozen={str(p):sha(p)for p in (Path(__file__),root/'checker/x87.py',root/'checker/exact.py',Path(trace),Path(gaps),Path(library),packet/'manifest.json')}
    counts=Counter();listing=[]
    for batch in batches:
        doc=read(packet/batch['path']);calls=doc['witnesses'];fn=batch['function']
        require(doc['function']==fn and len(calls)==batch['calls'],'batch function/count differs')
        seen=set();ids=set()
        for row in calls:
            require(set(row)<=set(('id','args','result_index','input_transport')),'an answer or unexpected field is in a request')
            key=(json.dumps(row['args'],sort_keys=True),json.dumps(row.get('result_index')),row.get('input_transport','cell_value2_bulk'))
            require(key not in seen and row['id']not in ids,'duplicate request in a batch')
            require(row.get('input_transport','cell_value2_bulk')in ('cell_value2_bulk','formula_literals'),'unknown transport')
            seen.add(key);ids.add(row['id']);counts[batch['origin']]+=1
            listing.append(dict(origin=batch['origin'],function=fn,**row))
    for name,n in (('round11',64),('structural-r11-v2',52),('round12-pmt-groups',64)):
        require(counts[name]==n,'unanswered packet obligation omitted')
        source=root/'packets'/name;copied=packet/'prior'/name
        src={str(p.relative_to(source)):sha(p)for p in source.rglob('*')if p.is_file()}
        dst={str(p.relative_to(copied)):sha(p)for p in copied.rglob('*')if p.is_file()}
        require(src==dst,'immutable original packet changed')
    with (packet/'READING-LIST.tsv').open(newline='')as file:
        tsv=list(csv.DictReader(file,delimiter='\t'))
    require(len(tsv)==len(listing)==manifest['calls'],'reading list does not cover every batch')
    for i,(a,b)in enumerate(zip(tsv,listing),1):
        require(a['ordinal']==str(i)and all(a[k]==b[k]for k in ('origin','function','id')),'reading list order/identity differs')
        require(a['transport']==b.get('input_transport','cell_value2_bulk')and json.loads(a['arguments_json'])==b['args'],'reading list changed arguments/transport')
    source={r['id']:r for r in records(trace)}
    gap_ids=[r['id']for r in records(gaps)]
    byid={r['id']:r for r in listing if r['origin']=='r13'}
    require(len(byid)==counts['r13'],'duplicate R13 ID')
    selection=read(packet/'selection-r13.json')
    pmt=[r for r in selection if r['function']=='PMT']
    require([r['source_id']for r in pmt]==gap_ids and len(pmt)==24,'original gap rechecks incomplete')
    exact=Exact(Instructions(library));checked=Counter();covered=set()
    for s in selection:
        request=byid[s['id']];original=source[s['source_id']];covered.add(s['id'])
        require(request['function']==s['function'],'selection function differs')
        args=request['args'];d=s['derivation'];r=binary64(original['args'][0])
        base=bits64(B.round(E.round(1+r)));n=int(original['args'][1],16)
        if s['function']=='PMT':require(args==original['args'],'changed original gap inputs')
        elif s['function']=='LN':require(args==[f'0x{base:016x}'],'wrong stored base LN input')
        elif s['function']=='POWER':
            require(args[0]==f'0x{base:016x}'and int(args[1],16)in (n,n^(1<<63)),'wrong POWER base/exponent')
        elif s['function']=='EXP':
            require(d['cw']in (0x123f,0x133f)and int(d['periods'],16)in (n,n^(1<<63)),'wrong exponent derivation')
            cw=d['cw'];zero=encode(0)
            rom=exact.instructions.run('fldln2',zero,zero,cw)
            log=exact.instructions.run('fyl2xp1',rom,encode(r),cw)
            t=exact.operation('mul',encode(binary64(d['periods'])),log,cw,False)
            require(t.hex()==d['resident_t'],'native exponent transcript differs')
            require(args==[f'0x{bits64(B.round(decode(t)),bool(t[9]&128)):016x}'],'EXP input does not match fixed native derivation')
        else:raise ValueError('unknown R13 audit function')
        checked['input_derivations']+=1
    require(covered==set(byid),'an R13 request has no selection/derivation')
    seen=set()
    for prediction in records(packet/'predictions/native-elementary.jsonl'):
        request=byid[prediction['id']];graph=prediction['graph'];cw=prediction['cw']
        require(request['function']==prediction['function']and prediction['function']in ('LN','POWER','EXP'),'unexpected standalone prediction')
        require(graph==standalone_graph(prediction['function'],request['args'],prediction['reciprocal']),'standalone function graph differs')
        key=(prediction['id'],prediction['reciprocal'],cw)
        require(key not in seen,'duplicate elementary hypothesis');seen.add(key)
        values=prediction['predictions_by_store_mask']
        require(len(values)==1<<graph['sites'],'missing standalone store mask')
        for mask,proposed in enumerate(values):
            raw=exact.graph(graph,[int(x,16)for x in request['args']],cw,mask)[:10]
            require(proposed==f'0x{bits64(B.round(decode(raw)),bool(raw[9]&128)):016x}','native standalone prediction differs')
            checked['standalone_predictions']+=1
    expected={(r['id'],reciprocal,cw)for r in byid.values()if r['function']!='PMT'
              for reciprocal in ((False,)if r['function']=='LN'else(False,True))for cw in (0x123f,0x133f)}
    require(seen==expected,'incomplete standalone hypothesis coverage')
    require(counts==manifest['origin_calls'],'manifest origin counts differ')
    require(frozen=={p:sha(p)for p in frozen},'sources/evidence changed')
    result=dict(complete=True,calls=len(listing),origin_calls=dict(counts),checked=dict(checked),
        native_microcode_replays=dict(exact.instructions.counts),inputs_sha256=frozen,
        packet_file_sha256=files,new_excel_answers=0,
        limits='Prepared questions and hypothesis predictions only. Original180 obligations preserved exactly, fixed derivations checked. No answer, semantic pass, surviving PMT model or CPU portability inferred.')
    Path(output).write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == 'rate-newton':
        print(json.dumps(check_rate_newton(*sys.argv[2:])))
        raise SystemExit(0)
    if len(sys.argv) > 1 and sys.argv[1] == 'annuity-questions':
        print(json.dumps(check_annuity_questions(*sys.argv[2:])))
        raise SystemExit(0)
    if len(sys.argv) > 1 and sys.argv[1] == 'score-readings':
        print(json.dumps(score_readings(*sys.argv[2:])))
        raise SystemExit(0)
    if len(sys.argv) > 1 and sys.argv[1] == 'readings':
        r = check_readings(*sys.argv[2:])
        print(json.dumps({k:v for k,v in r.items() if k not in
                          ('batches', 'result_files_sha256', 'packet_files_sha256')}))
        raise SystemExit(0)
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('packet','trace','gaps','library','output'):parser.add_argument(key)
    a=parser.parse_args();r=check(a.packet,a.trace,a.gaps,a.library,a.output)
    print(json.dumps({k:v for k,v in r.items()if k not in ('packet_file_sha256','inputs_sha256')}))
