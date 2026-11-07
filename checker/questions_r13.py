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


def delivery_files(root):
    """Check an immutable owner delivery without importing its claimed verdicts."""
    root = Path(root).resolve()
    manifest = read(root/'delivery-sha256.json')
    total = 0
    for relative, expected in manifest.items():
        path = (root/relative).resolve()
        require(path.is_relative_to(root), 'delivery path escapes root')
        require(path.is_file() and sha(path) == expected, 'delivery file changed: '+relative)
        total += path.stat().st_size
    return dict(files=len(manifest), bytes=total,
                manifest_sha256=sha(root/'delivery-sha256.json'))


def check_workbook_delivery(root, output):
    """Preserve and verify a corpus-fuzzer workbook delivery (R24 and later).

    Checks original question identities, canonical observations, raw worker
    values, actual workbook bytes and operand echoes. Reported prediction
    counts in a delivery are not trusted or copied into observation counts.
    """
    import zipfile
    import xml.etree.ElementTree as ET
    root=Path(root)
    files=delivery_files(root)
    packet=root/'packet'
    question_paths=sorted(packet.rglob('workbook-questions.json'))
    require(question_paths,'missing workbook questions')
    answers=records(root/'answers.jsonl')
    summaries=[]
    seen=set()
    for question_path in question_paths:
        document=read(question_path)
        sheet=document['sheet']
        require(sheet not in seen,'duplicate question worksheet')
        seen.add(sheet)
        selected=[a for a in answers if a.get('observation',{}).get('sheet')==sheet
                  or a.get('raw_cell',{}).get('sheet')==sheet
                  or a.get('condition','').endswith('-'+sheet)]
        summaries.append(check_workbook_readings(root,question_path,selected))
    require(sum(s['formula_cells']for s in summaries)==len(answers),'uncounted workbook answers')
    if len(summaries)==1:
        result=dict(summaries[0],files=files)
    else:
        keys=('roles','statuses','observed_calls_by_function','intended_operand_echoes')
        totals={k:Counter()for k in keys}
        for summary in summaries:
            for k in keys:totals[k].update(summary[k])
        result=dict(complete=True,files=files,formula_cells=len(answers),**totals,
                    workbooks=summaries,
                    observed_formula_spelling_changes=sum(s['observed_formula_spelling_changes']for s in summaries),
                    limitations=summaries[0]['limitations'])
    with Path(output).open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    return result


def check_workbook_readings(root, question_path, answers):
    """Check one sheet in either a single-book or multi-book delivery."""
    import re
    import zipfile
    import xml.etree.ElementTree as ET
    document=read(question_path)
    sheet=document['sheet'];questions=document['questions']
    workbooks=list((question_path.parent/'workbooks').glob('*/oracle.xlsx'))
    require(len(workbooks)==1,'expected one original question workbook')
    workbook=workbooks[0]
    captures=[p for p in root.glob('*.capture.json') if read(p).get('source')=='sha256:'+sha(workbook)]
    require(len(captures)==1,'expected one canonical workbook capture')
    capture=read(captures[0])
    require(capture['source']=='sha256:'+sha(workbook),'capture is for different workbook bytes')
    require(capture['acquisition']['kind']=='live_capture','not a live capture')
    observations={(r['sheet'],r['row'],r['column']):r for r in capture['observations']}
    require(len(observations)==len(capture['observations']),'duplicate capture observation')
    formulas={(r['sheet'],r['row'],r['column']):r['formula']for r in capture['formulas']}
    require(len(answers)==len(questions),'missing or additional workbook answers')
    ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(workbook)as archive:
        relationships={r.attrib['Id']:r.attrib['Target']for r in ET.fromstring(archive.read('xl/_rels/workbook.xml.rels'))}
        sheets=ET.fromstring(archive.read('xl/workbook.xml')).find('m:sheets',ns)
        selected=[s for s in sheets if s.attrib['name']==sheet]
        require(len(selected)==1,'missing question worksheet')
        target=relationships[selected[0].attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']]
        part=target.lstrip('/')if target.startswith('/')else'xl/'+target
        require('..' not in Path(part).parts,'worksheet relationship escapes package')
        cells={c.attrib['r']:c.find('m:f',ns).text for c in ET.fromstring(archive.read(part)).findall('.//m:c',ns)
               if c.find('m:f',ns)is not None}
    roles,statuses,echoes,functions=Counter(),Counter(),Counter(),Counter()
    formula_changes=[]
    seen=set()
    errors={2000:'null',2007:'div0',2015:'value',2023:'ref',2029:'name',2036:'num',2042:'na'}
    for question,answer in zip(questions,answers):
        require(question==answer['question'],'changed workbook question identity')
        key=(sheet,question['row'],question['column'])
        require(key not in seen,'duplicate workbook question')
        seen.add(key);roles[question['role']]+=1;statuses[answer['status']]+=1
        require(question['column']==4,'unsupported question column')
        require(cells['D'+str(question['row'])]==question['formula'],'workbook formula differs from question')
        if answer['status']!='observed':
            require(key not in observations,'unobserved question has a capture value')
            continue
        require(answer['observation']==observations[key],'answer differs from canonical capture')
        worker_formula=answer['raw_cell']['formula']
        worker_formula=''.join(chr(int(worker_formula[i:i+4],16))for i in range(0,len(worker_formula),4))
        require(formulas[key]==worker_formula,'canonical formula differs from raw worker')
        if formulas[key]!='='+question['formula']:
            formula_changes.append(dict(id=question['id'],requested=question['formula'],observed=formulas[key]))
        raw=answer['raw_cell']['value']
        if 'error_code'in raw:
            require(raw['error_code']in errors,'unknown raw error')
            raw={'error':errors[raw['error_code']]}
        require(raw==answer['observation']['excel_raw'],'raw worker differs from capture')
        if question['role']=='operand_echo'and question.get('intended_bits'):
            echoes['equal'if raw=={'number':question['intended_bits'].removeprefix('0x')}else'different']+=1
        if question['role']=='question':
            names=re.findall(r'([A-Za-z_][A-Za-z_0-9.]*)\s*\(',question['formula'])
            name=question.get('function')or next((n.removeprefix('_xlfn.')for n in names
                         if n.upper()not in ('INDEX','TRANSPOSE')), 'arithmetic_control')
            functions[name]+=1
    require(set(observations)<=seen,'uncounted canonical observations')
    result=dict(complete=True,formula_cells=len(questions),roles=roles,statuses=statuses,
        observed_calls_by_function=functions,intended_operand_echoes=echoes,
        observed_formula_spelling_changes=len(formula_changes),first_formula_changes=formula_changes[:8],
        workbook=str(workbook),workbook_sha256=sha(workbook),capture=str(captures[0]),capture_sha256=sha(captures[0]),
        capture_failures=capture.get('session',{}).get('failures',[]),
        limitations='Delivery observations only. No prediction is an Excel answer. Array/decimal neighbor intent requires its effective operand echoes. Product scores require a new application replay.')
    return result


def check_r25_deliveries(r23, eng12, checkout, output):
    """Preserve entry observations and distinguish reduced from original answers.

    No witness, capture or expected value is authored by this intake. The five
    refuted reductions are checked against their actual constant dependencies.
    """
    import collections
    import posixpath
    import zipfile
    import xml.etree.ElementTree as ET
    roots = dict(r23=Path(r23), eng12=Path(eng12))
    files = {name:delivery_files(path) for name,path in roots.items()}
    entry = read(roots['r23']/'r23-entry-v2.capture.json')
    observations = {(r['sheet'],r['row'],r['column']):r for r in entry['observations']}
    require(len(observations) == len(entry['observations']), 'duplicate entry observation')
    questions = read(roots['r23']/'packet/workbook-questions.json')['questions']
    answers = records(roots['r23']/'answers.jsonl')
    require(len(questions) == len(answers) == 968, 'R23 question count changed')
    roles, statuses, echoes = collections.Counter(), collections.Counter(), collections.Counter()
    seen = set()
    for question, answer in zip(questions, answers):
        require(question == answer['question'], 'R23 question identity changed')
        key = ('dcztvv_r23',question['row'],question['column'])
        require(key not in seen, 'duplicate R23 question')
        seen.add(key)
        roles[question['role']] += 1
        statuses[answer['status']] += 1
        if answer['status'] == 'observed':
            require(answer['observation'] == observations[key], 'R23 answer differs from capture')
            raw_value = answer['raw_cell']['value']
            if 'error_code' in raw_value:
                errors = {2000:'null',2007:'div0',2015:'value',2023:'ref',2029:'name',2036:'num',2042:'na'}
                require(raw_value['error_code'] in errors, 'unknown raw Excel error')
                raw_value = {'error':errors[raw_value['error_code']]}
            require(raw_value == answer['observation']['excel_raw'],
                    'R23 canonical answer differs from raw worker')
            intended = question.get('intended_bits')
            if question['role'] == 'operand_echo' and intended:
                actual = answer['observation']['excel_raw']
                echoes['equal' if actual == {'number':intended.removeprefix('0x')} else 'different'] += 1
        else:
            require(answer['status'] == 'entry_refused' and key not in observations,
                    'unsupported R23 answer status')
    require(set(observations) <= seen, 'uncounted R23 observations')

    manifest = read(roots['eng12']/'packet/corpus-cases.json')
    cases = {r['sha256']:r for r in manifest['cases']}
    confirmations = records(roots['eng12']/'independent-confirmations.jsonl')
    require(len(cases) == len(confirmations) == 1146, 'ENG12 case count changed')
    ns = {'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    def cells(path):
        result = {}
        with zipfile.ZipFile(path) as archive:
            rels = {r.attrib['Id']:r.attrib['Target'] for r in
                    ET.fromstring(archive.read('xl/_rels/workbook.xml.rels'))}
            for sheet in ET.fromstring(archive.read('xl/workbook.xml')).find('m:sheets',ns):
                target = rels[sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']]
                part = posixpath.normpath(target.lstrip('/') if target.startswith('/') else 'xl/'+target)
                require(part.startswith('xl/'), 'worksheet relationship escapes package folder')
                result[sheet.attrib['name']] = {c.attrib['r']:c for c in
                    ET.fromstring(archive.read(part)).findall('.//m:c',ns)}
        return result
    original = cells(Path(checkout)/'crates/calc-func/tests/fixtures/fn2_questions/oracle.xlsx')
    statuses12, refuted, seen12 = collections.Counter(), [], set()
    for confirmation in confirmations:
        identity = confirmation['sha256']
        require(identity not in seen12 and identity in cases, 'duplicate/unknown ENG12 identity')
        seen12.add(identity)
        case = cases[identity]
        package = roots['eng12']/'packet/corpus-cases'/case['case']
        require(sha(package) == identity, 'ENG12 reduced package changed')
        require(confirmation['inherited_reference'] == case['source_reference_observation'],
                'ENG12 inherited answer changed')
        statuses12[confirmation['status']] += 1
        if confirmation.get('capture'):
            capture = roots['eng12']/confirmation['capture']
            require(sha(capture) == confirmation['capture_sha256'], 'ENG12 capture changed')
            record = read(capture)
            target = (confirmation['sheet'],confirmation['row'],confirmation['column'])
            rows = [r for r in record['observations'] if (r['sheet'],r['row'],r['column']) == target]
            if confirmation.get('excel_raw') is not None:
                require(len(rows) == 1 and rows[0]['excel_raw'] == confirmation['excel_raw'],
                        'ENG12 answer differs from canonical capture')
        if confirmation['status'] != 'refuted_frozen_mismatch':
            continue
        reduced = cells(package)
        observer = reduced[confirmation['sheet']]['D1'].find('m:f',ns).text
        producer_sheet = observer.split('!')[0].strip("'")
        producer = reduced[producer_sheet]['J20']
        require(producer.find('m:f',ns) is None, 'refuted dependency was not frozen')
        refuted.append(dict(case=identity, observer_sheet=confirmation['sheet'],
            observer_formula=observer, producer_sheet=producer_sheet,
            original_producer_formula=original[producer_sheet]['J20'].find('m:f',ns).text,
            reduced_constant_type=producer.attrib.get('t','n'),
            reduced_constant=producer.find('m:v',ns).text,
            inherited_reference=confirmation['inherited_reference'],
            measured_reduced_answer=confirmation['excel_raw'],
            cause='The reduction froze the differing producer to the DNA value. '
                  'Excel confirms that constant in the reduced file; the original call remains open.'))
    require(len(refuted) == 5, 'ENG12 refutation count changed')
    result = dict(complete=True, files=files,
        r23=dict(roles=roles,statuses=statuses,intended_operand_echoes=echoes,
                 capture_sha256=sha(roots['r23']/'r23-entry-v2.capture.json')),
        eng12=dict(statuses=statuses12, refuted=refuted),
        limits='Measurements of the authored inputs only. Original answers remain unchanged; '
               'refutations do not repair the original corpus. Capture failures and missing answers fail. '
               'R23 effective operands require the associated echoes; intended bits are not observations.')
    with Path(output).open('x') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    return result


def frozen_prediction_path(results,name,expected_hash):
    """Locate one original worker model, including retained nested retries."""
    require(Path(name).name==name and name.endswith('.json') and
            not any(c in name for c in '*?[]'), 'invalid frozen worker name')
    found=list(Path(results).rglob('worker/frozen/'+name))
    require(len(found)==1, 'missing or ambiguous original frozen worker')
    require(sha(found[0])==expected_hash, 'frozen worker prediction changed')
    return found[0]


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
        require(doc['capture_provenance']['calculation_settings'].get('Date1904') is False,
                'scorer supports the recorded 1900 date system only')
        prediction = output/f'prediction-{i:02d}.jsonl'
        subprocess.run([str(executable), 'invoke-inputs', str(path), str(prediction), '1900'], check=True,
                       timeout=300)  # ubs:ignore — explicit compiled executable and argv, no shell; binary hash bound below
        actuals = records(prediction)
        require(len(actuals) == len(doc['witnesses']), 'missing calculated result')
        for original, actual in zip(doc['witnesses'], actuals):
            require(original['id'] == actual['id'] and original['args'] == actual['args'],
                    'calculated identity differs')
            env = doc['capture_provenance']['environment']
            rows.append(dict(id=original['id'], function=doc['function'], args=original['args'],
                input_transport='cell_value2_bulk',
                result_index=original.get('result_index'), expected=original['expected_bits'],
                actual=actual['calc_func_prediction'],
                reading=dict(bitness=env['excel_bitness'], cpu=env['cpu_id'],
                             excel_build=env['excel_build'], platform=env['platform']),
                source=f"readings/excel-20527/r{intake['packet_round']}-"+path.name,
                capture_provenance=doc['capture_provenance'], source_sha256=batch['reading_sha256']))
            status = 'matched' if rows[-1]['expected'] == rows[-1]['actual'] else 'failed'
            counts[status] += 1
        require(sha(path) == batch['reading_sha256'], 'reading changed during scoring')
    require(sha(executable) == binary_hash and len(rows) == intake['calls'], 'scoring execution changed')
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
    received_index = read(results/'packet-sha256.json')
    # The R15 runner exported the manifest's payload map; the other runners
    # exported the whole directory map. Whole-tree byte equality below still
    # verifies manifest.json itself against the original research packet.
    if manifest.get('calls') == 236:
        require({k: v for k, v in files.items() if k != 'manifest.json'} == received_index,
                'received R15 payload manifest differs')
    else:
        require(files == received_index, 'received packet manifest differs')
    require(files == {p.relative_to(received).as_posix(): sha(p)
                      for p in received.rglob('*') if p.is_file()}, 'received packet bytes differ')
    if manifest.get('format') == 'dnacalc-research-questions-r14-variance/1':
        packet_round, expected_calls, expected_batches = 14, 44, 6
        require(manifest['requests'] == expected_calls and manifest['counts'] ==
                {'VAR.S': 10, 'VAR.P': 6, 'STDEV.S': 6, 'STDEV.P': 6, 'DEVSQ': 10, 'AVERAGE': 6},
                'R14 packet scope differs')
        batch_list = [dict(function=f, path=f'requests-{f}.json', calls=count)
                      for f, count in manifest['counts'].items()]
        file_hashes = manifest['files']
    elif manifest.get('format') == 'dnacalc-research-questions-r16/1':
        packet_round, expected_calls, expected_batches = 16, 250, 13
        require(manifest['scalar_calls'] == expected_calls, 'R16 packet scope differs')
        batch_list = read(packet/'batches.json')
        file_hashes = manifest['files_sha256']
    elif manifest.get('format') == 'dnacalc-research-questions-r18/1':
        packet_round, expected_calls, expected_batches = 18, 44, 2
        require(manifest['scalar_calls'] == expected_calls and
                manifest['calls_by_function'] == {'YIELD': 29, 'CONCAT': 15},
                'R18 packet scope differs')
        batch_list = read(packet/'batches.json')
        file_hashes = manifest['files_sha256']
    elif manifest.get('format') == 'dnacalc-r20-discrete-questions/1':
        packet_round,expected_calls,expected_batches=20,44,2
        require(manifest['calls']==expected_calls and manifest['functions']==
                {'BINOM.DIST':36,'HYPGEOM.DIST':8}, 'R20 packet scope differs')
        batch_list=read(packet/'batches.json');file_hashes=manifest['files_sha256']
    elif manifest.get('format') == 'dnacalc-research-questions-r22/1':
        packet_round,expected_calls,expected_batches=22,49,1
        require(manifest['scalar_calls']==expected_calls, 'R22 scalar packet scope differs')
        batch_list=read(packet/'batches.json');file_hashes=manifest['files_sha256']
        require(batch_list==[dict(function='BESSELI',path='requests-BESSELI.json',calls=49)],
                'R22 scalar batch scope differs')
    elif manifest.get('calls') == 236 and manifest.get('calls_by_function') == {
            'PMT': 24, 'IPMT': 26, 'PPMT': 26, 'CUMPRINC': 38, 'CUMIPMT': 38,
            'FV': 30, 'PV': 30, 'NPER': 8, 'RATE': 16}:
        packet_round, expected_calls, expected_batches = 15, 236, 9
        batch_list = read(packet/'batches.json')
        file_hashes = manifest['files_sha256']
    else:
        packet_round, expected_calls, expected_batches = 13, 234, 13
        batch_list = read(packet/'batches.json')
        file_hashes = manifest['files_sha256']
    for name, expected in file_hashes.items():
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
    packet_predictions = Counter()
    for batch in batch_list:
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
        model_path = frozen_prediction_path(results,model_name,model['frozen_sha256'])
        frozen_worker = read(model_path)
        require(frozen_worker['function'] == function, 'frozen worker prediction function differs')
        require(frozen_worker['inputs_sha256'] == request_hash and frozen_worker['input_document'] == request,
                'frozen worker input document differs')
        worker = model_path.parent.parent/'read'/(model_name[:-5]+'-1')
        raw, submitted = read(worker/'result.json'), read(worker/'request.json')
        require(raw['fatal'] is None and not raw['failures'], 'worker failed')
        require(submitted['function'] == function, 'submitted function differs')
        require(raw['automation_security_before'] == raw['automation_security_after'],
                'worker security settings changed')
        for key in ('automation_security_before', 'automation_security_after'):
            require(provenance[key] == raw[key], 'security provenance differs')
        settings = raw['calculation_settings']
        require(provenance['calculation_settings'] == settings, 'settings provenance differs')
        # Environment is provenance, never a whitelist of qualifying captures.
        # Check the reader's transcription against every field the worker exposes.
        environment = provenance['environment']
        require(all(key in raw for key in ('excel_version', 'excel_build',
                'excel_file_version', 'excel_machine', 'os_version', 'cpu')),
                'missing worker environment provenance')
        for canonical_key, worker_key in (
                ('excel_file_version', 'excel_file_version'), ('excel_version', 'excel_version'),
                ('excel_build', 'excel_build'), ('cpu_id', 'cpu'),
                ('excel_operating_system', 'os_version'), ('platform', 'platform'),
                ('excel_bitness', 'excel_bitness')):
            if worker_key in raw:
                require(environment[canonical_key] == raw[worker_key], 'environment provenance differs')
        bitness = {'x86-64': '64-bit', 'x86': '32-bit', 'arm64': '64-bit'}.get(raw.get('excel_machine'))
        if bitness is not None:
            require(environment['excel_bitness'] == bitness, 'environment bitness differs')
        count = batch['calls']
        frozen_packet_predictions = None
        if packet_round in (15, 16, 18, 20, 22):
            frozen_packet_predictions = records(packet/f'application-predictions-{function}.jsonl')
            require(len(frozen_packet_predictions) == count, 'frozen packet prediction coverage differs')
        require(all(len(rows) == count for rows in (request['witnesses'], doc['witnesses'],
                doc['input_audits'], submitted['rows'], raw['rows'], frozen_worker['predictions'])), 'batch coverage differs')
        # The reader groups mixed-arity requests before submitting them. Keep
        # its original runtime order, and join to the frozen request by unique ID.
        original_ids = [r['id'] for r in request['witnesses']]
        require(len(set(original_ids)) == count, 'duplicate original request ID')
        runtime = [doc['witnesses'], doc['input_audits'], submitted['rows'], raw['rows']]
        runtime_ids = [[r['id'] for r in rows] for rows in runtime]
        require(all(ids == runtime_ids[0] for ids in runtime_ids), 'runtime request/answer order differs')
        def indexed(rows):
            ids = [r['id'] for r in rows]
            require(len(set(ids)) == count and set(ids) == set(original_ids), 'duplicate or missing batch identity')
            return {r['id']:r for r in rows}
        canonical_by_id, audit_by_id, call_by_id, answer_by_id = map(indexed, runtime)
        frozen_by_id = indexed(frozen_worker['predictions'])
        require([r['id'] for r in frozen_worker['predictions']] == original_ids, 'frozen prediction request order differs')
        packet_by_id = indexed(frozen_packet_predictions) if frozen_packet_predictions is not None else None
        for asked in request['witnesses']:
            qid = asked['id']
            canonical, audit, call, answer = (rows[qid] for rows in
                (canonical_by_id, audit_by_id, call_by_id, answer_by_id))
            require(qid not in seen, 'duplicate packet ID'); seen.add(qid)
            require(qid == canonical['id'] == audit['id'] == call['id'] == answer['id'],
                    'request/answer order differs')
            require(asked['args'] == canonical['args'] and
                    asked.get('result_index') == canonical.get('result_index'), 'typed identity differs')
            frozen_row = frozen_by_id[qid]
            require(frozen_row['id'] == qid and frozen_row['args'] == asked['args']
                    and frozen_row.get('result_index') == asked.get('result_index')
                    and frozen_row['calc_func_prediction'] == canonical['calc_func_prediction'],
                    'canonical/frozen worker prediction differs')
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
            if frozen_packet_predictions is not None:
                prediction = packet_by_id[qid]
                require(prediction['id'] == qid and prediction['args'] == asked['args']
                        and prediction.get('not_excel_evidence') is True,
                        'frozen packet prediction identity differs')
                packet_predictions['total'] += 1
                packet_predictions['matched' if prediction['calc_func_prediction'] == expected else 'failed'] += 1
        total += count
        batches.append(dict(**batch, reading=path.relative_to(results).as_posix(),
                            reading_sha256=sha(path), request_sha256=sha(request_path),
                            worker_root=model_path.parent.parent.relative_to(results).as_posix(),
                            original_request_ids=original_ids, runtime_ids=runtime_ids[0],
                            reader_reordered=original_ids != runtime_ids[0]))
    require(not exports and total == expected_calls and len(batches) == expected_batches, 'missing/excess batch')
    require(frozen == {p.relative_to(results).as_posix(): sha(p)
                      for p in results.rglob('*') if p.is_file()}, 'result tree changed')
    result = dict(complete=True, packet_round=packet_round, calls=total, batches=batches, result_files_sha256=frozen,
                  packet_files_sha256=files, frozen_windows_application=dict(predictions),
                  frozen_packet_application=dict(packet_predictions),
                  limits='Existing reader exports only; Value2 transport audited. Effective callee '
                         'operands remain unobservable. Native CPU measurements and predictions '
                         'are never Excel answers. Repeated captures are observations, not distinct inputs. '
                         'Mixed-arity runtime grouping is retained; unique-ID bijections join it to the frozen original request, with runtime ordering checked independently.')
    with Path(output).open('x') as f:
        json.dump(result, f, indent=2); f.write('\n')
    return result


def check_round16_questions(packet, readings, workbook, output):
    """Independently re-derive every requested tuple and catalogue formula."""
    import itertools
    import math
    import struct
    import zipfile
    import xml.etree.ElementTree as ET
    packet,readings,workbook=map(Path,(packet,readings,workbook))
    manifest=read(packet/'manifest.json')
    checker_hash=sha(Path(__file__))
    input_hashes={str(p):sha(p)for p in (workbook,*sorted(readings.glob('*.json')))}
    require(manifest['format']=='dnacalc-research-questions-r16/1','wrong packet format')
    files={p.relative_to(packet).as_posix():sha(p)for p in packet.rglob('*')if p.is_file()and p.name!='manifest.json'}
    require(files==manifest['files_sha256'],'changed/missing/extra packet file')
    for source,expected_hash in manifest['sources_sha256'].items():
        require(sha(readings/Path(source).name)==expected_hash,'changed R14 source reading')
    bits=lambda x:'0x'+struct.pack('>d',float(x)).hex()
    number=lambda x:struct.unpack('>d',bytes.fromhex(x[2:]))[0]
    lit=lambda x:format(float(x),'.17g')
    r14=read(readings/'VAR.S.json')['witnesses']
    require(len(r14)==10,'missing R14 anchors')
    triples=[list(map(number,r['args']))for r in r14]
    expected={}
    def ask(f,values):expected.setdefault(f,set()).add(tuple(map(bits,values)))
    names=('VAR.S','VAR.P','STDEV.S','STDEV.P','DEVSQ','AVERAGE')
    for t,f in itertools.product(triples,names):ask(f,t)
    limit=float(2**27)
    grid=[limit+v for v in(-2,-1,-.75,-.5,-.25,0,.25,.5,.75,1,2)]+[math.nextafter(limit,0),math.nextafter(limit,math.inf)]
    for x,n in itertools.product(grid,(0,1,2,3,4)):ask('BESSELJ',(x,n))
    for x,n in itertools.product((limit-1,limit-.5,limit,limit+.5),(0,1,2)):ask('BESSELJ',(-x,n))
    for x,n in itertools.product((limit-1,limit-.5,math.nextafter(limit,0),limit,math.nextafter(limit,math.inf),limit+1),(0,1,2)):ask('BESSELY',(x,n))
    for x,n in itertools.product((709,709.5,math.nextafter(709.7,0),709.7,math.nextafter(709.7,math.inf),709.78,710,713,714),(0,1,2,3,10,100)):ask('BESSELI',(x,n))
    for x,n in itertools.product((709.7,math.nextafter(709.7,math.inf),710),(0,1,2)):ask('BESSELI',(-x,n))
    for x in(709,709.5,709.7,math.nextafter(709.7,math.inf),709.78,math.log(float.fromhex('0x1.fffffffffffffp+1023')),710):ask('EXP',(x,))
    phase_grid=(limit-1,limit-.5,math.nextafter(limit,0),limit,math.nextafter(limit,math.inf),limit+.5)
    for x,f in itertools.product(phase_grid,('SIN','COS')):ask(f,(x,))
    geo=list(itertools.permutations((1,5,3)))+[(15,),(1,15),(3,5),(15,1,1),(1,15,1),(1,1,15),(1,3,5,1)]
    for values in geo:ask('GEOMEAN',values)
    found={};ids=set();counts=Counter()
    for batch in read(packet/'batches.json'):
        f=batch['function'];doc=read(packet/batch['path']);require(doc['function']==f,'wrong request function')
        require(f not in found and f in expected,'duplicate/unknown function batch');found[f]=set()
        predictions=records(packet/f'application-predictions-{f}.jsonl')
        require(len(predictions)==len(doc['witnesses'])==batch['calls'],'missing prediction/request')
        for index,(row,pred)in enumerate(zip(doc['witnesses'],predictions)):
            require(set(row)=={'id','args','input_transport'}and row['input_transport']=='cell_value2_bulk','unexpected request answer/transport field')
            require(row['id']==f"r16-{f.lower().replace('.','-')}-{index:03d}"and row['id']not in ids,'changed request identity')
            ids.add(row['id']);args=tuple(row['args']);require(args not in found[f],'duplicate scalar question');found[f].add(args)
            require(pred['id']==row['id']and pred['args']==row['args']and pred['not_excel_evidence']is True,'prediction identity/evidence status changed')
            counts[f]+=1
    require(found==expected and len(ids)==manifest['scalar_calls']==250 and dict(counts)==manifest['calls_by_function'],'scalar derivation/count differs')
    selected=read(packet/'selection.json');require(len(selected)==250 and {r['id']for r in selected}==ids,'missing selection provenance')
    lookup={r['id']:(f,r['args'])for batch in read(packet/'batches.json')for r in read(packet/batch['path'])['witnesses']for f in[batch['function']]}
    for r in selected:require(lookup[r['id']]==(r['function'],r['requested_args']),'selection identity changed')
    formulas=[]
    def add(text):formulas.append(text);return f'D{len(formulas)}'
    def surfaces(f,values,refs):
        literals=list(map(lit,values));span=f'{refs[0]}:{refs[-1]}'
        for body in(','.join(refs),span,'{'+','.join(literals)+'}','{'+ ';'.join(literals)+'}',f'TRANSPOSE({span})',','.join(literals)):
            add(f'{f}({body})')
    for i,values in enumerate(triples):
        refs=[add(lit(v))for v in values]
        for f in names if i<6 else('VAR.S','DEVSQ'):surfaces(f,values,refs)
    require(len(formulas)==294,'internal variance question census')
    for x in phase_grid:
        base=add(lit(x))
        for phase in(1,3):
            p=add(f'{base}-{phase}*PI()/4')
            for f in('SIN','COS'):add(f'{f}({p})');add(f'{f}({base}-{phase}*PI()/4)')
    require(len(formulas)==360,'internal phase question census')
    used=set()
    def regression(delta,s1=1,s2=1,affine=False):
        rows=[[x*s1,(2*x+int(affine)+(delta if i==3 else 0))*s2]for i,x in enumerate((1,2,5,4))]
        x='{'+ ';'.join(','.join(map(lit,r))for r in rows)+'}'
        for const,f in itertools.product(('FALSE','TRUE'),('LINEST','LOGEST')):
            for r,c in((1,1),(1,2),(4,2)):
                text=f'INDEX({f}({{1;2;4;7}},{x},{const},TRUE),{r},{c})'
                if text not in used:add(text);used.add(text)
    for delta in(0,2**-30,2**-26,2**-24,2**-22,2**-20,1e-6,2**-18,2**-16,2**-12):regression(delta)
    for delta,scales in itertools.product((0,2**-24,2**-20,2**-16),((2**-20,2**-20),(2**20,2**20),(1,2**-20),(1,2**20))):regression(delta,*scales)
    regression(0,affine=True)
    for x,const,f in itertools.product(('{1,3;2,5;5,11;4,9}','{1,0,1;0,1,1;1,1,2;2,1,3}','{1,0;2,0;5,0;4,0}'),('FALSE','TRUE'),('LINEST','LOGEST')):
        add(f'INDEX({f}({{1;2;4;7}},{x},{const},TRUE),4,2)')
    require(len(formulas)==696,'internal rank question census')
    for values in geo:surfaces('GEOMEAN',values,[add(lit(v))for v in values])
    require(len(formulas)==810,'internal full question census')
    questions=read(packet/'catalogue-questions.json');require(questions['sheet']=='dcj9k_r16'and len(questions['questions'])==810,'wrong catalogue question inventory')
    predictions=records(packet/'application-predictions-catalogue.jsonl');require(len(predictions)==810,'missing catalogue prediction')
    for row,(q,pred,text)in enumerate(zip(questions['questions'],predictions,formulas),1):
        require(q['id']==f'r16-sheet-{row:04d}'and q['row']==row and q['column']==4 and q['formula']==text,'catalogue derivation changed')
        require(not any(k in q for k in('expected','expected_bits','calc_func_prediction')),'question contains an answer')
        require(pred['id']==q['id']and pred['row']==row and pred['column']==4 and pred['formula']==text and pred['not_excel_evidence']is True,'catalogue prediction identity changed')
    lines=[s for s in(packet/'sheet-dcj9k_r16.txt').read_text().splitlines()if s and not s.startswith('#')]
    require(lines==formulas and sha(packet/'sheet-dcj9k_r16.txt')==manifest['catalogue_section_sha256'],'section mirror differs')
    namespace={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    relns='{http://schemas.openxmlformats.org/officeDocument/2006/relationships}'
    with zipfile.ZipFile(workbook)as z:
        rels={r.attrib['Id']:r.attrib['Target']for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
        selected=[s for s in ET.fromstring(z.read('xl/workbook.xml')).find('m:sheets',namespace)if s.attrib['name']=='dcj9k_r16']
        require(len(selected)==1,'missing registered question sheet')
        xml=ET.fromstring(z.read('xl/'+rels[selected[0].attrib[relns+'id']]))
        actual=[(c.attrib['r'],c.find('m:f',namespace).text)for c in xml.findall('.//m:c',namespace)if c.find('m:f',namespace)is not None]
        require(actual==[(f'D{i}',text)for i,text in enumerate(formulas,1)],'workbook question formula differs')
    require(files=={p.relative_to(packet).as_posix():sha(p)for p in packet.rglob('*')if p.is_file()and p.name!='manifest.json'},'packet changed during check')
    require(sha(Path(__file__))==checker_hash and all(sha(p)==h for p,h in input_hashes.items()),'question checker/input changed')
    result=dict(complete=True,scalar_calls=250,catalogue_formula_cells=810,new_excel_answers=0,
        files_sha256=files,workbook_sha256=sha(workbook),checker_sha256=checker_hash,input_sha256=input_hashes,
        limits='Every scalar input and worksheet formula re-derived from R14 anchors and fixed stated grids. Predictions are calculation outputs only. No measured Excel answer or candidate validation is created. Worksheet literal echoes and capture settings must be inspected when answers return.')
    with Path(output).open('x')as out:json.dump(result,out,indent=2);out.write('\n')
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


def scalar_key(name,args):
    return name,json.dumps(args,sort_keys=True,separators=(',',':'))


def scalar_reply_proof(receipts,cache):
    """One reader for the original scalar request/readback/answer lineage."""
    import gzip
    original={};observations={};raw_counts=Counter();input_counts=Counter()
    for receipt in receipts:
        path=Path(cache)/receipt['file'];require(sha(path)==receipt['sha256'],'raw coupon reply hash changed')
        if str(path)in original:continue
        original[str(path)]=receipt['sha256']
        with gzip.open(path,'rt')as stream:raw=json.load(stream)
        request,reply=raw['request'],raw['reply'];name=request['function']
        require(name==receipt['function']and reply['ok']and reply['version']==1,'invalid coupon worker reply')
        require(reply['automation_security_before']==reply['automation_security_after'],'worker policy changed')
        require(all(k in reply for k in('excel_build','machine','timestamp','request_sha256','session_receipt_sha256')),'missing coupon provenance')
        args=[[{k:v for k,v in a.items()if k!='mode'}for a in row]for row in request['rows']]
        require(reply['inputs']==args and len(reply['rows'])==len(args),'coupon input/readback changed')
        require(len(args)==receipt['calls'],'scalar receipt call count changed')
        require(all(receipt[k]==reply[k]for k in('excel_build','machine','timestamp','request_sha256')),'scalar receipt provenance changed')
        for inputs,answer in zip(args,reply['rows']):
            token=scalar_key(name,inputs)
            require(token not in observations or observations[token]==answer,'repeated coupon answer differs')
            observations[token]=answer;raw_counts[name]+=1
            input_counts[token]+=1
    return original,observations,raw_counts,input_counts


def check_yield_coupon(study, cache, output, day_policy='direct'):
    """Replay original single-coupon calls and the normalized formula independently."""
    import struct
    from fractions import Fraction as Q
    from exact import operation
    study=Path(study);doc=read(study);require(doc['complete'],'incomplete coupon study')
    policies=('direct','published-period-boundaries','nasd-published-period','numerator-published-period')
    require(day_policy in policies,'unknown coupon day policy')
    require(doc['format']=='dnacalc-yield-direct-day-counts/1','direct day controls required')
    sources={str(p):sha(p)for p in(Path(__file__),Path(__file__).with_name('exact.py'))}
    original,observations,raw_counts,input_counts=scalar_reply_proof(doc['parent_receipts']+doc['receipts'],cache)
    yield_inputs=Counter({key:count for key,count in input_counts.items()if key[0]=='YIELD'})
    key=scalar_key
    def encoded(x):return struct.pack('>d',float(x)).hex()
    def number(x):return dict(type='number',bits=encoded(x))
    def observed(name,args):
        token=key(name,args);require(token in observations,'missing independent coupon observation')
        return observations[token]
    def numeric(name,args):
        result=observed(name,args);require(result['type']=='number','nonnumeric day factor')
        return binary64(result['bits'])
    proposed=Counter(key('YIELD',[number(row[k])for k in('settlement','maturity','rate','price','redemption','frequency','basis')])
                     for row in doc['records'])
    require(proposed==yield_inputs,'YIELD observation cohort omitted or duplicated')
    counts=Counter();calls=[];failures=[]
    for index,row in enumerate(doc['records']):
        args=[number(row[k])for k in('settlement','maturity','rate','price','redemption','frequency','basis')]
        require(observed('YIELD',args)==row['answer'],'derived YIELD answer changed')
        factors=[number(row[k])for k in('settlement','maturity','frequency','basis')]
        for name,answer in row['factors'].items():require(observed(name,factors)==answer,'derived coupon factor changed')
        prev,next_=[numeric(name,factors)for name in('COUPPCD','COUPNCD')]
        settlement,maturity=map(Q,(row['settlement'],row['maturity']))
        require(next_==maturity and prev<=settlement<maturity,'not a final coupon interval')
        spans={}
        for label,start,end in(('a',prev,settlement),('e',prev,next_),('dsr',settlement,next_)):
            span=(numeric('DAYS360',[number(start),number(end),dict(type='boolean',value=row['basis']==4)])
                  if row['basis']in(0,4)else end-start)
            require(row['direct_day_counts'][label]==number(span),'derived direct span changed')
            spans[label]=span
        a,e,dsr=[spans[k]for k in('a','e','dsr')]
        if day_policy!='direct':
            a=numeric('COUPDAYBS',factors)
            if row['basis']not in(2,3):e=numeric('COUPDAYS',factors)
            if day_policy=='published-period-boundaries'and settlement==prev:dsr=e
            if day_policy in('nasd-published-period','numerator-published-period')and row['basis']==0:
                import datetime as dt,calendar
                first,last=[dt.date(1899,12,30)+dt.timedelta(days=int(x))for x in(settlement,maturity)]
                ends=lambda d:d.month==2 and d.day==calendar.monthrange(d.year,2)[1]
                if ends(first)and ends(last):dsr+=30-last.day
                elif day_policy=='numerator-published-period'and ends(first)and last.day==31:dsr+=1
        require(e>0 and dsr>=0,'unexpected direct coupon day domain')
        counts['rows']+=1
        if dsr==0:
            counts['zero_day_divisor']+=1
            if day_policy!='direct':
                matched=row['answer']==dict(type='number',bits='0000000000000000')
                for word in('strict','stored_pc64'):counts[word+'_matches']+=matched
                calls.append(dict(id='yield-coupon-'+str(index),function='YIELD',args=['0x'+v['bits']for v in args],
                    numeric_cell_inputs=True,date_system_1900=True,actual='0x0000000000000000',
                    prediction_source='Declared zero remaining-day publication policy; not an Excel reading.'))
                if not matched:failures.append(dict(index=index,answer=row['answer'],reason='zero-day publication differs'))
            else:failures.append(dict(index=index,answer=row['answer'],reason='zero direct remaining-day divisor; no numeric prediction'))
            continue
        for word in('strict','stored_pc64'):
            op=lambda x,y,kind:operation(Q(x),Q(y),kind,word)
            rate,price,red,frequency=[binary64(encoded(row[k]))for k in('rate','price','redemption','frequency')]
            coupon=op(rate,frequency,'div');p=op(price,100,'div');red=op(red,100,'div')
            accrual=op(op(a,e,'div'),coupon,'mul');dirty=op(p,accrual,'add')
            numerator=op(op(red,coupon,'add'),-dirty,'add')
            ratio=op(numerator,dirty,'div');factor=op(op(frequency,e,'mul'),dsr,'div')
            result=encoded(op(ratio,factor,'mul'))
            matched=row['answer'].get('type')=='number'and row['answer']['bits']==result
            counts[word+'_matches']+=matched
            if not matched:failures.append(dict(index=index,word=word,answer=row['answer'],prediction=result))
            if word=='stored_pc64':
                calls.append(dict(id='yield-coupon-'+str(index),function='YIELD',args=['0x'+v['bits']for v in args],
                    numeric_cell_inputs=True,date_system_1900=True,actual='0x'+result,
                    prediction_source='Independent exact normalized coupon formula; not an Excel reading.'))
    require(all(sha(p)==h for p,h in sources.items()),'coupon checker source changed')
    call_path=Path(str(output)+'.calls.jsonl')
    with call_path.open('x')as stream:
        for call in calls:stream.write(json.dumps(call)+'\n')
    result=dict(format='dnacalc-yield-coupon-independent/1',complete=True,counts=counts,failures=failures,
        source_sha256=sources,inputs_sha256={str(study):sha(study),**original},raw_calls=raw_counts,day_policy=day_policy,
        calls_sha256=sha(call_path),
        limits='Independent rational arithmetic and original scalar readbacks. Direct day spans are '
        'measured separately from the YIELD outputs. The optional boundary policy publishes zero '
        'on a zero remaining-day span; the original direct formula leaves those as failures. '
        'Does not prove literal private stores, non-final coupon behavior or Excel callee transport.')
    with Path(output).open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(counts),flush=True);return result


def check_duration_redemption(study,cache,output):
    """Independent raw-readback and exact-prefix proof for redemption division."""
    import itertools,struct,datetime as dt
    from fractions import Fraction as Q
    from exact import operation,two
    study=Path(study);doc=read(study);require(doc['complete'],'incomplete duration study')
    zero=doc['format']=='dnacalc-duration-zero-coupon/1'
    require(zero or doc['format']=='dnacalc-duration-two-coupon-transfer/1','unknown duration study')
    sources={str(p):sha(p)for p in(Path(__file__),Path(__file__).with_name('exact.py'))}
    original,observations,raw_counts,input_counts=scalar_reply_proof(doc['receipts'],cache)
    def encoded(x):return struct.pack('>d',float(x)).hex()
    def number(x):return dict(type='number',bits=encoded(x))
    def observed(name,args):
        token=scalar_key(name,args);require(token in observations,'missing original duration control')
        return observations[token]
    proposed=Counter();states={};prepared=[]
    for row in doc['records']:
        f=row['frequency'];n=row['n']if zero else 2;q=binary64(row['base_strict']if zero else row['base_bits'])
        y=binary64(row['yield_bits']);rate=Q(0)if zero else binary64(row['coupon_rate_bits'])
        require(f in(1,2,4)and n>=2,'duration frequency/count outside scope')
        require(row['settlement']==43831,'changed on-coupon settlement')
        months=2020*12+n*(12//f);year,month=divmod(months,12)
        maturity=(dt.date(year,month+1,1)-dt.date(1899,12,30)).days
        require(row['maturity']==maturity,'changed on-coupon maturity')
        require(B.round(1+y/f)==q,'duration base prefix changed')
        if zero:
            den0=num0=Q(0);powers=row['powers'];power=powers['strict'];answer=row['duration']
            proposed[scalar_key('SUM',[number(1),number(y/f)])]+=1
            require(observed('SUM',[number(1),number(y/f)])==row['public_sum'],'base SUM reading changed')
            require(number(B.round(E.round(1+y/f)))==row['public_sum'],'base SUM arithmetic differs')
            for base in('strict','stored_pc64'):
                args=[number(binary64(row['base_'+base])),number(n)]
                require(observed('POWER',args)==powers[base],'discount POWER reading changed')
                if base=='strict'or row['base_stored_pc64']!=row['base_strict']:proposed[scalar_key('POWER',args)]+=1
            states[row['settlement'],maturity,f]=n
        else:
            require(1<q<2 and q.numerator%5==0 and q.denominator<=1<<25,'controlled dyadic base changed')
            square=q*q;scale=two(row['scale_exponent']);coupon=square*scale
            require(B.round(square)==square==binary64(row['square_bits']),'controlled square is not exact')
            require(B.round(100*rate)==100*rate and 100*rate/f==coupon,'controlled coupon prefix is not exact')
            require(y==(q-1)*f,'controlled yield prefix is not exact')
            den0=(q+1)*scale;num0=(q+2)*scale
            require(B.round(den0)==den0==binary64(row['coupon_denominator_prefix']),'controlled coupon PV sum differs')
            require(B.round(num0)==num0==binary64(row['coupon_numerator_prefix']),'controlled weighted sum differs')
            power=row['power'];answer=row['answer'];args=[number(q),number(2)]
            require(observed('POWER',args)==power==number(square),'exact public square differs')
            proposed[scalar_key('POWER',args)]+=1
        require(power['type']=='number','nonnumeric discount factor')
        args=[number(row['settlement']),number(maturity),number(rate),number(y),number(f),number(0)]
        price_args=args[:4]+[number(100)]+args[4:]
        for name,inputs,expected in(('DURATION',args,answer),('MDURATION',args,row['modified_duration']),('PRICE',price_args,row['price'])):
            proposed[scalar_key(name,inputs)]+=1
            require(observed(name,inputs)==expected,'original '+name+' reading changed')
        prepared.append(dict(row=row,n=n,f=f,q=q,den0=den0,num0=num0,discount=binary64(power['bits']),
                             answer=answer,args=args,price_args=price_args))
    for(s,m,f),n in states.items():
        args=[number(s),number(m),number(f),number(0)]
        for name,expected in(('COUPNUM',n),('COUPDAYBS',0)):
            proposed[scalar_key(name,args)]+=1
            require(observed(name,args)==number(expected),'on-coupon public schedule differs')
    require(proposed==input_counts,'duration observation cohort omitted or duplicated')
    programs=[dict(words=list(words),matches=0,first_failure=None)for words in itertools.product(('strict','stored_pc64'),repeat=5)]
    price_programs=[dict(words=list(words),matches=0,first_failure=None)for words in itertools.product(('strict','stored_pc64'),repeat=2)]
    counts=Counter();calls={profile:[]for profile in('before','denominator','redemption','modified')}
    def prediction(case,words):
        rd,rn,dw,nw,fw=words
        red=operation(Q(100),case['discount'],'div',rd)
        weighted=operation(Q(case['n']*100),case['discount'],'div',rn)
        den=operation(case['den0'],red,'add',dw);num=operation(case['num0'],weighted,'add',nw)
        return operation(num,den,'div',fw)/case['f']
    for index,case in enumerate(prepared):
        row=case['row'];answer=case['answer'];counts['rows']+=1
        for program in programs:
            value=encoded(prediction(case,program['words']));good=answer.get('type')=='number'and answer['bits']==value
            program['matches']+=good
            if not good and program['first_failure']is None:program['first_failure']=dict(index=index,answer=answer,prediction=value)
        for program in price_programs:
            rd,dw=program['words'];red=operation(Q(100),case['discount'],'div',rd)
            value=encoded(operation(case['den0'],red,'add',dw));good=row['price'].get('type')=='number'and row['price']['bits']==value
            program['matches']+=good
            if not good and program['first_failure']is None:program['first_failure']=dict(index=index,answer=row['price'],prediction=value)
        for word in('strict','stored_pc64'):
            if answer.get('type')=='number':
                value=encoded(operation(binary64(answer['bits']),case['q'],'div',word))
                counts['mduration_'+word+'_matches']+=row['modified_duration'].get('type')=='number'and row['modified_duration']['bits']==value
        for profile in calls:
            rd='strict'if profile=='before'else'stored_pc64'
            rn='stored_pc64'if profile in('redemption','modified')else'strict'
            duration=prediction(case,(rd,rn,'strict','strict','strict'))
            mdword='stored_pc64'if profile=='modified'else'strict'
            md=operation(B.round(duration),case['q'],'div',mdword)
            price=operation(case['den0'],operation(Q(100),case['discount'],'div',rd),'add','strict')
            for name,args,result in(('DURATION',case['args'],duration),('MDURATION',case['args'],md),('PRICE',case['price_args'],price)):
                calls[profile].append(dict(id=name.lower()+'-'+str(index),function=name,
                    args=['0x'+v['bits']for v in args],numeric_cell_inputs=True,date_system_1900=True,
                    actual='0x'+encoded(result),prediction_source='Independent exact '+profile+' arithmetic; not an Excel reading.'))
    call_hashes={}
    for profile,rows in calls.items():
        path=Path(str(output)+'.'+profile+'.calls.jsonl')
        with path.open('x')as stream:
            for row in rows:stream.write(json.dumps(row)+'\n')
        call_hashes[profile]=sha(path)
    require(all(sha(p)==h for p,h in sources.items()),'duration checker source changed')
    result=dict(format='dnacalc-duration-redemption-independent/1',complete=True,counts=counts,programs=programs,
        complete_programs=[p for p in programs if p['matches']==len(prepared)],price_programs=price_programs,
        complete_price_programs=[p for p in price_programs if p['matches']==len(prepared)],calls_sha256=call_hashes,
        source_sha256=sources,inputs_sha256={str(study):sha(study),**original},raw_calls=raw_counts,
        limits='Original scalar observations and exact on-coupon prefixes only. Redemption numerator and denominator '
        'words are independent profile coordinates. Remaining sum/ratio words are not inferred when they remain '
        'indistinguishable. No general coupon loop or fractional discount helper is identified.')
    with Path(output).open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(counts=counts,complete=result['complete_programs'],price=result['complete_price_programs'])),flush=True)
    return result


if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='duration-redemption':
        check_duration_redemption(*sys.argv[2:])
        raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='yield-coupon':
        check_yield_coupon(*sys.argv[2:])
        raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='workbook-delivery':
        print(json.dumps(check_workbook_delivery(*sys.argv[2:])))
        raise SystemExit(0)
    if len(sys.argv) > 1 and sys.argv[1] == 'round16-questions':
        print(json.dumps(check_round16_questions(*sys.argv[2:])))
        raise SystemExit(0)
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
