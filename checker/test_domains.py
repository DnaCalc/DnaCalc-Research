"""Constructed controls for adaptive questions; none are Excel observations."""
import copy
import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch
import json
import tempfile

SPEC = importlib.util.spec_from_file_location('domain_bisection',Path(__file__).parents[1]/'racer/domain_bisection.py')
domain = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(domain)


class ReadingRetries(unittest.TestCase):
    def test_multi_workbook_delivery_keeps_every_raw_answer_bound(self):
        import hashlib
        import zipfile
        from checker.questions_r13 import check_workbook_delivery
        with tempfile.TemporaryDirectory()as scratch:
            root=Path(scratch);answers=[]
            for sheet in ('first','second'):
                packet=root/'packet'/sheet;book=packet/'workbooks'/sheet/'oracle.xlsx'
                book.parent.mkdir(parents=True)
                question=dict(id=sheet,row=1,column=4,formula='EXP(0)',role='question')
                (packet/'workbook-questions.json').write_text(json.dumps(dict(sheet=sheet,questions=[question])))
                ns='http://schemas.openxmlformats.org/spreadsheetml/2006/main'
                with zipfile.ZipFile(book,'w')as z:
                    z.writestr('xl/workbook.xml',f'<workbook xmlns="{ns}" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="{sheet}" r:id="r1"/></sheets></workbook>')
                    z.writestr('xl/_rels/workbook.xml.rels','<Relationships><Relationship Id="r1" Target="worksheets/sheet1.xml"/></Relationships>')
                    z.writestr('xl/worksheets/sheet1.xml',f'<worksheet xmlns="{ns}"><sheetData><row r="1"><c r="D1"><f>EXP(0)</f></c></row></sheetData></worksheet>')
                observation=dict(sheet=sheet,row=1,column=4,excel_raw={'number':'3ff0000000000000'})
                capture=dict(source='sha256:'+hashlib.sha256(book.read_bytes()).hexdigest(),
                             acquisition={'kind':'live_capture'},observations=[observation],
                             formulas=[dict(sheet=sheet,row=1,column=4,formula='=EXP(0)')])
                (root/(sheet+'.capture.json')).write_text(json.dumps(capture))
                answers.append(dict(question=question,condition='r27-'+sheet,status='observed',
                                    observation=observation,raw_cell=dict(sheet=sheet,
                                    formula=''.join(f'{ord(c):04x}'for c in '=EXP(0)'),
                                    value=observation['excel_raw'].copy())))
            def bind():
                (root/'answers.jsonl').write_text(''.join(json.dumps(a)+'\n'for a in answers))
                (root/'delivery-sha256.json').write_text(json.dumps({str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in root.rglob('*')if p.is_file()and p.name!='delivery-sha256.json'}))
            bind()
            result=check_workbook_delivery(root,root/'valid.json')
            self.assertEqual((result['formula_cells'],len(result['workbooks'])),(2,2))
            self.assertEqual(result['observed_calls_by_function'],{'EXP':2})
            answers[1]['raw_cell']['value']={'number':'4000000000000000'};bind()
            with self.assertRaisesRegex(ValueError,'raw worker differs'):
                check_workbook_delivery(root,root/'invalid.json')
            self.assertFalse((root/'invalid.json').exists())

    def test_delivery_files_reject_changed_bytes_and_escaping_paths(self):
        import hashlib
        from checker.questions_r13 import delivery_files
        with tempfile.TemporaryDirectory() as scratch:
            root=Path(scratch);payload=root/'answer.json';payload.write_bytes(b'constructed answer')
            manifest=root/'delivery-sha256.json'
            manifest.write_text(json.dumps({'answer.json':hashlib.sha256(payload.read_bytes()).hexdigest()}))
            self.assertEqual(delivery_files(root)['files'],1)
            payload.write_bytes(b'changed answer')
            with self.assertRaisesRegex(ValueError,'file changed'):
                delivery_files(root)
            manifest.write_text(json.dumps({'../outside.json':'0'*64}))
            with self.assertRaisesRegex(ValueError,'escapes root'):
                delivery_files(root)

    def test_r22_rejects_wrong_scalar_batch(self):
        import hashlib
        import shutil
        from checker.questions_r13 import check_readings
        with tempfile.TemporaryDirectory() as scratch:
            root=Path(scratch);packet=root/'packet';packet.mkdir()
            batch=packet/'batches.json'
            batch.write_text(json.dumps([dict(function='OTHER',path='requests-BESSELI.json',calls=49)]))
            digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
            (packet/'manifest.json').write_text(json.dumps(dict(format='dnacalc-research-questions-r22/1',
                scalar_calls=49,files_sha256={'batches.json':digest(batch)})))
            result=root/'result';result.mkdir();shutil.copytree(packet,result/'packet')
            (result/'packet-sha256.json').write_text(json.dumps({p.name:digest(p)for p in packet.iterdir()}))
            with self.assertRaisesRegex(ValueError,'R22 scalar batch scope'):
                check_readings(packet,result,root/'out.json')

    def test_nested_retry_is_unique_and_hash_bound(self):
        import hashlib
        from checker.questions_r13 import frozen_prediction_path
        with tempfile.TemporaryDirectory() as scratch:
            root=Path(scratch);path=root/'retry/worker/frozen/model.json'
            path.parent.mkdir(parents=True);path.write_bytes(b'constructed model')
            digest=hashlib.sha256(path.read_bytes()).hexdigest()
            self.assertEqual(frozen_prediction_path(root,'model.json',digest),path)
            with self.assertRaisesRegex(ValueError,'prediction changed'):
                frozen_prediction_path(root,'model.json','0'*64)
            with self.assertRaisesRegex(ValueError,'invalid frozen worker name'):
                frozen_prediction_path(root,'../model.json',digest)
            duplicate=root/'worker/frozen/model.json';duplicate.parent.mkdir(parents=True)
            duplicate.write_bytes(path.read_bytes())
            with self.assertRaisesRegex(ValueError,'ambiguous'):
                frozen_prediction_path(root,'model.json',digest)


class MatchBitIntegrity(unittest.TestCase):
    def test_inclusive_counts_losses_and_changed_identity(self):
        import sys
        from racer.audit import bank_match_bits
        import racer.corpus_cases as corpus_cases
        with tempfile.TemporaryDirectory() as scratch,patch.dict(sys.modules,{'corpus_cases':corpus_cases}):
            root=Path(scratch);binary=root/'judge';data=root/'rows.json'
            binary.write_text('''#!/usr/bin/env python3
import json,sys
from pathlib import Path
rows=json.loads(Path(__file__).with_name('rows.json').read_text())
with open(sys.argv[sys.argv.index('--trace')+1],'w') as out:
    for row in rows:out.write(json.dumps(row)+'\\n')
matches=sum(r['actual']==r['expected'] for r in rows)
print(f'xtask witness: constructed control; {matches} of {len(rows)+1} rows matched, {len(rows)+1-matches} failed (mismatched, refused or unread; 1 unread)')
''')
            binary.chmod(0o755)
            rows=[dict(function='CONTROL',args=['0x3ff0000000000000'],id=str(i),
                expected='0x3ff0000000000000',actual=value,source='constructed control',
                reading={},result_index=None) for i,value in enumerate(['0x3ff0000000000000','0x4000000000000000'])]
            data.write_text(json.dumps(rows))
            before=bank_match_bits(binary,root,root/'before')
            self.assertEqual((before['total'],before['matched'],before['failed'],before['unread']),(3,1,2,1))
            rows[0]['actual'],rows[1]['actual']=rows[1]['actual'],rows[0]['actual']
            data.write_text(json.dumps(rows))
            after=bank_match_bits(binary,root,root/'after',root/'before')
            self.assertEqual((after['gains'],after['losses']),(1,1))
            rows[0]['args']=['0x4000000000000000'];data.write_text(json.dumps(rows))
            with self.assertRaisesRegex(ValueError,'input/answer/context stream changed'):
                bank_match_bits(binary,root,root/'corrupt',root/'before')
            self.assertFalse((root/'corrupt/summary.json').exists())


class DomainQuestions(unittest.TestCase):
    def plan(self,negative=False):
        return dict(brackets=[dict(id='synthetic',function='BESSELJ',order=domain.bits(0),negative=negative,
            low=domain.bits(2),high=domain.bits(4))],controls=[],maximum_rounds=64,
            maximum_endpoint_expansions=8,neighbor_radius=2,transition_assumption='SYNTHETIC single transition')

    def search(self,plan,threshold):
        def ask(calls,stage):
            return {(f,args):domain.bits(1) if abs(domain.number(args[0])) <= threshold else 'error:Num'
                    for f,args in calls}
        return domain.search(plan,ask,lambda *unused:None)

    def test_both_signs_find_adjacent_binary64_values(self):
        for negative in (False,True):
            result = self.search(self.plan(negative),3.14159)
            bracket = result['brackets'][0]
            self.assertEqual(bracket['low'],domain.bits(3.14159))
            self.assertEqual(int(bracket['high'],16)-int(bracket['low'],16),1)

    def test_a_numeric_high_endpoint_expands(self):
        result = self.search(self.plan(),6.0)
        self.assertEqual(result['brackets'][0]['low'],domain.bits(6))

    def test_missing_answers_stop_the_search(self):
        with self.assertRaisesRegex(ValueError,'incomplete adaptive batch'):
            domain.search(self.plan(),lambda calls,stage:{},lambda *unused:None)

    def test_neighborhood_reversal_stops_the_search(self):
        def ask(calls,stage):
            result = {(f,args):domain.bits(1) if abs(domain.number(args[0])) <= 3 else 'error:Num'
                      for f,args in calls}
            if stage == 'neighbors':
                for key in result:
                    if int(key[1][0],16) == int(domain.bits(3),16)+3:
                        result[key] = domain.bits(1)
            return result
        with self.assertRaisesRegex(ValueError,'nonmonotone'):
            domain.search(self.plan(),ask,lambda *unused:None)

    def fixture(self):
        args = [domain.bits(3),domain.bits(0)]
        request = dict(function='BESSELJ',_sha256='synthetic',witnesses=[dict(id='test',args=args)])
        raw_args = [dict(number=a[2:])for a in args]
        doc = dict(function='BESSELJ',unread=[],witnesses=[dict(id='test',args=args,expected_bits='error:Num')],
            capture_provenance=dict(frozen_model_predictions=dict(inputs_sha256='synthetic'),
                                    environment=dict(excel_input_plumbing='cell_value2_bulk')),
            input_audits=[dict(id='test',reading_failure=None,observed_raw=dict(error_code=2036),
                audit=dict(requested=raw_args,submitted_value2=raw_args,read_back=raw_args,
                    transport='cell_value2_bulk',argument_literals=None,formula='=BESSELJ(RC[1],RC[2])',
                    effective_callee_inputs='not directly observable'))])
        doc['transport_audit_readings'] = []
        doc['witnesses'][0]['calc_func_prediction'] = domain.bits(1)
        env = doc['capture_provenance']['environment']
        env.update(excel_version='SYNTHETIC',excel_build='SYNTHETIC',excel_file_version='SYNTHETIC',
                   cpu_id='SYNTHETIC',excel_operating_system='SYNTHETIC',excel_bitness='64-bit')
        provenance = doc['capture_provenance']
        provenance.update(automation_security_before=1,automation_security_after=1,calculation_settings={})
        frozen = dict(function='BESSELJ',input_document={k:v for k,v in request.items() if k != '_sha256'},
            inputs_sha256='synthetic',predictions=[dict(id='test',args=args,calc_func_prediction=domain.bits(1))])
        submitted = dict(function='BESSELJ',rows=[dict(id='test',args=raw_args,argument_literals=None,
                                                     formula='=BESSELJ(RC[1],RC[2])')])
        worker = dict(fatal=None,failures=[],automation_security_before=1,automation_security_after=1,
            calculation_settings={},excel_version='SYNTHETIC',excel_build='SYNTHETIC',excel_file_version='SYNTHETIC',
            cpu='SYNTHETIC',os_version='SYNTHETIC',excel_machine='x86-64',rows=[dict(id='test',input_kept=True,
                value=dict(error_code=2036),input_audit=copy.deepcopy(doc['input_audits'][0]['audit']))])
        return request,doc,frozen,submitted,worker

    def test_exact_canonical_audit_accepts(self):
        fixture = self.fixture()
        self.assertEqual(domain.verify_reading(*fixture),{tuple(fixture[0]['witnesses'][0]['args']):'error:Num'})

    def test_changed_operand_and_corrupt_result_and_duplicate_reject(self):
        request,doc,frozen,submitted,worker = self.fixture()
        for corruption in ('operand','result','duplicate','hash','missing'):
            changed = copy.deepcopy(doc)
            if corruption == 'operand':
                changed['input_audits'][0]['audit']['read_back'] = [dict(number='0000000000000000')]*2
            elif corruption == 'result':
                changed['witnesses'][0]['expected_bits'] = domain.bits(1)
            elif corruption == 'duplicate':
                changed['witnesses'].append(copy.deepcopy(changed['witnesses'][0]))
            elif corruption == 'hash':
                changed['capture_provenance']['frozen_model_predictions']['inputs_sha256'] = 'corrupted'
            else:
                changed['input_audits'] = []
            with self.subTest(corruption=corruption),self.assertRaises(ValueError):
                domain.verify_reading(request,changed,frozen,submitted,worker)

    def test_original_raw_corruption_cannot_be_hidden_by_consistent_canonical_audits(self):
        for corruption in ('result','submitted','frozen','environment'):
            request,doc,frozen,submitted,worker = self.fixture()
            if corruption == 'result':
                worker['rows'][0]['value'] = dict(number=domain.bits(1)[2:])
            elif corruption == 'submitted':
                submitted['rows'][0]['args'] = []
            elif corruption == 'frozen':
                frozen['predictions'][0]['calc_func_prediction'] = domain.bits(2)
            else:
                worker['excel_build'] = 'CORRUPTED'
            with self.subTest(corruption=corruption),self.assertRaises(ValueError):
                domain.verify_reading(request,doc,frozen,submitted,worker)

    def test_runner_joins_a_retry_with_a_different_worker_timestamp(self):
        with tempfile.TemporaryDirectory() as scratch:
            root = Path(scratch)
            packet,checkout,evidence,output = [root/name for name in ('packet','application','OxFunc-Evidence','output')]
            for directory in (packet,checkout,evidence/'readings'):
                directory.mkdir(parents=True)
            plan = self.plan()
            # Start next to the synthetic transition to keep this control small.
            plan['brackets'][0].update(low=domain.bits(3),high=f'0x{int(domain.bits(3),16)+8:016x}')
            domain.write(packet/'plan.json',plan)
            domain.write(packet/'manifest.json',dict(files_sha256={'plan.json':domain.sha(packet/'plan.json')}))
            counter = 0
            def cargo(arguments,**kwargs):
                nonlocal counter
                counter += 1
                request_path = Path(arguments[arguments.index('--inputs')+1])
                request = json.loads(request_path.read_text())
                _,doc,frozen,submitted,worker = self.fixture()
                for obj,key in ((doc,'witnesses'),(doc,'input_audits'),(frozen,'predictions'),(submitted,'rows'),(worker,'rows')):
                    obj[key] = []
                for row in request['witnesses']:
                    qid,args = row['id'],row['args']
                    numeric = int(args[0],16) <= int(domain.bits(3),16)+3
                    raw = dict(number=domain.bits(1)[2:]) if numeric else dict(error_code=2036)
                    expected = domain.bits(1) if numeric else 'error:Num'
                    operands = [dict(number=v[2:])for v in args]
                    audit = dict(requested=operands,submitted_value2=operands,read_back=operands,argument_literals=None,
                        formula='=BESSELJ(RC[1],RC[2])',transport='cell_value2_bulk',effective_callee_inputs='not directly observable')
                    doc['witnesses'].append(dict(id=qid,args=args,expected_bits=expected,calc_func_prediction=domain.bits(1)))
                    doc['input_audits'].append(dict(id=qid,audit=audit,observed_raw=raw,reading_failure=None))
                    frozen['predictions'].append(dict(id=qid,args=args,calc_func_prediction=domain.bits(1)))
                    submitted['rows'].append(dict(id=qid,args=operands,argument_literals=None,formula=audit['formula']))
                    worker['rows'].append(dict(id=qid,input_kept=True,value=raw,input_audit=audit))
                frozen.update(input_document=request,inputs_sha256=domain.sha(request_path))
                frozen_path = checkout/'target/witness/frozen'/f'BESSELJ-frozen-{counter}.json'
                frozen_path.parent.mkdir(parents=True,exist_ok=True)
                domain.write(frozen_path,frozen)
                model = doc['capture_provenance']['frozen_model_predictions']
                model.update(inputs_sha256=domain.sha(request_path),frozen_path=str(frozen_path),frozen_sha256=domain.sha(frozen_path))
                worker_root = checkout/'target/witness/read'
                failed = worker_root/f'BESSELJ-different-clock-{counter}-attempt-1'
                success = worker_root/f'BESSELJ-different-clock-{counter}-attempt-2'
                failed.mkdir(parents=True)
                success.mkdir()
                domain.write(failed/'request.json',submitted)
                domain.write(failed/'result.json',dict(worker,fatal='SYNTHETIC failed attempt'))
                domain.write(success/'request.json',submitted)
                domain.write(success/'result.json',worker)
                domain.write(evidence/'readings'/f'synthetic-{counter}.json',doc)
            with patch.object(domain,'__file__',str(packet/'domain_bisection.py')),patch.object(domain.subprocess,'run',cargo):
                domain.run(checkout,evidence,output)
            result = json.loads((output/'summary.json').read_text())
            self.assertEqual(int(result['brackets'][0]['low'],16),int(domain.bits(3),16)+3)
            self.assertTrue(json.loads((output/'run-status.json').read_text())['complete'])
            self.assertEqual(len(list((output/'worker/read').glob('*/result.json'))),2*counter)


if __name__ == '__main__':
    unittest.main()
