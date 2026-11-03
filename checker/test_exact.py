"""Deliberate independent arithmetic suite; these are not Excel readings."""
from fractions import Fraction as R
import itertools
import math
import random
import struct
import unittest
import copy
import contextlib
import io
import json
from pathlib import Path
import tempfile
from unittest.mock import patch

from exact import (BINARY64 as B, EXTENDED64 as E, EXTENDED53 as P53,
                   Lattice, Span, binary64, inverse, operation, two)
from check import validate_row, quotient, product_preimage, write_result, graph_shape, graph_trace


def extended_control(q):
    """Encode a normal, exactly representable constructed Ext80 control."""
    sign=q<0;q=abs(q)
    if not q:return 'ext80-le:00000000000000000000'
    exponent=q.numerator.bit_length()-q.denominator.bit_length()
    if q<two(exponent):exponent-=1
    significand=q*two(63-exponent)
    if significand.denominator!=1:raise ValueError('control is not representable in Ext80')
    return 'ext80-le:'+(int(significand).to_bytes(8,'little')+
        (exponent+16383+(int(sign)<<15)).to_bytes(2,'little')).hex()


class ExactArithmetic(unittest.TestCase):
    def test_variance_checker_rejects_complete_empty_scope(self):
        from checker.variance import profiles, check
        with tempfile.TemporaryDirectory()as scratch,contextlib.redirect_stdout(io.StringIO()):
            root=Path(scratch);proof=root/'proof';proof.mkdir();scoped=root/'empty.jsonl';scoped.write_text('')
            ps=list(profiles('AVERAGE'))
            (proof/'summary.json').write_text(json.dumps(dict(complete=True,function='AVERAGE',
                programs=75,observations=0,executions=0,survivors=75)))
            (proof/'profiles.jsonl').write_text(''.join(json.dumps(p)+'\n'for p in ps))
            (proof/'first-kills.jsonl').write_text('');(proof/'witnesses.jsonl').write_text('')
            (proof/'survivors.json').write_text(json.dumps(list(range(75))))
            receipt=root/'must-not-exist.json'
            with self.assertRaisesRegex(ValueError,'empty'):
                check(scoped,proof,receipt)
            self.assertFalse(receipt.exists())

    def test_devsq_memoization_uses_only_identical_root_expressions(self):
        from checker.variance import evaluate, root_profile_key, MODELS
        values=[[R(1)+two(-52),R(1),R(1)-two(-53)],
                [two(390)+two(338),two(390),two(390)-two(337)],
                [two(-390),two(-389),two(-388)]]
        for index,modes in enumerate(itertools.product(MODELS,repeat=4)):
            if index%7:continue
            for algorithm,order,finish in itertools.product((0,12),(0,1,2),MODELS):
                p=dict(id=index,algorithm=algorithm,linear_order=order,square_order=order,
                       models=[*modes,finish])
                a,l,s,m=root_profile_key(p,'DEVSQ')
                canonical=dict(algorithm=a,linear_order=l,square_order=s,models=m)
                for xs in values:self.assertEqual(evaluate(p,xs,'DEVSQ'),evaluate(canonical,xs,'DEVSQ'))
                self.assertEqual(root_profile_key(p,'VAR.S'),(algorithm,order,order,tuple(p['models'])))

    def test_variance_checker_nonempty_survivor_and_first_kill_routes(self):
        from checker.variance import profiles, check
        with tempfile.TemporaryDirectory()as scratch,contextlib.redirect_stdout(io.StringIO()):
            root=Path(scratch);ps=list(profiles('AVERAGE'))
            for matched in (False,True):
                proof=root/str(matched);proof.mkdir();scoped=proof/'input.jsonl'
                target='0x3ff8000000000000'if matched else'0x3ffc000000000000'
                raw=dict(id='SYNTHETIC-not-Excel',function='AVERAGE',result_index=None,
                    args=['0x3ff0000000000000','0x4000000000000000'],expected=target)
                scoped.write_text(json.dumps(dict(trace_index=0,raw=raw,values_bits=[0x3ff0000000000000,0x4000000000000000]))+'\n')
                (proof/'profiles.jsonl').write_text(''.join(json.dumps(p)+'\n'for p in ps))
                (proof/'summary.json').write_text(json.dumps(dict(complete=True,function='AVERAGE',programs=75,observations=1,executions=75,survivors=75 if matched else 0)))
                (proof/'survivors.json').write_text(json.dumps(list(range(75))if matched else []))
                (proof/'witnesses.jsonl').write_text(json.dumps(dict(order=0,trace_index=0,killed=0 if matched else 75,remaining=75 if matched else 0))+'\n')
                kills=[]
                if not matched:
                    for p in ps:
                        encoded=extended_control(R(3,2))if p['models'][1].startswith('resident')else'binary64:3ff8000000000000'
                        kills.append(dict(candidate_id=p['id'],order=0,trace_index=0,root=encoded,got='0x3ff8000000000000',expected=target))
                (proof/'first-kills.jsonl').write_text(''.join(json.dumps(k)+'\n'for k in kills))
                receipt=proof/'checked.json';check(scoped,proof,receipt)
                result=json.loads(receipt.read_text())
                self.assertTrue(result['complete'])
                self.assertEqual(result['counts']['program_executions'],75)

    def test_rank_schur_exact_dependencies_and_zero_column_policy(self):
        from checker.regression import rank_schur
        for threshold,normalizer,zero in itertools.product(
                ('zero','epsilon','epsilon_squared','sqrt_epsilon','gamma_n','gamma_n_squared'),
                ('column','largest','first','absolute'),('keep','drop')):
            p=dict(threshold=threshold,normalizer=normalizer,zero=zero)
            columns=[list(map(R,(1,2,5,4))),list(map(R,(2,4,10,8)))]
            self.assertEqual(rank_schur(p,(columns,False,0))[0],3)
            self.assertEqual(rank_schur(p,(columns,True,0))[0],2)
            columns=[list(map(R,(1,2,5,4))),[R(0)]*4]
            self.assertEqual(rank_schur(p,(columns,False,0))[0],2 if zero=='keep'else 3)

    def test_mpfr_log_bounds_retain_direction_and_exact_zero(self):
        from checker.reference import log_bounds
        self.assertEqual(log_bounds(0x3ff0000000000000),(R(0),R(0)))
        lo,hi=log_bounds(0x4000000000000000)
        self.assertLess(lo,hi)
        self.assertTrue(R(69,100)<lo<hi<R(7,10))

    def test_scope_does_not_compare_different_selected_results(self):
        import importlib.util
        spec=importlib.util.spec_from_file_location('scope_audit',Path(__file__).resolve().parents[1]/'racer/corpus_cases.py')
        audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
        common=dict(function='LINEST',args=[[['0x3ff0000000000000']]],source='SYNTHETIC',reading={})
        rows=[dict(common,id='coef',result_index=[1,1],expected='0x3ff0000000000000'),
              dict(common,id='df',result_index=[4,2],expected='0x4000000000000000')]
        with tempfile.TemporaryDirectory()as scratch,contextlib.redirect_stdout(io.StringIO()):
            root=Path(scratch)
            for name,values in [('bank',rows),('corpus',[]),('observations',[])]:
                (root/name).write_text(''.join(json.dumps(r)+'\n'for r in values))
            audit.variance_scope(root/'bank',root/'corpus',root/'observations',root/'out','LINEST')
            self.assertEqual(json.loads((root/'out/different-answers.json').read_text()),[])

    def test_corpus_freeze_reads_ordinary_committed_capture_and_excludes_intersection(self):
        import importlib.util
        spec=importlib.util.spec_from_file_location('corpus_audit',Path(__file__).resolve().parents[1]/'racer/corpus_cases.py')
        audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
        with tempfile.TemporaryDirectory()as scratch,contextlib.redirect_stdout(io.StringIO()):
            root=Path(scratch);folder=root/'workbooks'/'ordinary_oracle.xlsx';folder.mkdir(parents=True)
            capture=root/'original.capture.json'
            cells=[dict(sheet='S',row=i,column=1)for i in(1,2)]
            number=dict(number='4005555555555555')
            capture.write_text(json.dumps(dict(environment={},formulas=[dict(c,formula=f)for c,f in zip(cells,['=AVERAGE(1,2,5)','=\u0040AVERAGE(1,2,5)'])],observations=[dict(c,excel_raw=number)for c in cells])))
            (folder/'evaluations.json').write_text(json.dumps(dict(evaluations=[dict(c,returned=dict(value=number))for c in cells])))
            (folder/'triage.json').write_text(json.dumps(dict(cells=[dict(c,verdict='matched')for c in cells])))
            (root/'summary.json').write_text(json.dumps(dict(workbooks=[dict(name='ordinary/oracle.xlsx',capture='committed capture '+str(capture))])))
            audit.freeze(root,root/'frozen','AVERAGE')
            rows=[json.loads(line)for line in(root/'frozen/cases.jsonl').read_text().splitlines()]
            self.assertEqual(len(rows),1)
            self.assertEqual(rows[0]['formula'],'=AVERAGE(1,2,5)')
            excluded=[json.loads(line)for line in(root/'frozen/excluded.jsonl').read_text().splitlines()]
            self.assertEqual(len(excluded),1)
            self.assertIn('implicit intersection',excluded[0]['reason'])

    def test_bank_delta_orders_function_prefixes_and_counts_real_losses(self):
        import importlib.util
        spec=importlib.util.spec_from_file_location('r15_audit',Path(__file__).resolve().parents[1]/'racer/audit.py')
        audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
        def row(function,arg,actual='one',source='old'):
            return dict(function=function,id=source,args=[arg],result_index=None,
                        expected='one',actual=actual,source=source)
        before=[row('ACCRINT','a'),row('ACCRINTM','a')]
        after=[row('ACCRINT','a',source='new capture'),row('ACCRINT','b'),row('ACCRINTM','a',actual='different')]
        with tempfile.TemporaryDirectory()as scratch,contextlib.redirect_stdout(io.StringIO()):
            root=Path(scratch)
            for name,rows in [('before',before),('after',after)]:
                (root/name).write_text(''.join(json.dumps(r)+'\n'for r in rows))
            result=audit.bank_delta(root/'before',root/'after',root/'result.json')['counts']
            self.assertEqual(result['shared_inputs'],2);self.assertEqual(result['added_inputs'],1)
            self.assertEqual(result['changed_outputs'],1);self.assertEqual(result['shared_losses'],1)
            self.assertEqual(result['representative_metadata_changes'],1)
            (root/'bad').write_text(''.join(json.dumps(r)+'\n'for r in reversed(before)))
            with self.assertRaisesRegex(ValueError,'order/uniqueness'):
                audit.bank_delta(root/'bad',root/'after',root/'must-not-exist.json')
            self.assertFalse((root/'must-not-exist.json').exists())

    def test_racer_empty_scope_is_rejected_with_python_optimization(self):
        import subprocess
        import sys
        root=Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as scratch:
            path=Path(scratch);source=path/'empty.jsonl';source.write_text('')
            result=subprocess.run([sys.executable,'-O',str(root/'racer/bessel.py'),str(source),str(path/'out')],
                                  capture_output=True,text=True,cwd=root,check=False,timeout=30)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('empty recurrence cohort',result.stderr)
            self.assertFalse((path/'out/summary.json').exists())

    def test_bessel_consumer_forms_and_seed_publication_bounds(self):
        from checker.bessel import calculate, seed_bounds
        for bits in ('0x3ff0000000000000','0x3ff0000000000001','0xbfe0000000000000'):
            target=binary64(bits)
            lo,hi=seed_bounds(target)
            self.assertEqual(B.round(lo),target)
            self.assertEqual(B.round(hi),target)
            self.assertLessEqual(lo,target)
            self.assertGreaterEqual(hi,target)
        for form,model in itertools.product(range(4),('strict','stored_pc64','resident_pc64','stored_pc53','resident_pc53')):
            profile=dict(shape=form,models=[model]*3)
            for function,expected in [('BESSELJ',R(1)),('BESSELY',R(1)),('BESSELI',R(-1)),('BESSELK',R(3))]:
                self.assertEqual(calculate(profile,function,R(2),R(1),R(2)),expected)

    def test_forward_bessel_boxes_keep_exact_sequences_and_reject_wrong_later_order(self):
        from checker.bessel import forward_cohorts, forward_proof
        bits = lambda x: '0x'+struct.pack('>d', float(x)).hex()
        for function in ('BESSELJ', 'BESSELY'):
            # A synthetic recurrence at x=16, seeded by 1 and3. These are
            # arithmetic controls, not values of the Bessel functions.
            rows = [dict(function=function, args=[bits(16), bits(n)], expected=bits(v),
                         result_index=None, reading='synthetic', source='no Excel evidence')
                    for n, v in enumerate((R(1), R(3), -R(5, 8), -R(101, 32)))]
            groups = forward_cohorts(rows)
            self.assertEqual(len(groups), 1)
            for form, model in itertools.product(range(5), ('strict', 'stored_pc64', 'resident_pc64', 'stored_pc53', 'resident_pc53')):
                profile = dict(form=form, models=[model]*(4 if form == 4 else 3))
                self.assertFalse(forward_proof(rows, groups[0], profile)['excluded'])
                altered = copy.deepcopy(rows)
                altered[-1]['expected'] = bits(R(1, 2))
                proof = forward_proof(altered, groups[0], profile)
                self.assertTrue(proof['excluded'])
                self.assertEqual(proof['steps'][-1]['order'], 3)

    def test_raw_moment_closures_cover_exact_variances_and_do_not_create_other_answers(self):
        from checker.variance_closure import raw_moment_result
        bits = lambda x: '0x'+struct.pack('>d', float(x)).hex()
        for function, answer in (('VAR.S', R(1)), ('VAR.P', R(2, 3))):
            row = dict(function=function, args=[[[bits(1), bits(2), bits(3)]]], expected=bits(answer))
            for fused in (False, True):
                proof = raw_moment_result(row, fused)
                self.assertEqual(len(proof['forms']), 10 if fused else 6)
                self.assertTrue(all(form['reachable'] for form in proof['forms']))
                self.assertTrue(all(form['root_count'] > 0 for form in proof['forms']))
                wrong = dict(row, expected=bits(1000))
                self.assertTrue(all(not form['reachable'] for form in raw_moment_result(wrong, fused)['forms']))

    def test_forward_bessel_checker_rejects_an_empty_scope(self):
        from checker.bessel import check_forward
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source, directory = root/'empty.jsonl', root/'proof'
            source.write_text('')
            directory.mkdir()
            (directory/'groups.json').write_text('[]')
            with self.assertRaisesRegex(ValueError, 'empty'):
                check_forward(source, directory, root/'receipt.json')
            self.assertFalse((root/'receipt.json').exists())

    def test_miller_rescaling_keeps_its_represented_states(self):
        from checker.bessel import miller_result, decode_pair
        bits = lambda x: '0x'+struct.pack('>d', float(x)).hex()
        row = dict(args=[bits(R(1, 8)), bits(100)])
        for model in ('strict', 'resident_pc64', 'stored_pc64', 'resident_pc53', 'stored_pc53'):
            p = dict(margin=64, include_boundary=True, normalization=1, models=[model]*4)
            result = miller_result(row, p)
            self.assertGreater(result['scalings'], 0)
            self.assertGreater(decode_pair(result['published']), 0)

    def test_miller_site_decoder_and_coherent_subset(self):
        from checker.bessel import miller_profile, miller_result, MODELS
        first = miller_profile(0, True, True)
        self.assertEqual(first, dict(id=0, margin=0, include_boundary=False, models=['strict']*6, normalization=0))
        last = miller_profile(6093750-1, True, True)
        self.assertEqual((last['margin'], last['include_boundary'], last['models'], last['normalization']),
                         (64, True, ['resident_pc53']*6, 2))
        row = dict(args=['0x3ffe666666666666', '0x4000000000000000'])
        for ratio, recurrence, reduction, finish in itertools.product(MODELS, repeat=4):
            old = dict(margin=8, include_boundary=False, models=[ratio, recurrence, reduction, finish], normalization=1)
            new = dict(old, models=[ratio, recurrence, recurrence, recurrence, reduction, finish])
            self.assertEqual(miller_result(row, old), miller_result(row, new))

    def test_miller_checker_rejects_duplicate_and_missing_programs(self):
        from checker.bessel import miller_profile, miller_result, check_miller
        import gzip
        raw = dict(function='BESSELJ', result_index=None,
                   args=['0x3fc0000000000000', '0x4000000000000000'], expected='0x3ff0000000000000')
        proof = miller_result(raw, miller_profile(0, True, True))
        self.assertNotEqual(proof['published'], ['1', '1'])
        record = dict(candidate_id=0, order=0, trace_index=0, proof=proof)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root/'input.jsonl'
            source.write_text(json.dumps(raw)+'\n')
            (root/'scope.json').write_text('[0]\n')
            (root/'survivors.json').write_text('[]\n')
            for count, message in ((2, 'duplicate or unknown'), (1, 'incomplete Miller accounting')):
                with gzip.open(root/'first-kills.jsonl.gz', 'wt') as out:
                    out.write((json.dumps(record)+'\n')*count)
                with self.assertRaisesRegex(ValueError, message):
                    check_miller(source, root, root/f'invalid-{count}.json', True, True)

    def test_raw_moment_receipt_rejects_changed_digest_answer_and_coverage(self):
        from checker.variance_closure import raw_moment_result, raw_moment_check
        from collections import Counter
        row = dict(function='VAR.S', args=[['0x3ff0000000000000', '0x4000000000000000',
                                         '0x4008000000000000']], expected='0x3ff0000000000000')
        for fused in (False, True):
            proof = raw_moment_result(row, fused)
            counts = Counter(observations=1)
            for form, root in enumerate(proof['forms']):
                counts['form_'+str(form)+'_unreachable'] += not root['reachable']
                counts['root_values'] += root['root_count']
            summary = dict(complete=True, fused=fused, counts=counts, by_function={'VAR.S': counts})
            # Per-function counts do not contain the global root total.
            summary['by_function'] = {'VAR.S': {k:v for k,v in counts.items() if k != 'root_values'}}
            with tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                source = root/'input.jsonl'
                source.write_text(json.dumps(row)+'\n')
                (root/'summary.json').write_text(json.dumps(summary))
                closure = root/'closures.jsonl'
                closure.write_text(json.dumps(proof)+'\n')
                with contextlib.redirect_stdout(io.StringIO()):
                    raw_moment_check(source, root, root/'valid.json', fused)
                mutations = []
                changed = copy.deepcopy(proof)
                changed['forms'][0]['roots_sha256'] = '0'*64
                mutations.append(changed)
                changed = copy.deepcopy(proof)
                changed['raw']['expected'] = '0x4000000000000000'
                mutations.append(changed)
                mutations.append(None)
                for index, changed in enumerate(mutations):
                    closure.write_text('' if changed is None else json.dumps(changed)+'\n')
                    with self.assertRaises(ValueError):
                        raw_moment_check(source, root, root/f'invalid-{index}.json', fused)

    def test_geomean_arbitrary_site_products_and_inverse_root_bounds(self):
        from racer.geomean import product_closure
        from checker.geomean import product_result
        row = dict(function='GEOMEAN', result_index=None,
                   args=['0x3ff0000000000000', '0x4000000000000000', '0x4010000000000000'],
                   expected='0x4000000000000000')
        self.assertEqual(product_closure(row), product_result(row))
        self.assertTrue(product_result(row)['reachable'])
        wrong = dict(row, expected='0x4008000000000000')
        self.assertFalse(product_result(wrong)['reachable'])
        for values in ([0.1, 0.2], [0.1, 0.2, 0.3], [2.0**80, 2.0**-80, 1.1]):
            raw = dict(row, args=['0x'+struct.pack('>d', v).hex() for v in values])
            self.assertEqual(product_closure(raw), product_result(raw))

    def test_geomean_product_scope_and_missing_tree_control(self):
        from racer.geomean import product_closure_run
        from checker.geomean import check_product_closure
        raw = dict(function='GEOMEAN', result_index=None,
                   args=['0x'+struct.pack('>d', x).hex() for x in (0.1, 0.2, 0.3)],
                   expected='0x3ff0000000000000')
        outside = dict(raw, args=['0x'+struct.pack('>d', x).hex() for x in (1e-40, 1.0)])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, proof = root/'inputs.jsonl', root/'proof'
            source.write_text(json.dumps(raw)+'\n'+json.dumps(outside)+'\n')
            with contextlib.redirect_stdout(io.StringIO()):
                product_closure_run(source, proof)
                check_product_closure(source, proof, root/'valid.json')
            counts = json.loads((proof/'summary.json').read_text())['counts']
            self.assertEqual(counts['scoped'], 1)
            self.assertEqual(counts['outside_modest_magnitudes'], 1)
            closure = proof/'closures.jsonl'
            record = json.loads(closure.read_text())
            record['products']['7'].pop()
            closure.write_text(json.dumps(record)+'\n')
            with self.assertRaises(ValueError):
                check_product_closure(source, proof, root/'invalid.json')

    def test_geomean_root_double_rounding_and_endpoint_parity(self):
        from racer.geomean import root_publication_domains as produce
        from checker.geomean import root_publication_domains as verify
        even, odd = R(1), R(1)+two(-52)
        midpoint = R(1)+two(-53)
        for count in (2, 3):
            for root in (midpoint, midpoint+two(-65), midpoint+two(-64), midpoint+two(-64)+two(-66)):
                product = root**count
                for target in (even, odd):
                    a, b = produce([product], count, target, True), verify([product], count, target, True)
                    self.assertEqual(a, b)
                    self.assertEqual(a['binary64'], a['extended53_then_binary64'])
                above = verify([product], count, even, True)
                if root == midpoint:
                    self.assertTrue(above['binary64']['reachable'])
                elif root <= midpoint+two(-64):
                    self.assertFalse(above['binary64']['reachable'])
                    self.assertTrue(above['extended64_then_binary64']['reachable'])
                else:
                    self.assertFalse(above['extended64_then_binary64']['reachable'])
            opposite = verify([(midpoint+two(-65))**count], count, odd, True)
            self.assertTrue(opposite['binary64']['reachable'])
            self.assertFalse(opposite['extended64_then_binary64']['reachable'])

    def test_geomean_staged_root_receipt_rejects_changed_stage(self):
        from racer.geomean import product_closure_run
        from checker.geomean import check_product_closure
        raw = dict(function='GEOMEAN', result_index=None,
                   args=['0x3ff0000000000000', '0x4010000000000000'], expected='0x4000000000000000')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, proof = root/'inputs.jsonl', root/'proof'
            source.write_text(json.dumps(raw)+'\n')
            with contextlib.redirect_stdout(io.StringIO()):
                product_closure_run(source, proof, True)
                check_product_closure(source, proof, root/'valid.json', True)
            path = proof/'closures.jsonl'
            record = json.loads(path.read_text())
            record['root_roundings']['extended64_then_binary64']['reachable'] = False
            path.write_text(json.dumps(record)+'\n')
            with self.assertRaisesRegex(ValueError, 'incomplete product closure'):
                check_product_closure(source, proof, root/'invalid.json', True)

    def test_stdev_four_input_closures_and_sqrt_publication(self):
        from racer.variance_closure import raw_moment_close
        from checker.variance_closure import raw_moment_result
        # The first two controls have exact sample/population variance1.
        for function, values in [('STDEV.S', [-1.5, 0.5, 0.5, 0.5]),
                                 ('STDEV.P', [-1.0, -1.0, 1.0, 1.0]),
                                 ('STDEV', [0.1, 0.2, 0.3, 0.4])]:
            row = dict(function=function, result_index=None,
                       args=['0x'+struct.pack('>d', v).hex() for v in values], expected='0x3ff0000000000000')
            a, b = raw_moment_close(row, True, True), raw_moment_result(row, True, True)
            self.assertEqual(a, b)
            self.assertEqual(len(a['forms']), 9)
            if function != 'STDEV':
                self.assertTrue(all(form['reachable'] for form in a['forms']))
            wrong = dict(row, expected='0x408f400000000000')
            self.assertFalse(any(form['reachable'] for form in raw_moment_result(wrong, True, True)['forms']))

    def test_stdev_root_receipt_rejects_changed_boundary_and_missing_form(self):
        from racer.variance_closure import raw_moment_run
        from checker.variance_closure import raw_moment_check
        row = dict(function='STDEV.P', result_index=None, not_excel_evidence=True,
                   args=['0xbff0000000000000']*2+['0x3ff0000000000000']*2,
                   expected='0x3ff0000000000000')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, proof = root/'input.jsonl', root/'proof'
            source.write_text(json.dumps(row)+'\n')
            with contextlib.redirect_stdout(io.StringIO()):
                raw_moment_run(source, proof, True, True)
                raw_moment_check(source, proof, root/'valid.json', True, True)
            path = proof/'closures.jsonl'
            original = json.loads(path.read_text())
            changed = copy.deepcopy(original)
            boundary = changed['forms'][0]['root_roundings']['extended64_then_binary64']['span']
            boundary['lower_closed'] = not boundary['lower_closed']
            for index, record in enumerate((changed, dict(original, forms=original['forms'][:-1]))):
                path.write_text(json.dumps(record)+'\n')
                with self.assertRaisesRegex(ValueError, 'raw-moment source, root set'):
                    raw_moment_check(source, proof, root/f'invalid-{index}.json', True, True)

    def test_binomial_consumers_obey_exact_dyadic_mass_ratios(self):
        from checker.binomial import step
        for form,model in itertools.product(range(6),('strict','stored_pc64','resident_pc64','stored_pc53','resident_pc53')):
            profile=dict(shape=form,models=[model]*5)
            middle=step(profile,R(1,4),R(2),R(1),R(1,2),R(1,2))
            self.assertEqual(middle,R(1,2))
            self.assertEqual(step(profile,middle,R(1),R(2),R(1,2),R(1,2)),R(1,4))

    def test_cdf_prefix_boxes_cover_exact_sums_and_exclude_other_targets(self):
        from checker.binomial import cdf_bounds
        rows=[dict(expected='0x3fd0000000000000'),dict(expected='0x3fe0000000000000'),
              dict(expected='0x3fe8000000000000')]
        group=dict(masses=[0,1],target=2)
        for order,model,clip in itertools.product(range(3),('strict','stored_pc64','resident_pc64','stored_pc53','resident_pc53'),(False,True)):
            profile=dict(order=order,model=model,clip=clip)
            self.assertFalse(cdf_bounds(rows,group,profile)['excluded'])
            wrong=copy.deepcopy(rows);wrong[2]['expected']='0x3feccccccccccccd'
            self.assertTrue(cdf_bounds(wrong,group,profile)['excluded'])

    def test_regression_coefficients_recover_exact_lines(self):
        from checker.regression import evaluate
        for algorithm,linear,squared,model in itertools.product(range(3),range(3),range(3),('strict','stored_pc64','resident_pc64','stored_pc53','resident_pc53')):
            profile=dict(algorithm=algorithm,linear_order=linear,square_order=squared,models=[model]*5)
            x=list(map(R,(1,2,3)))
            self.assertEqual(evaluate(profile,(x,list(map(R,(3,5,7))),True),1),(R(2),'0x4000000000000000'))
            self.assertEqual(evaluate(profile,(x,list(map(R,(3,5,7))),True),2),(R(1),'0x3ff0000000000000'))
            self.assertEqual(evaluate(profile,(x,list(map(R,(2,4,6))),False),2),(R(0),'0x0000000000000000'))

    def test_average_and_devsq_use_their_own_finishes(self):
        from checker.variance import profiles, evaluate, outcome
        self.assertEqual(len(list(profiles('AVERAGE'))), 75)
        for profile in profiles('AVERAGE'):
            self.assertEqual(evaluate(profile, [R(-4), R(1), R(3)], 'AVERAGE'), R(0))
            self.assertEqual(outcome(evaluate(profile, [R(-4), R(-1), R(-1)], 'AVERAGE')),
                             '0xc000000000000000')
        for profile in profiles():
            if len(set(profile['models'])) == 1:
                self.assertEqual(evaluate(profile, [R(-3), R(0), R(3)], 'DEVSQ'), R(18))

    def test_variance_published_identities_and_grammar(self):
        from checker.variance import profiles, evaluate, outcome
        ps = list(profiles())
        self.assertEqual(len(ps), 350000)
        self.assertEqual({p['id'] for p in ps}, set(range(350000)))
        # Exact symmetric controls exercise every algorithm/order with all
        # five coherent models, including the two different online updates.
        for p in ps:
            if len(set(p['models'])) != 1:
                continue
            if p['algorithm'] >= 7:
                self.assertEqual(evaluate(p, [R(-3),R(0),R(3)], 'VAR.S'), R(9))
                self.assertEqual(evaluate(p, [R(-3),R(0),R(3)], 'VAR.P'), R(6))
                self.assertEqual(evaluate(p, [R(-3),R(0),R(3)], 'STDEV.S'), R(3))
                continue
            self.assertEqual(evaluate(p, [R(1),R(2),R(3)], 'VAR.S'), R(1))
            self.assertEqual(outcome(evaluate(p, [R(0),R(2),R(4)], 'STDEV.S')), '0x4000000000000000')
            self.assertEqual(outcome(evaluate(p, [R(2),R(4),R(6)], 'VAR.P')), '0x4005555555555555')

    def test_variance_rejects_changed_input_and_full_width_root(self):
        from checker.variance import decode_scoped, verify_root
        row = dict(values_bits=[0,0], raw=dict(args=['0x0000000000000000','0x3ff0000000000000'],
                   result_index=None, expected='0x3fe0000000000000'))
        with self.assertRaisesRegex(ValueError,'input bits'):
            decode_scoped(row)
        with self.assertRaisesRegex(ValueError,'full-width root'):
            verify_root(R(1)+two(-60), 'ext80-le:0000000000000080ff3f', 'resident_pc64')

    def test_lossless_pc53_equivalence_guards_signed_zero_and_exponent_spills(self):
        from timing import lossless_pc53_nodes
        def encoded(x):return 'binary64:'+struct.pack('>d',float(x)).hex()
        graph=dict(nodes=[dict(op={'arg':0},model='strict'),dict(op={'arg':1},model='strict'),
            dict(op={'add':[0,1]},model='stored_pc53')],output=2)
        nodes=[extended_control(1+two(-60)),encoded(-1),encoded(two(-60))]
        self.assertEqual(lossless_pc53_nodes(graph,nodes),1)
        wrong=nodes.copy();wrong[2]=extended_control(two(-60))
        with self.assertRaises(ValueError):lossless_pc53_nodes(graph,wrong)
        zero=[encoded(1),encoded(-1),encoded(0)]
        with self.assertRaisesRegex(ValueError,'zero sign'):lossless_pc53_nodes(graph,zero)
        graph['nodes'][2]['op']={'div':[0,1]}
        for divisor in (R(3),two(100)):
            raw=operation(two(-1022),divisor,'div','resident_pc53');stored=B.round(raw)
            self.assertNotEqual(raw,stored)
            nodes=[encoded(two(-1022)),encoded(divisor),encoded(stored)]
            with self.assertRaises(ValueError):lossless_pc53_nodes(graph,nodes)

    def test_timing_and_factor_preimages_contain_generated_helper(self):
        from timing import preimage
        def bits(x):
            return struct.unpack('>Q',struct.pack('>d',float(x)))[0]
        rate,pv,timing,helper=R(1,8),R(7,4),R(1),-R(7,32)
        for models in itertools.product(('strict','stored_pc64','resident_pc64'),repeat=6):
            tf=operation(R(1),rate*timing,'add',models[0])
            ratio=operation(rate*timing,tf,'div',models[1])
            factor=operation(R(1),-ratio,'add',models[2])
            q=operation(pv,helper,'div',models[3])
            output=B.round(operation(operation(q,factor,'mul',models[4]),rate,'mul',models[5]))
            row=dict(leaves=list(map(bits,(rate,pv,0,timing,0))),raw=dict(expected=f'0x{bits(output):016x}'))
            self.assertTrue(preimage(row,models).contains(helper),models)
        for model,op in itertools.product(('strict','stored_pc64','resident_pc64'),('mul','div')):
            output=B.round(operation(pv,helper,op,model))
            row['raw']['expected']=f'0x{bits(output):016x}'
            interval=preimage(row,[model],op)
            self.assertTrue(interval.contains(helper),(model,op))
            wrong=helper+R(1,100)
            self.assertFalse(interval.contains(wrong),(model,op))

    def test_due_coefficient_preimages_and_exact_cache(self):
        from timing import preimage
        import sys
        sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'racer'))
        import intersect as consumer
        from intervals import number, Interval
        def bits(x):return struct.unpack('>Q',struct.pack('>d',float(x)))[0]
        rate=R.from_float(0.03125);pv=R.from_float(1234.56789)
        helper=E.round(-R(7,31));names=('strict','stored_pc64','resident_pc64')
        for form,sites in enumerate((6,6,5)):
            for models in itertools.product(names,repeat=sites):
                timing=R(1);rt=rate*timing
                if form==0:
                    tf=operation(R(1),rt,'add',models[0])
                    correction=operation(operation(rate,rt,'mul',models[1]),tf,'div',models[2])
                    factor=operation(rate,-correction,'add',models[3])
                elif form==1:
                    tf=operation(R(1),rt,'add',models[0])
                    ratio=operation(rt,tf,'div',models[1])
                    factor=operation(rate,operation(R(1),-ratio,'add',models[2]),'mul',models[3])
                else:
                    denominator=operation(operation(R(1),rate,'div',models[0]),timing,'add',models[1])
                    factor=operation(R(1),denominator,'div',models[2])
                q=operation(pv,helper,'div',models[-2])
                output=B.round(operation(q,factor,'mul',models[-1]))
                row=dict(leaves=list(map(bits,(rate,pv,0,timing,0))),raw=dict(expected=f'0x{bits(output):016x}'))
                checked=preimage(row,models,due_form=form)
                self.assertTrue(checked.contains(helper),(form,models))
                ids=[consumer.MODELS[m]for m in models]
                expr=consumer.due_rate_expr(form,ids)
                generic=consumer.inverse_expr(expr,[pv,rate,None,None,timing,R(1),R(-1)],
                    consumer.B64.preimage(Interval(output,output)))
                cached=consumer.row_preimage(row,({'family':'due_rate_factor'},expr,ids))
                self.assertEqual(generic,cached,(form,models))
                self.assertEqual(checked.json(),cached.json(),(form,models))

    def test_resident_sign_change_preserves_low_bits_and_has_no_free_constant(self):
        graph = dict(nodes=[
            dict(op={'arg':0}, model='strict'),
            dict(op={'const':{'bits':0xbff0000000000000}}, model='strict'),
            dict(op={'mul':[1,0]}, model='resident_pc64')], output=2)
        value = 1+two(-60)
        self.assertEqual(graph_trace(graph,[value])[-1],-value)
        stored = copy.deepcopy(graph);stored['nodes'][2]['model']='strict'
        self.assertEqual(graph_trace(stored,[value])[-1],-1)
        with self.assertRaises(ValueError): graph_shape(graph)
        self.assertEqual(graph_shape(graph,True),('mul',('minus_one',),('arg',0)))
        altered=copy.deepcopy(graph);altered['nodes'][1]['op']['const']['bits']=0x4000000000000000
        with self.assertRaises(ValueError): graph_shape(altered,True)

    def test_exhaustive_small_lattices(self):
        checks = 0
        for fmt in (Lattice(2, -2, 3), Lattice(3, -3, 4), Lattice(4, -2, 2)):
            table = [(k, fmt.value(k)) for k in range(-fmt.last-1, fmt.last+2)]
            for (k, x), (n, y) in zip(table, table[1:]):
                self.assertEqual(n, k+1)
                self.assertLess(x, y)
                midpoint = (x+y)/2
                for t in (x, midpoint-(y-x)/16, midpoint, midpoint+(y-x)/16, y):
                    # Brute finite table with the two ideal overflow neighbors.
                    distance = min(abs(t-v) for _, v in table)
                    nearest = [(i, v) for i, v in table if abs(t-v) == distance]
                    index, result = min(nearest, key=lambda iv: iv[0] % 2)
                    if abs(index) > fmt.last:
                        result = None
                    self.assertEqual(fmt.round(t), result)
                    if abs(k) <= fmt.last and x <= t < y:
                        self.assertEqual(fmt.floor_coordinate(t), k)
                    checks += 1
            # Independently enumerate lattice members and their midpoint cells.
            probes = [v for _, v in table]
            probes += [(a[1]+b[1])/2 for a, b in zip(table, table[1:])]
            for a, b in zip(probes[::3], probes[2::3]):
                lo, hi = sorted((a, b))
                if lo == hi:
                    continue
                for lc, hc in itertools.product((False, True), repeat=2):
                    span = Span(lo, hi, lc, hc)
                    members = [(k, v) for k, v in table if abs(k) <= fmt.last and span.contains(v)]
                    bounds = fmt.coordinates(span)
                    self.assertEqual(bounds, (members[0][0], members[-1][0]) if members else None)
                    pre = fmt.pullback(span)
                    if members:
                        first, last = members[0][0], members[-1][0]
                        self.assertEqual(pre, Span((fmt.value(first-1)+fmt.value(first))/2,
                            (fmt.value(last)+fmt.value(last+1))/2, first % 2 == 0, last % 2 == 0))
                    else:
                        self.assertIsNone(pre)
                    checks += 1
        print('exhaustive small-lattice checks:', checks)

    def test_binary64_against_raw_encoding_and_neighbors(self):
        rng = random.Random(0x524e395f43484543)
        values = [0, 1, 2, (1 << 52)-1, 1 << 52, 0x3ff0000000000000,
                  0x7fefffffffffffff] + [rng.randrange(0x7ff0000000000000) for _ in range(2048)]
        for bits in values:
            for signed in (bits, bits | (1 << 63)):
                x = struct.unpack('>d', signed.to_bytes(8, 'big'))[0]
                q = binary64(signed)
                self.assertEqual(q, R.from_float(x))
                self.assertEqual(B.round(q), q)
                k = B.floor_coordinate(q)
                self.assertEqual(B.value(k), q)
                for direction, delta in ((math.inf, 1), (-math.inf, -1)):
                    neighbor = math.nextafter(x, direction)
                    if math.isfinite(neighbor):
                        self.assertEqual(B.value(k+delta), R.from_float(neighbor))

    def test_deliberate_double_rounding_and_mixed_precision(self):
        one = Span(R(1), R(1))
        self.assertEqual(E.pullback(B.pullback(one)),
            Span(1-two(-54)-two(-65), 1+two(-53)+two(-64)))
        b = two(-53)+two(-65)
        self.assertEqual(B.round(1+b), 1+two(-52))
        self.assertEqual(B.round(E.round(1+b)), 1)
        resident = E.round(1+two(-60))
        self.assertEqual(operation(resident, R(-1), 'add', 'strict'), 0)
        self.assertEqual(operation(resident, R(-1), 'add', 'resident_pc64'), two(-60))
        self.assertEqual(operation(resident, R(-1), 'add', 'resident_pc53'), two(-60))
        wide = P53.round(two(1023)+two(1023))
        self.assertEqual(wide, two(1024))
        self.assertEqual(P53.round(wide/2), two(1023))
        self.assertIsNone(B.round(two(1024)))
        self.assertEqual(E.round(resident)-1, two(-60))
        self.assertEqual(B.round(resident)-1, 0)

    def test_open_endpoints_negative_inversion_underflow_overflow(self):
        odd = 1+two(-52)
        span = B.pullback(Span(odd, odd))
        self.assertFalse(span.left)
        self.assertFalse(span.right)
        self.assertEqual(B.round(span.low), R(1))
        self.assertEqual(B.round(span.high), 1+two(-51))
        reflected = inverse(Span(R(1), R(1)), R(-3), 'mul', True, 'resident_pc64')
        self.assertEqual(reflected, E.pullback(Span(R(1), R(1))).linear(R(-1, 3)))
        self.assertIsNone(B.coordinates(Span(R(1), 1+two(-52), False, False)))
        self.assertEqual(B.round(two(-1075)), 0)
        self.assertEqual(B.round(3*two(-1075)), two(-1073))
        self.assertEqual(B.round(-two(-1075)), 0)  # rational zero; sign is out of scope
        maximum = B.value(B.last)
        threshold = (maximum+two(1024))/2
        self.assertIsNone(B.round(threshold))
        self.assertEqual(B.round(threshold-two(970)), maximum)
        self.assertEqual(B.coordinates(Span(-two(1100), two(1100))), (-B.last, B.last))
        with self.assertRaises(ValueError):
            binary64(0x7ff0000000000000)
        with self.assertRaises(ValueError):
            Span(R(-1), R(1)).reciprocal(R(1))

    def test_original_argument_binding_and_defaults(self):
        args = ['0x3fc0000000000000', '0x4000000000000000', '0x3ff0000000000000',
                '0x8000000000000000', '0x4000000000000000']
        row = dict(raw=dict(function='PMT', result_index=None, args=args),
                   leaves=[int(args[0], 16), int(args[2], 16), int(args[3], 16),
                           0x3ff0000000000000, 0xbfd0000000000000], neg_log=0xbfd0000000000000)
        validate_row(row)
        for index in range(4):
            altered = copy.deepcopy(row)
            altered['leaves'][index] ^= 1
            with self.assertRaises(ValueError):
                validate_row(altered)
        for arity in (3, 4):
            altered = copy.deepcopy(row)
            altered['raw']['args'] = args[:arity]
            altered['leaves'][3] = 0
            if arity == 3:
                altered['leaves'][2] = 0
            validate_row(altered)
        altered = copy.deepcopy(row)
        altered['raw']['result_index'] = [1, 1]
        with self.assertRaises(ValueError):
            validate_row(altered)

    def test_quotient_trace_expected_then_actual_columns(self):
        # Constructed arithmetic control: 1.5+oneULP cannot equal DR(q*1.5)
        # for binary64 q. This is not a recorded Excel observation.
        rate_bits, expected_bits = '3ff8000000000000', '3ff8000000000001'
        rate, expected = binary64(rate_bits), binary64(expected_bits)
        span = product_preimage(expected).linear(1/rate)
        self.assertIsNone(B.coordinates(span))
        cert = dict(id='synthetic-parser-gap', rate_bits=rate_bits, expected_bits=expected_bits,
            old_actual_bits=rate_bits, old_quotient_bits='3ff0000000000000',
            lower=span.json()['lower'], upper=span.json()['upper'], closed=span.left,
            neighbor_magnitudes=['3ff0000000000000', '3ff0000000000001'])
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            (base/'raw.tsv').write_text('Q\tsynthetic-parser-gap\t'+rate_bits+
                '\t3ff0000000000000\t'+expected_bits+'\t'+rate_bits+'\n')
            (base/'certs.jsonl').write_text(json.dumps(cert)+'\n')
            with contextlib.redirect_stdout(io.StringIO()):
                quotient(base/'raw.tsv', base/'certs.jsonl', base/'result.json')
            result = json.loads((base/'result.json').read_text())
            self.assertEqual(result['counts']['unreachable_by_binary64_quotient'], 1)
            self.assertEqual(result['counts']['original_quotient_product_replays'], 1)

    def test_current_product_shell_binds_routes_and_rejects_changed_inputs(self):
        from timing import product_shells
        args=['0x3ff8000000000000','0x4000000000000000','0x3ff0000000000000',
              '0x0000000000000000','0x0000000000000000']
        raw=dict(id='SYNTHETIC-shell-gap',function='PMT',result_index=None,args=args,
            expected='0x3ff8000000000001')
        row=dict(raw=raw,leaves=[int(args[0],16),int(args[2],16),0,0,0xbff0000000000000],
            neg_log=0xbff0000000000000)
        route=dict(id=raw['id'],args=args,expected=raw['expected'],authored_echoes=args,
            **{k:'0x3ff8000000000000'for k in ('workbook','registered','direct_body','saved_reopened')})
        with tempfile.TemporaryDirectory()as temp:
            base=Path(temp);scoped=base/'scoped.jsonl';routes=base/'routes.json'
            scoped.write_text(json.dumps(row)+'\n');routes.write_text(json.dumps([route]))
            with contextlib.redirect_stdout(io.StringIO()):product_shells(scoped,routes,base/'receipt.json')
            result=json.loads((base/'receipt.json').read_text())
            self.assertEqual(result['counts'],dict(rows=1,application_matches=0,strict_unreachable=1,stored_pc64_unreachable=1))
            for field in ('args','authored_echoes','direct_body'):
                altered=copy.deepcopy(route)
                if isinstance(altered[field],list):altered[field][0]='0x3ff0000000000000'
                else:altered[field]='0x3ff0000000000000'
                routes.write_text(json.dumps([altered]));output=base/(field+'.json')
                with self.assertRaises(ValueError):product_shells(scoped,routes,output)
                self.assertFalse(output.exists())

    def test_helper_component_preserves_resident_bits_and_binds_stores(self):
        from timing import helper_component
        extended=extended_control
        eps=two(-60);low=-R(1,4)+eps
        h=dict(t='resident_pc64',exp_store=False,difference='resident_pc64',
            correction=True,log_store=True,association=0,
            arithmetic=['resident_pc64','resident_pc64'],em_store=False)
        base=dict(t=extended(-R(1,4)),u=extended(R(3,4)+eps),em=extended(low),
            exp_raw=extended(R(3,4)+eps),log_raw=extended(low))
        self.assertEqual(helper_component(h,base,R(2),R(1,8),{}),low)
        # Five distinct stores remove the retained low bits in this constructed
        # control. The elementary values are deliberate inputs, not Excel data.
        for change in ({'exp_store':True},{'difference':'strict'},{'em_store':True},
                       {'log_store':False},{'arithmetic':['strict','resident_pc64']}):
            selected={**h,**change};comp=copy.deepcopy(base);comp['em']=extended(-R(1,4))
            if change.get('exp_store'):comp['u']='binary64:3fe8000000000000'
            self.assertEqual(helper_component(selected,comp,R(2),R(1,8),{}),-R(1,4),change)
        for field in ('t','u','em'):
            altered=copy.deepcopy(base);altered[field]=extended(R(1,2))
            with self.assertRaises(ValueError):helper_component(h,altered,R(2),R(1,8),{})
        elementary={};helper_component(h,base,R(2),R(1,8),elementary)
        altered=copy.deepcopy(base);altered['exp_raw']=extended(R(3,4))
        with self.assertRaises(ValueError):helper_component(h,altered,R(2),R(1,8),elementary)

    def test_direct_helper_reduction_scaling_stores_and_tampering(self):
        from timing import direct_component,series_log1p
        self.assertEqual(series_log1p(two(-80)),two(-80))
        for bad in (R(0),R(-1,8),R(1,2)):
            with self.assertRaises(ValueError):series_log1p(bad)
        eps=two(-60);raw=-R(1,8)+eps;rom=R(0xb8aa3b295c17f0bc,1<<63)
        h=dict(models=['resident_pc64']*4,f2xm1_store=False)
        for n,log,k in ((R(2),R(1,8),R(0)),(R(1),R(1),R(-1))):
            t=-n*log;x=E.round(t*rom);reduced=E.round(x-k)
            power=two(int(k));scaled=raw*power;integer=power-1;em=scaled+integer
            values=[t,x,k,reduced,raw,raw,scaled,power,integer,em]
            comp=dict(nodes=list(map(extended_control,values)),em=extended_control(em))
            self.assertEqual(direct_component(h,comp,n,log,{}),em)
            for index in (0,2,3,6,7,9):
                altered=copy.deepcopy(comp);altered['nodes'][index]=extended_control(values[index]+R(1,2))
                with self.assertRaises(ValueError):direct_component(h,altered,n,log,{})
        # F2XM1 store loses the explicit resident low bits before reconstruction.
        selected={**h,'f2xm1_store':True};values[5]=B.round(raw)
        values[6]=values[5]*power;values[9]=values[6]+integer
        comp=dict(nodes=list(map(extended_control,values)),em=extended_control(values[9]))
        self.assertEqual(direct_component(selected,comp,n,log,{}),values[9])

    def test_source_drift_refuses_a_success_receipt(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp)/'receipt.json'
            with patch('check.LOADED_SOURCE_SHA256', {'exact.py': 'changed'}):
                with self.assertRaisesRegex(ValueError, 'source changed'):
                    write_result(output, {})
            self.assertFalse(output.exists())


class NormalReferenceTests(unittest.TestCase):
    def test_complement_lattice_discriminator(self):
        from normal import complement_sets
        p,w=two(-54),two(-66)
        exactly=complement_sets(Span(p,p))
        self.assertEqual(exactly[:3],[['0x3ff0000000000000']]*3)
        self.assertEqual(exactly[3:],[['0x3fefffffffffffff','0x3ff0000000000000']]*2)
        above=complement_sets(Span(p,p+w,False,True))
        self.assertEqual(above,[['0x3fefffffffffffff'],['0x3ff0000000000000'],['0x3ff0000000000000'],['0x3fefffffffffffff'],['0x3fefffffffffffff']])

    def test_root_midpoints_and_extended_root_store(self):
        from normal import sqrt_round, encoding
        middle = R(1)+two(-53)
        self.assertEqual(sqrt_round(middle*middle, 'resident_pc64'), middle)
        self.assertEqual(sqrt_round(middle*middle, 'stored_pc64'), R(1))
        for model in ('strict','resident_pc64','resident_pc53'):
            root = sqrt_round(R(2), model)
            fmt = B if model == 'strict' else E if model.endswith('pc64') else P53
            k = fmt.floor_coordinate(root)
            lower, upper = (fmt.value(k-1)+root)/2, (root+fmt.value(k+1))/2
            self.assertLess(lower*lower,R(2))
            self.assertGreater(upper*upper,R(2))
        self.assertEqual(encoding(two(-1023)), '0x0000000000000000')
        self.assertEqual(encoding(two(-1022)), '0x0010000000000000')

    def test_concrete_transcendental_bindings_reject_altered_input_or_result(self):
        from normal import Replay
        nodes = [dict(op='input',bits='0000000000000000',value='binary64:0000000000000000'),
                 dict(op='exp',args=[0],provider=3,concrete_leaf='3:binary64:0000000000000000',value='binary64:3ff0000000000000')]
        helper = {'3:binary64:0000000000000000':dict(provider=3,input='binary64:0000000000000000',output='binary64:3ff0000000000000')}
        m = Replay(nodes,{},helper);m.node('input',R(0),bits='0000000000000000');m.exp(0,3)
        for key,value in [('input','binary64:0000000000000001'),('output','binary64:3ff0000000000001'),('provider',2)]:
            altered = copy.deepcopy(helper);altered[next(iter(altered))][key] = value
            m = Replay(nodes,{},altered);m.node('input',R(0),bits='0000000000000000')
            with self.assertRaises(ValueError):m.exp(0,3)

    def test_directed_mathematical_reference_uses_two_precisions(self):
        from reference import normal_bounds
        self.assertEqual(normal_bounds('0000000000000000')['rounded_binary64'],'0x3fe0000000000000')
        for bits in ('3ff0000000000000','bff0000000000000','bfe537a6f4de9bd4','c042c00000000000'):
            a,b = normal_bounds(bits,192),normal_bounds(bits,384)
            self.assertEqual(a['rounded_binary64'],b['rounded_binary64'])
            self.assertLessEqual(R(*map(int,a['lower'])),R(*map(int,b['lower'])))
            self.assertGreaterEqual(R(*map(int,a['upper'])),R(*map(int,b['upper'])))


if __name__ == '__main__':
    unittest.main()
