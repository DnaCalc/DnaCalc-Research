"""Constructed numerical/provenance controls; none are Excel observations."""
import copy
from fractions import Fraction as Q
import json
import math
from pathlib import Path
import tempfile
import unittest

from racer import fingerprint as f


class FingerprintControls(unittest.TestCase):
    def test_backward_bond_newton_keeps_all_real_branches(self):
        from racer.fable_followup import inverse_newton_interval
        import mpmath as mp
        with mp.workdps(70):
            r=mp.mpf(1)/4
            # At the double root the third preimage is negative; keeping
            # only positive predecessors loses the observed negative basin.
            intervals=inverse_newton_interval((mp.mpf(2),mp.mpf(2)),r,2)
            for seed in [mp.mpf(2),mp.mpf(-4)]:
                self.assertTrue(any(a<=seed<=b for a,b in intervals))
            for n in [2,4,8]:
                for seed in [mp.mpf(-5),mp.mpf('-1.25'),mp.mpf('.125'),
                             mp.mpf('1.25'),mp.mpf(5)]:
                    target=((n+1)*seed-r*seed**(n+1))/n
                    intervals=inverse_newton_interval((target,target),r,n)
                    self.assertTrue(any(a<=seed<=b for a,b in intervals), (n,seed))
                    self.assertTrue(all(a<=b for a,b in intervals))

    def test_corrected_exponent_inverse_keeps_the_constructed_normalizer(self):
        from racer.gamma_stages import corrected_exponent_normalizer
        from exact import operation, EXTENDED64 as E
        # The EXP leaf is arbitrary: inversion tests only the measured-leaf
        # interface, independent of any exp implementation or coefficient.
        for start in (Q(-448),Q(-450)+Q(1,2**40)):
            normalizer=Q(13)+Q(1,2**55)
            ss=operation(start,normalizer,'add','stored_pc64')
            dropped=operation(start-ss,normalizer,'add','stored_pc64')
            factor=operation(Q(1),dropped,'add','stored_pc64')
            exp=Q(1,2**600)
            z=operation(exp,factor,'mul','stored_pc64')
            answer=operation(z,Q(9),'div','stored_pc64')
            span=corrected_exponent_normalizer(start,ss,exp,answer,[('div',Q(9))])
            self.assertIsNotNone(span)
            self.assertTrue(span.contains(normalizer))
            self.assertFalse(span.contains(normalizer+Q(1,2**30)))
            for product_model in ('stored_pc64','resident_pc64'):
                for scale_model in ('stored_pc64','resident_pc64'):
                    z=operation(exp,factor,'mul',product_model)
                    scaled=operation(z,Q(9),'div',scale_model)
                    answer=operation(scaled,Q(5,4),'mul','stored_pc64')
                    suffix=[('div',Q(9),scale_model),('mul',Q(5,4))]
                    span=corrected_exponent_normalizer(start,ss,exp,answer,suffix,
                        product_model=product_model)
                    self.assertIsNotNone(span)
                    self.assertTrue(span.contains(normalizer))
        with self.assertRaises(ValueError):
            corrected_exponent_normalizer(Q(-1),Q(1),Q(1),Q(1),[])

    def test_stirling_proof_offline_never_asks_for_a_missing_reply(self):
        from unittest.mock import patch
        from types import SimpleNamespace
        from racer import gamma_stages as gamma
        class NoLiveClient:
            def __init__(self, **options):pass
            def ask(self, *args, **options):
                raise AssertionError('offline proof attempted a live Excel question')
        client=SimpleNamespace(Client=NoLiveClient,number=lambda raw:dict(type='number',bits=raw))
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            with patch.object(gamma,'gammaln_stirling_shape',return_value={}), \
                 patch.object(gamma,'Native',return_value=object()), \
                 patch.object(gamma,'load_service',return_value=(client,None)):
                with self.assertRaises(FileNotFoundError):
                    gamma.gammaln_stirling_proof('.',root/'proof.json',root/'cache',
                        'unused.so',count=2,offline=True)

    def test_final_inverse_subtraction_grid_keeps_constructed_members(self):
        from racer.normal_quantile import final_subtraction_lattice
        from exact import BINARY64 as B, EXTENDED64 as E
        for t in (.52,.58,.64,.7,.706):
            # Binary64 C and R in adjacent binades, independently rounded
            # by rational arithmetic; both products use the two roundings.
            c=Q(8.);r=Q(7.25);u=c-r
            y=B.round(E.round(Q(t)*u))
            answer=float(B.round(E.round(y*Q(math.sqrt(2.)))))
            self.assertTrue(final_subtraction_lattice(t,answer,-50))

    def test_earlier_cancellation_grid_keeps_members_and_rejects_half_steps(self):
        from racer.normal_quantile import cancellation_grid
        from exact import operation
        for t in (.52,.58,.64,.7,.706):
            factors=[('mul',t*t),('div',1-t*t),('mul',t),('mul',math.sqrt(2.))]
            u=Q(3,4)
            for op,constant in factors:
                u=operation(u,Q(constant),op,'stored_pc64')
            self.assertTrue(cancellation_grid(t,float(u),-50,factors))
        self.assertTrue(cancellation_grid(.5,1.,-50,[]))
        self.assertFalse(cancellation_grid(.5,math.nextafter(1.,2.),-50,[]))
        with self.assertRaises(ValueError):
            cancellation_grid(.5,1.,-50,[('div',0.)])

    def test_phase_reader_preserves_fractional_cycles_at_large_phases(self):
        from decimal import Decimal as D, localcontext
        import numpy as np
        track=f.PhaseTrack('constructed',1.3,.2,0.)
        with localcontext() as ctx:
            ctx.prec=70
            track.coefficients[0]=D(2)**128+D('.125')
            self.assertEqual(track.evaluate(1.3) % 1,D('.125'))
        # Actual and predicted frequencies differ; the known phase shifts
        # must be read with the correct sign and fitted across the ladder.
        dx=math.ulp(1.3)
        actual=D('.200003')/D(dx)
        for center in (0,6000,-6000,14000,-14000,24000,-24000):
            xs=[1.3+(center+k-512)*dx for k in range(1024)]
            phases=[float((actual*(D(x)-D(1.3))+D('.125')) % 1) for x in xs]
            residual=np.cos(2*np.pi*np.array(phases))
            reading=track.observe(xs,residual)
            self.assertEqual(reading['status'],'read')
            self.assertLess(reading['max_cycles'],.002)
        self.assertLess(abs(float((track.evaluate(1.3,1)-actual)*D(dx))),1e-7)
        self.assertIn('aliases',track.record()['limits'])

    def test_short_noisy_phase_windows_cannot_invent_high_derivatives(self):
        from decimal import Decimal as D
        import numpy as np
        track=f.PhaseTrack('noisy control',1.3,.2,0.)
        dx=math.ulp(1.3)
        for center in (0,6000,-6000,14000,-14000,24000,-24000,50000,-50000):
            xs=[1.3+(center+k-512)*dx for k in range(1024)]
            phases=[float((D('.2')*(D(x)-D(1.3))/D(dx)) % 1) for x in xs]
            residual=np.cos(2*np.pi*(np.array(phases)+.02*math.sin(center)))
            track.observe(xs,residual)
            self.assertEqual(track.active_degree,1)
        self.assertEqual(track.coefficients[2:],[D(0)]*6)
        before=copy.deepcopy(track.record())
        track.demodulate(xs,residual)
        self.assertEqual(track.record(),before)

    def test_phase_window_keeps_lattice_bits_and_rejects_binade_crossings(self):
        import numpy as np
        answers=[dict(type='number',bits=f'{0x3ff0000000000000+3*k:016x}') for k in range(1024)]
        self.assertLess(float(np.max(np.abs(f.window_residual(answers)))),1e-10)
        answers[0]=dict(type='number',bits=f.bits(.99))
        with self.assertRaisesRegex(ValueError,'binade'):
            f.window_residual(answers)
        answers[0]=dict(type='number',bits=f.bits(-1.))
        with self.assertRaisesRegex(ValueError,'binade'):
            f.window_residual(answers)

    def test_study_bindings_cannot_replace_the_loaded_producer(self):
        source=next(iter(f.PRODUCER_BINDINGS))
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/'study.json'
            with self.assertRaisesRegex(ValueError,'conflicting analysis source binding'):
                f.write(output,dict(format='dnacalc-test/1',source_bindings={source:'0'*64}))
            self.assertFalse(output.exists())
            extra=Path(directory)/'extra.py'
            extra.write_text('# constructed source control\n')
            original=f.sha(extra)
            extra.write_text('# changed source control\n')
            with self.assertRaisesRegex(ValueError,'analysis source changed'):
                f.write(output,dict(format='dnacalc-test/1',source_bindings={str(extra):original}))
            self.assertFalse(output.exists())
            f.write(output,dict(format='dnacalc-test/1',source_bindings={str(extra):f.sha(extra)}))
            self.assertEqual(json.loads(output.read_text())['producer_bindings'],
                             {**f.PRODUCER_BINDINGS,str(extra):f.sha(extra)})

    def test_quotient_summary_cannot_repeat_a_measured_group(self):
        from checker.geomean import quotient_summary_coverage
        models=[dict(quotient=q,published=p,groups=[dict(product_power=1)])
                for q in ('divide','stored_reciprocal')for p in (False,True)]
        groups={(m['quotient'],m['published'],1)for m in models}
        quotient_summary_coverage(models,groups)
        changed=copy.deepcopy(models);changed[0]['groups']*=2
        with self.assertRaisesRegex(ValueError,'repeated'):
            quotient_summary_coverage(changed,groups)
        with self.assertRaisesRegex(ValueError,'model order'):
            quotient_summary_coverage(models[::-1],groups)

    def test_independent_primitive_grammar_rejects_changed_store_metadata(self):
        from checker.geomean import check_stage_graph_metadata,primitive_stage_grammar
        doc=dict(LN_ROM80=f.primitive_stage_graph('LN'),LN_constant53=f.primitive_stage_graph('LN','binary64'),
            EXP_ROM80_resident_reduction=f.primitive_stage_graph('EXP'),
            EXP_constant53_resident_reduction=f.primitive_stage_graph('EXP','binary64'),
            EXP_ROM80_stored_reduction=f.primitive_stage_graph('EXP',stored_reduction=True),
            EXP2=f.primitive_stage_graph('EXP2'),LOG2=f.primitive_stage_graph('LOG2'))
        check_stage_graph_metadata(doc)
        bad=copy.deepcopy(doc);bad['EXP2']['nodes'][-2]['site']=0;bad['EXP2']['sites']=1
        with self.assertRaisesRegex(ValueError,'graph metadata'):
            check_stage_graph_metadata(bad)
        bad=copy.deepcopy(doc);bad['LN_ROM80']['outputs']=[0,0]
        with self.assertRaisesRegex(ValueError,'graph metadata'):
            check_stage_graph_metadata(bad)
        for kind in ('EXP','EXP2'):
            for store in ('fraction','f2xm1_return'):
                graph=f.primitive_stage_graph(kind)
                next(n for n in graph['nodes']if n['op']==('sub'if store=='fraction'else'f2xm1'))['site']=0
                graph['sites']=1
                self.assertEqual(graph,primitive_stage_grammar(kind,spill=store))

    def test_unknown_seed_cannot_change_on_an_exact_identity_multiply(self):
        for seed in (Q(1),Q(1,2),Q(3,2)):
            self.assertIsNotNone(f.multiplication_seed_interval(seed,seed,Q(1)))
            neighbor=Q(math.nextafter(float(seed),0.))
            self.assertIsNone(f.multiplication_seed_interval(seed,neighbor,Q(1)))
        self.assertIsNotNone(f.multiplication_seed_interval(Q(3,2),Q(3),Q(2)))
        with self.assertRaisesRegex(ValueError,'positive multiplier'):
            f.multiplication_seed_interval(Q(1),Q(1),Q(0))

    def test_rank_packet_refuses_a_bracket_wider_than_its_bisection_budget(self):
        rows=[dict(id='low',changed_operand='4020000000000001',df=3),
              dict(id='high',changed_operand=f'{0x4020000000000000+2**37:016x}',df=2)]
        doc=dict(observations=rows,necessary_intervals=[dict(normalizer='column',reverse=False,
            lower_witness='low',upper_witness='high',function='LINEST',constant=False)])
        with tempfile.TemporaryDirectory()as directory:
            path=Path(directory)/'intake.json';path.write_text(json.dumps(doc))
            with self.assertRaisesRegex(ValueError,'too wide'):
                f.rank_bisection_questions(path)

    def test_independent_stage_check_rejects_dropped_pairs_and_missing_arities(self):
        from checker.geomean import stage_observation_coverage
        number=lambda x:dict(type='number',bits=f.bits(x))
        answers={('GEOMEAN',json.dumps([number(x)]+[number(1.)]*(n-1),sort_keys=True)):number(2.)
                 for x in (1.3,7.)for n in (2,4,8)}
        doc=dict(primitive_checks=[dict(program='declared')],
                 root_intervals=[dict(input=f.bits(x),primitive=k)for x in (1.3,7.)for k in ('EXP','EXP2')])
        stage_observation_coverage(doc,answers,['declared'])
        removed=copy.deepcopy(doc);removed['root_intervals']=removed['root_intervals'][2:]
        with self.assertRaisesRegex(ValueError,'missing, extra or duplicated'):
            stage_observation_coverage(removed,answers,['declared'])
        missing=dict(answers);missing.pop(next(k for k in missing if len(json.loads(k[1]))==8))
        with self.assertRaisesRegex(ValueError,'missing original root arity'):
            stage_observation_coverage(doc,missing,['declared'])
        with self.assertRaisesRegex(ValueError,'primitive programmes'):
            stage_observation_coverage(doc,answers,['declared','omitted'])

    def test_batched_edges_finish_at_observed_adjacent_inputs(self):
        class Service:
            def ask(self,name,rows):
                return dict(rows=[dict(type='boolean',value=int(r[0],16)>=r[1])for r in rows])
        brackets=[dict(low='0000000000000001',high='0000000000000064',threshold=t)for t in (23,79)]
        edges=f.bisect_observed_edges(Service(),'constructed',brackets,
            lambda s,x:[x,s['threshold']],lambda s,x,a:a['value'])
        for edge,t in zip(edges,(23,79)):
            self.assertEqual((int(edge['low'],16),int(edge['high'],16)),(t-1,t))
            seen={r['input']:r['side']for r in edge['trace']}
            self.assertIs(seen[edge['low']],False)
            self.assertIs(seen[edge['high']],True)
        with self.assertRaisesRegex(ValueError,'cannot distinguish'):
            f.bisect_observed_edges(Service(),'constructed',brackets,
                lambda s,x:[x,s['threshold']],lambda s,x,a:None)

    def test_extended_preimage_keeps_tie_parity_and_empty_images(self):
        from checker.x87 import encode,decode,bits64
        class Identity:
            def evaluate(self,graph,inputs,cw,extended):return extended*2
            def publish(self,raw):return bits64(f.BINARY64.round(decode(raw)))
        for x in (1.3,1.5):
            cell=f.positive_extended_cell(Identity(),None,f.bits(x),Q(x))
            raw=encode(Q(x));middle=int.from_bytes(raw[:8],'little')
            radius=1023 if int(f.bits(x),16)&1 else 1024
            self.assertEqual(int.from_bytes(bytes.fromhex(cell['low'])[:8],'little'),middle-radius)
            self.assertEqual(int.from_bytes(bytes.fromhex(cell['high'])[:8],'little'),middle+radius)
            self.assertEqual(int(cell['predecessor_output'],16),int(f.bits(x),16)-1)
            self.assertEqual(int(cell['successor_output'],16),int(f.bits(x),16)+1)
        class SkipsOdd(Identity):
            def publish(self,raw):return super().publish(raw)&~1
        self.assertIsNone(f.positive_extended_cell(SkipsOdd(),None,f.bits(1.3),Q(1.3)))

    def test_decreasing_preimage_keeps_original_neighbor_outputs(self):
        from checker.x87 import encode,decode,bits64
        class Reciprocal:
            def evaluate(self,graph,inputs,cw,extended):
                return encode(f.EXTENDED64.round(1/decode(extended)))*2
            def publish(self,raw):return bits64(f.BINARY64.round(decode(raw)))
        for x in (1.,1.3):
            model=Reciprocal();answer=f.bits(1/x)
            cell=f.positive_extended_cell(model,None,answer,Q(x),decreasing=True)
            self.assertIsNotNone(cell)
            self.assertGreater(int(cell['predecessor_output'],16),int(answer,16))
            self.assertLess(int(cell['successor_output'],16),int(answer,16))
            for endpoint in ('low','high'):
                self.assertEqual(model.publish(model.evaluate(None,[],0,bytes.fromhex(cell[endpoint]))[:10]),int(answer,16))

    def test_resident_output_span_preserves_open_edges_and_decreasing_images(self):
        from checker.exact import Span
        from checker.x87 import encode,decode,bits64
        class Model:
            def __init__(self,reciprocal=False):self.reciprocal=reciprocal
            def evaluate(self,graph,inputs,cw,extended):
                x=decode(extended)
                return encode(f.EXTENDED64.round(1/x)if self.reciprocal else x)*2
            def publish(self,raw):return bits64(f.BINARY64.round(decode(raw)))
        center=Q(1.3);k=f.EXTENDED64.floor_coordinate(center)
        for reciprocal in (False,True):
            model=Model(reciprocal)
            y=f.EXTENDED64.round(1/center)if reciprocal else center
            j=f.EXTENDED64.floor_coordinate(y)
            span=Span(f.EXTENDED64.value(j-7),f.EXTENDED64.value(j+11),False,True)
            cell=f.positive_extended_cell(model,None,None,center,decreasing=reciprocal,output_span=span)
            feasible=[i for i in range(k-100,k+101)
                      if span.contains(decode(model.evaluate(None,[],0,encode(f.EXTENDED64.value(i)))[:10]))]
            self.assertEqual(decode(bytes.fromhex(cell['low'])),f.EXTENDED64.value(min(feasible)))
            self.assertEqual(decode(bytes.fromhex(cell['high'])),f.EXTENDED64.value(max(feasible)))
            self.assertFalse(span.contains(decode(bytes.fromhex(cell['predecessor_output_extended']))))
            self.assertFalse(span.contains(decode(bytes.fromhex(cell['successor_output_extended']))))
        delta=f.EXTENDED64.value(k+1)-center
        empty=Span(center+delta/4,center+3*delta/4)
        self.assertIsNone(f.positive_extended_cell(Model(),None,None,center,output_span=empty))

    def test_dft_reports_a_known_line_and_not_a_store_count(self):
        values=[math.cos(2*math.pi*7*i/64) for i in range(64)]
        powers=[abs(z)**2 for z in f.fft(values)]
        self.assertEqual(sorted(range(64),key=lambda i:-powers[i])[:2],[7,57])
        answers=[dict(type='number',bits=f.bits(1+v*2**-45))for v in values]
        spectrum=f.spectrum(answers)
        self.assertEqual(spectrum['lines'][0]['bin'],7)
        self.assertIsNone(spectrum['resolved_stored_stages'])
        self.assertEqual(f.spectrum([dict(type='error',code='#NUM!')]*64)['status'],'non-numeric window')

    def test_spectrum_reports_a_root_window_without_normalizing_by_zero_ulp(self):
        # GAMMALN at its exact roots has a zero midpoint. Dividing surrounding
        # finite answers by ulp(0) used to overflow before the DFT ran.
        for values in ((-2**-50, 0., 2**-50), (0., 0., 0.),
                       (-2**-50, 2**-50)):
            result=f.spectrum([dict(type='number',bits=f.bits(x))for x in values])
            self.assertEqual(result['status'],'window touches or crosses zero')
            self.assertEqual(result['lines'],[])
            self.assertIsNone(result['resolved_stored_stages'])

    def test_mathematical_publication_discriminates_a_double_rounding_tie(self):
        # Just above a binary64 midpoint, but inside an extended64 rounding
        # cell centered at that midpoint. This distinguishes RN53(real) from
        # RN53(RN64(real)) without using any captured expected value.
        midpoint=Q(1)+Q(1,2**53)
        real=midpoint+Q(1,2**66)
        single,double=f.publication((real,real))
        self.assertEqual(single,'3ff0000000000001')
        self.assertEqual(double,'3ff0000000000000')
        with self.assertRaisesRegex(ValueError,'cannot certify'):
            f.publication((midpoint-Q(1,2**66),real))

    def test_period_two_line_is_retained_at_nyquist(self):
        answers=[dict(type='number',bits=f.bits(1+(-1 if i%2 else 1)*2**-50))for i in range(64)]
        spectrum=f.spectrum(answers)
        self.assertEqual(spectrum['lines'][0]['bin'],32)
        self.assertEqual(spectrum['lines'][0]['cycles_per_input'],.5)
        self.assertIsNone(spectrum['resolved_stored_stages'])

    def test_mpfr_references_enclose_exact_identities(self):
        number=lambda x:dict(type='number',bits=f.bits(x))
        for name,args,exact in [('EXP',[0.],Q(1)),('LN',[1.],Q(0)),
                                ('POWER',[2.,3.],Q(8)),('GEOMEAN',[1.,8.],None),
                                ('SQRT',[4.],Q(2)),('STANDARDIZE',[3.,1.,2.],Q(1))]:
            bounds=f.elementary_bounds(name,list(map(number,args)))
            if exact is not None:self.assertLessEqual(bounds[0],exact);self.assertGreaterEqual(bounds[1],exact)
            self.assertIsNotNone(f.publication(bounds))
        args=list(map(number,[2.,4.,.5,0.]))
        self.assertEqual(f.elementary_bounds('BINOM.DIST',args),(Q(3,8),Q(3,8)))
        self.assertIsNone(f.elementary_bounds('LINEST',[number(1.)]))

    def test_exact_complement_is_not_prepublished_to_binary64(self):
        argument=lambda x:[dict(type='number',bits=f.bits(x))]
        x=math.nextafter(.13,1.)
        self.assertNotEqual(Q(1)-Q(x),Q(1-x))
        exact=f.elementary_bounds('LOG1P',argument(-x))
        published=f.elementary_bounds('LN',argument(1-x))
        self.assertNotEqual(exact,published)
        self.assertEqual(f.elementary_bounds('LOG1P',argument(-.5)),
                         f.elementary_bounds('LN',argument(.5)))
        self.assertEqual(f.elementary_bounds('LOG1P',argument(0.)),(Q(0),Q(0)))
        self.assertIsNone(f.elementary_bounds('LOG1P',argument(-1.)))

    def test_readbacks_and_cached_requests_are_mandatory(self):
        rows=[[dict(type='number',bits=f.bits(1.))]]
        reply=dict(version=1,id='1'*32,ok=True,inputs=rows,rows=[dict(type='number',bits=f.bits(2.))],
                   automation_security_before=1,automation_security_after=1,
                   excel_build='constructed',machine='constructed',timestamp='constructed',
                   request_sha256='constructed',session_receipt_sha256='constructed')
        f.accept(reply,'CONTROL',rows)
        for field in ('inputs','rows','session_receipt_sha256','automation_security_before','automation_security_after','version','id'):
            changed=copy.deepcopy(reply);changed.pop(field)
            with self.assertRaises(ValueError):f.accept(changed,'CONTROL',rows)
        changed=copy.deepcopy(reply)
        changed.pop('automation_security_before');changed.pop('automation_security_after')
        with self.assertRaises(ValueError):f.accept(changed,'CONTROL',rows)
        changed=copy.deepcopy(reply);changed['inputs'][0][0]['bits']=f.bits(3.)
        with self.assertRaises(ValueError):f.accept(changed,'CONTROL',rows)
        with self.assertRaises(ValueError):f.accept(changed,'CONTROL',rows,allow_unread_inputs=True)
        changed['rows']=[dict(type='unread',reason='input-changed')]
        with self.assertRaises(ValueError):f.accept(changed,'CONTROL',rows)
        f.accept(changed,'CONTROL',rows,allow_unread_inputs=True)
        changed['rows'][0]['reason']='other'
        with self.assertRaises(ValueError):f.accept(changed,'CONTROL',rows,allow_unread_inputs=True)
        class Client:
            def ask(self,function,arguments):return reply
        with tempfile.TemporaryDirectory()as root:
            service=f.CachedService(Client(),root)
            self.assertEqual(service.ask('CONTROL',rows),reply)
            cached=f.CachedService(None,root,offline=True)
            self.assertEqual(cached.ask('CONTROL',rows),reply)
            with self.assertRaises(FileNotFoundError):cached.ask('OTHER',rows)

    def test_product_stage_checker_rejects_a_changed_mathematical_product(self):
        from checker.geomean import check_intermediate_stores
        source=f.ROOT/'evidence/round26/product-stages.json'
        data=json.loads(source.read_text())
        data['records'][0]['products']['stored_each']='0000000000000000'
        with tempfile.TemporaryDirectory()as root:
            input_path=Path(root)/'changed.json';output=Path(root)/'receipt.json'
            input_path.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError,'incorrect independently computed'):
                check_intermediate_stores(input_path,f.ROOT/'evidence/round26/fingerprints',output)
            self.assertFalse(output.exists())


if __name__=='__main__':
    unittest.main()
