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
    def test_dft_reports_a_known_line_and_not_a_store_count(self):
        values=[math.cos(2*math.pi*7*i/64) for i in range(64)]
        powers=[abs(z)**2 for z in f.fft(values)]
        self.assertEqual(sorted(range(64),key=lambda i:-powers[i])[:2],[7,57])
        answers=[dict(type='number',bits=f.bits(1+v*2**-45))for v in values]
        spectrum=f.spectrum(answers)
        self.assertEqual(spectrum['lines'][0]['bin'],7)
        self.assertIsNone(spectrum['resolved_stored_stages'])
        self.assertEqual(f.spectrum([dict(type='error',code='#NUM!')]*64)['status'],'non-numeric window')

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
