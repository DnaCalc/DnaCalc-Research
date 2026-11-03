"""Adversarial native-ISA, spill, and certificate coverage controls (not Excel)."""
import copy
import ctypes as C
import importlib.util
from fractions import Fraction as Q
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(1,str(Path(__file__).resolve().parents[1]/'racer'))
from x87_programs import Builder, Native, ROOT, build_library
from x87 import Exact, Instructions, decode, encode, subtract_box, published_bits
from exact import BINARY64 as B, two
from compare_x87_probe import compare


class NativeControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.scratch=tempfile.TemporaryDirectory(prefix='r13-x87-controls-')
        root=Path(cls.scratch.name)
        build_library(root/'candidate.so')
        subprocess.run(['cc','-std=c11','-O2','-Wall','-Wextra','-Werror','-shared',
            '-fPIC',str(ROOT/'checker/x87_isa.c'),'-o',str(root/'checker.so')],check=True,timeout=120)
        cls.native=Native(root/'candidate.so')
        cls.exact=Exact(Instructions(root/'checker.so'))

    @classmethod
    def tearDownClass(cls):
        cls.scratch.cleanup()

    def test_pc53_consumes_full_resident_operand(self):
        g=Builder(2);out=g.op('sub',1,0)
        graph=g.graph(out,out);resident=encode(1+two(-60))
        actual=self.native.evaluate(graph,[0x3ff0000000000000],0x123f,extended=resident)[:10]
        self.assertEqual(decode(actual),two(-60))
        self.assertEqual(actual,self.exact.graph(graph,[0x3ff0000000000000],0x123f,0,resident)[:10])
        self.assertEqual(B.round(1+two(-60))-1,0)

    def test_explicit_store_destroys_the_extra_operand_bits(self):
        g=Builder(1);a=g.op('add',g.one,0);b=g.op('sub',a,g.one)
        graph=g.graph(b,b);inputs=[0x3c30000000000000] # 2^-60
        resident=self.native.evaluate(graph,inputs,0x133f,mask=0)[:10]
        stored=self.native.evaluate(graph,inputs,0x133f,mask=1)[:10]
        self.assertEqual(decode(resident),two(-60));self.assertEqual(decode(stored),0)

    def test_double_rounding_and_phase_word_are_observable(self):
        g=Builder(1);g.phase=1;a=g.op('add',g.one,0)
        graph=g.graph(a,a);inputs=[0x3ca0008000000000] # 2^-53 + 2^-66
        ext=self.native.evaluate(graph,inputs,0x123f,mask=1,inner_cw=0x133f)[:10]
        low=self.native.evaluate(graph,inputs,0x133f,mask=1,inner_cw=0x123f)[:10]
        self.assertEqual(decode(ext),1)
        self.assertEqual(decode(low),1+two(-52))
        for prep in (0x123f,0x133f):
            for inner in (0x123f,0x133f):
                for mask in (0,1):
                    self.assertEqual(self.native.evaluate(graph,inputs,prep,mask,inner_cw=inner),
                        self.exact.graph(graph,inputs,prep,mask,inner_cw=inner))

    def test_microcode_ignores_pc_and_retains_low_bits(self):
        one=encode(Q(1));arg=encode(Q(5,16))
        for kind,a,b in (('f2xm1',arg,one),('fyl2xp1',one,arg),('fyl2x',one,encode(Q(3,2)))):
            values=[self.native.leaf(kind,a,b,cw)for cw in (0x123f,0x133f)]
            self.assertEqual(values[0],values[1])
            self.assertNotEqual(int.from_bytes(values[0][:8],'little')&2047,0)
            for cw in (0x123f,0x133f):
                self.assertEqual(values[0],self.exact.instructions.run(kind,a,b,cw))

    def test_scale_keeps_64_bit_operand_under_pc53(self):
        a,b=encode(1+two(-60)),encode(Q(2))
        for cw in (0x123f,0x133f):
            value=self.native.leaf('fscale',a,b,cw)
            self.assertEqual(decode(value),4+two(-58))
            self.assertEqual(value,self.exact.operation('fscale',a,b,cw,False))

    def test_round_to_integer_ties_and_negative_zero(self):
        for a,expected in ((Q(-1,2),Q(0)),(Q(3,2),Q(2)),(Q(5,2),Q(2)),(Q(-5,2),Q(-2))):
            for cw in (0x123f,0x133f):
                value=self.native.leaf('frndint',encode(a),encode(Q(0)),cw)
                self.assertEqual(decode(value),expected)
                self.assertEqual(value,self.exact.operation('frndint',encode(a),encode(Q(0)),cw,False))

    def test_malformed_native_requests_abort(self):
        g=Builder(1);a=g.op('add',g.one,0);graph=g.graph(a,a)
        with self.assertRaisesRegex(ValueError,'incomplete extended'):
            self.native.evaluate(graph,[0],0x133f,extended=b'bad')
        bad=copy.deepcopy(graph);bad['nodes'][-1]['a']=1000
        with self.assertRaisesRegex(ValueError,'malformed'):
            self.native.evaluate(bad,[0],0x133f)
        with self.assertRaisesRegex(ValueError,'malformed'):
            self.native.evaluate(graph,[0],0x133f,inner_cw=0x113f)

    def test_native_full_traces_agree_with_rationals(self):
        g=Builder(1);rom=g.op('fldln2',store=False);log=g.op('fyl2xp1',rom,0)
        g.phase=1;out=g.op('mul',log,g.one);graph=g.graph(out,out)
        for prep in (0x123f,0x133f):
            for inner in (0x123f,0x133f):
                for mask in range(1<<graph['sites']):
                    self.assertEqual(self.native.evaluate(graph,[0x3fb999999999999a],prep,mask,trace=True,inner_cw=inner),
                        self.exact.graph(graph,[0x3fb999999999999a],prep,mask,inner_cw=inner,trace=True))

    def test_fused_period_log1p_uses_native_instruction_for_every_store(self):
        from x87_programs import helper_templates
        h=next(h for h in helper_templates()if h['path']==5 and h['mode']=='growth'and h['reconstruction']=='plain')
        args=[0x3f847ae147ae147b,0x4042800000000000,0x3ff0000000000000,0,0]
        for prep in (0x123f,0x133f):
            for inner in (0x123f,0x133f):
                for mask in range(1<<h['graph']['sites']):
                    self.assertEqual(self.native.evaluate(h['graph'],args,prep,mask,inner_cw=inner),
                        self.exact.graph(h['graph'],args,prep,mask,inner_cw=inner))


class CoverageControls(unittest.TestCase):
    def fixture(self):
        value=bytes(20)
        return ([15,15],[value,value],[dict(template=0),dict(template=0)],
                [dict(mode='growth')],{('growth',value):5},
                dict(helpers=[0,1],helper=value.hex(),consumers=hex(10)),set())

    def test_complete_box_subtracts_every_assignment_once(self):
        args=self.fixture();self.assertEqual(subtract_box(*args),4)
        self.assertEqual(args[0],[5,5])

    def test_omitted_consumer_is_rejected(self):
        args=self.fixture();args[5]['consumers']=hex(2)
        with self.assertRaisesRegex(ValueError,'incomplete'):subtract_box(*args)

    def test_matched_consumer_cannot_be_a_kill(self):
        args=self.fixture();args[5]['consumers']=hex(11)
        with self.assertRaisesRegex(ValueError,'includes a match'):subtract_box(*args)

    def test_unequal_helpers_cannot_share_a_box(self):
        args=self.fixture();args[1][1]=b'\x01'+bytes(19)
        with self.assertRaisesRegex(ValueError,'equality fails'):subtract_box(*args)

    def test_duplicate_helper_does_not_inflate_coverage(self):
        args=self.fixture();args[5]['helpers']=[0,0]
        with self.assertRaisesRegex(ValueError,'duplicate'):subtract_box(*args)

    def test_an_already_counted_helper_cannot_reappear(self):
        args=self.fixture();args[6].add(0)
        with self.assertRaisesRegex(ValueError,'duplicate helper'):subtract_box(*args)


class GrammarControls(unittest.TestCase):
    def test_fused_family_multiplies_inside_fyl2xp1(self):
        from x87_programs import helper_templates
        selected=[h for h in helper_templates()if h['path']==5]
        self.assertEqual(len(selected),17)
        for h in selected:
            nodes=h['graph']['nodes'];i=next(i for i,n in enumerate(nodes)if n['op']=='fyl2xp1')
            instruction=nodes[i]
            self.assertEqual(instruction['b'],0)
            multiplier=instruction['a']
            if h['mode']=='discount'and h['reconstruction']!='reciprocal':
                self.assertEqual(nodes[multiplier],dict(op='neg',a=1,b=0))
            else:self.assertEqual(multiplier,1)
            self.assertTrue(all(n['op']in ('input','const','neg')for n in nodes[:i]))

    def test_publication_cache_keeps_every_operand_bit(self):
        below=encode(1+two(-53)-two(-63))
        above=encode(1+two(-53)+two(-63))
        self.assertNotEqual(below,above)
        self.assertEqual(published_bits(below),0x3ff0000000000000)
        self.assertEqual(published_bits(above),0x3ff0000000000001)
        self.assertEqual(published_bits(encode(Q(0),True)),0x8000000000000000)

    def test_every_consumer_satisfies_the_annuity_identity_exactly(self):
        from x87_consumers import templates
        from exact import binary64
        for c in templates(Builder):
            for timing in (Q(0),Q(1)):
                growth=Q(3,2)if c['mode']=='growth'else Q(2,3)
                args=[Q(1,100),Q(10),Q(7),Q(3),timing,growth,growth-1]
                values=[]
                for n in c['graph']['nodes']:
                    kind=n['op']
                    if kind=='input':value=args[n['a']]
                    elif kind=='const':value=binary64(int(n['bits'],16))
                    elif kind=='neg':value=-values[n['a']]
                    else:
                        a,b=values[n['a']],values[n['b']]
                        value={'add':lambda:a+b,'sub':lambda:a-b,'mul':lambda:a*b,'div':lambda:a/b}[kind]()
                    values.append(value)
                amount=7*growth+3 if c['mode']=='growth'else 7+3*growth
                expected=amount*Q(1,100)/((growth-1)*(1+timing*Q(1,100)))
                if c['mode']=='growth':expected=-expected
                self.assertEqual(values[c['graph']['outputs'][0]],expected,c)
                self.assertEqual(sum(n.get('phase',0)==1 for n in c['graph']['nodes']),1)

    def test_wrong_standalone_function_graph_is_not_accepted(self):
        path=ROOT/'checker/questions_r13.py'
        spec=importlib.util.spec_from_file_location('r13_packet_checker',path)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        ln=module.standalone_graph('LN',['0x3ff8000000000000'],False)
        exp=module.standalone_graph('EXP',['0x3ff8000000000000'],False)
        self.assertNotEqual(ln,exp)
        with self.assertRaisesRegex(ValueError,'logical'):
            module.standalone_graph('LN',['0x3ff8000000000000'],'False')


class ProbeControls(unittest.TestCase):
    def exercise(self,change):
        with tempfile.TemporaryDirectory(prefix='r13-probe-control-')as tmp:
            p=Path(tmp);zero='0'*20
            row=f'leaf-00000\tf2xm1\t133f\t{zero}\t{zero}'
            (p/'request.tsv').write_text('id\tkind\tcw\ta\tb\n'+row+'\n')
            (p/'reference.tsv').write_text('id\tkind\tcw\ta\tb\tresult\n'+row+'\t'+zero+'\n')
            actual=row+'\t'+zero+'\n'
            if change=='result':actual=row+'\t'+'01'+zero[2:]+'\n'
            elif change=='input':actual=row.replace('133f','123f')+'\t'+zero+'\n'
            elif change=='duplicate':actual+=actual
            (p/'actual.tsv').write_text('id\tkind\tcw\ta\tb\tresult\n'+actual)
            return compare(p/'request.tsv',p/'reference.tsv',p/'actual.tsv',p/'receipt.json')

    def test_changed_native_output_remains_a_mismatch(self):
        result=self.exercise('result');self.assertEqual(result['mismatches'],1)
        self.assertFalse(result['all_certificate_leaves_reproduced'])

    def test_changed_native_operand_is_not_a_semantic_comparison(self):
        with self.assertRaisesRegex(ValueError,'changed native'):self.exercise('input')

    def test_duplicate_native_measurement_does_not_count_twice(self):
        with self.assertRaisesRegex(ValueError,'duplicate'):self.exercise('duplicate')


class StructuralVarianceControls(unittest.TestCase):
    def test_compensation_recovers_the_small_addend_before_averaging(self):
        from variance import evaluate
        for reduction, large in [(5, 2**53), (6, 2**53), (7, 2**54), (8, 2**54)]:
            profile = dict(algorithm=0, linear_order=reduction, square_order=0,
                           models=['strict']*5)
            result = evaluate(profile, [Q(large), Q(1), Q(-large)], 'AVERAGE')
            self.assertEqual(result, B.round(Q(1, 3)))

    def test_extra_pass_and_merge_forms_keep_translation_and_sample_divisor(self):
        from variance import evaluate, MODELS
        values = [Q(2**40+1), Q(2**40+3), Q(2**40+5)]
        for algorithm in range(14, 22):
            for mode in MODELS:
                for order in (0, 1):
                    profile = dict(algorithm=algorithm, linear_order=order,
                                   square_order=order, models=[mode]*5)
                    self.assertEqual(evaluate(profile, values, 'VAR.S'), Q(4),
                                     (algorithm, mode, order))


if __name__=='__main__':unittest.main()
